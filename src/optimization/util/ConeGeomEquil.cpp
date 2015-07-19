/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
// TODO: Consider making MinAbsNonzero global
#include "../../lapack_like/equilibrate/Util.hpp"

// The following are extensions to cones of the routine
// 'StackedGeomEquil'

namespace El {

template<typename Real>
inline Real DampScaling( Real alpha )
{
    const Real tol = Pow(Epsilon<Real>(),Real(0.33));
    if( alpha == Real(0) )
        return 1;
    else
        return Max(alpha,tol);
}

template<typename F>
void ConeGeomEquil
(       Matrix<F>& A, 
        Matrix<F>& B, 
        Matrix<Base<F>>& dRowA, 
        Matrix<Base<F>>& dRowB, 
        Matrix<Base<F>>& dCol, 
  const Matrix<Int>& orders,  
  const Matrix<Int>& firstInds,
  bool progress )
{
    DEBUG_ONLY(CSE cse("ConeGeomEquil"))
    LogicError("This routine is not yet written");
}

template<typename F>
void ConeGeomEquil
(       AbstractDistMatrix<F>& APre, 
        AbstractDistMatrix<F>& BPre,
        AbstractDistMatrix<Base<F>>& dRowAPre, 
        AbstractDistMatrix<Base<F>>& dRowBPre,
        AbstractDistMatrix<Base<F>>& dColPre,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff,
  bool progress )
{
    DEBUG_ONLY(CSE cse("ConeGeomEquil"))
    typedef Base<F> Real;

    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    auto APtr     = ReadWriteProxy<F,MC,MR>(&APre,control);
    auto BPtr     = ReadWriteProxy<F,MC,MR>(&BPre,control);
    auto dRowAPtr = WriteProxy<Real,MC,STAR>(&dRowAPre,control); 
    auto dRowBPtr = WriteProxy<Real,MC,STAR>(&dRowBPre,control); 
    auto dColPtr  = WriteProxy<Real,MR,STAR>(&dColPre,control);
    auto& A = *APtr;
    auto& B = *BPtr;
    auto& dRowA = *dRowAPtr;
    auto& dRowB = *dRowBPtr;
    auto& dCol  = *dColPtr;

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6; 
    const Real relTol = Real(9)/Real(10);

    // TODO: Incorporate damping
    //const Real damp = Real(1)/Real(1000);
    //const Real sqrtDamp = Sqrt(damp);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbsA = MaxAbs( A );
    auto maxAbsB = MaxAbs( B );
    const Real maxAbsVal = Max(maxAbsA.value,maxAbsB.value);
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsValA = MinAbsNonzero( A, maxAbsVal );
    const Real minAbsValB = MinAbsNonzero( B, maxAbsVal );
    const Real minAbsVal = Min(minAbsValA,minAbsValB);
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && A.Grid().Rank() == 0 )
        Output("Original ratio is ",maxAbsVal,"/",minAbsVal,"=",ratio);

    DistMatrix<Real,MC,STAR> rowScaleA(A.Grid()),
                             rowScaleB(A.Grid());
    DistMatrix<Real,MR,STAR> colScale(A.Grid()), colScaleB(B.Grid());
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        StackedGeometricColumnScaling( A, B, colScale ); 
        EntrywiseMap( colScale, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically equilibrate the rows
        // ----------------------------------
        GeometricRowScaling( A, rowScaleA );
        EntrywiseMap( rowScaleA, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        // TODO: Get rid of GeometricRowScaling since we need intrusive change
        GeometricRowScaling( B, rowScaleB );
        ConeAllReduce( rowScaleB, orders, firstInds, mpi::MAX, cutoff );
        EntrywiseMap( rowScaleB, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleB, dRowB );
        DiagonalSolve( LEFT, NORMAL, rowScaleB, B );

        auto newMaxAbsA = MaxAbs( A );
        auto newMaxAbsB = MaxAbs( B );
        const Real newMaxAbsVal = Max(newMaxAbsA.value,newMaxAbsB.value);
        const Real newMinAbsValA = MinAbsNonzero( A, newMaxAbsVal );
        const Real newMinAbsValB = MinAbsNonzero( B, newMaxAbsVal );
        const Real newMinAbsVal = Min(newMinAbsValA,newMinAbsValB);
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && A.Grid().Rank() == 0 )
            Output("New ratio is ",newMaxAbsVal,"/",newMinAbsVal,"=",newRatio);
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }
    SetIndent( indent );

    // Scale each column so that its maximum entry is 1 or 0
    // =====================================================
    colScaleB.AlignWith( colScale );
    ColumnMaxNorms( A, colScale );
    ColumnMaxNorms( B, colScaleB );
    const Int nLocal = A.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real maxScale = 
            Max(colScale.GetLocal(jLoc,0),colScaleB.GetLocal(jLoc,0));
        maxScale = DampScaling<Real>(maxScale); 
        colScale.SetLocal(jLoc,0,maxScale);
    }
    DiagonalScale( LEFT, NORMAL, colScale, dCol );
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
    DiagonalSolve( RIGHT, NORMAL, colScale, B );
}

template<typename F>
void ConeGeomEquil
(       SparseMatrix<F>& A, 
        SparseMatrix<F>& B,
        Matrix<Base<F>>& dRowA, 
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  bool progress )
{
    DEBUG_ONLY(CSE cse("ConeGeomEquil"))
    typedef Base<F> Real;
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbsA = MaxAbs( A );
    auto maxAbsB = MaxAbs( B );
    const Real maxAbsVal = Max(maxAbsA.value,maxAbsB.value);
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsValA = MinAbsNonzero( A, maxAbsVal );
    const Real minAbsValB = MinAbsNonzero( B, maxAbsVal );
    const Real minAbsVal = Min(minAbsValA,minAbsValB);
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        Output("Original ratio is ",maxAbsVal,"/",minAbsVal,"=",ratio);

    const Real sqrtDamp = Sqrt(damp);
    Matrix<Real> rowScaleA(mA,1), rowScaleB(mB,1), colScale(n,1),
                 maxAbsValsA, maxAbsValsB, minAbsValsA, minAbsValsB;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        ColumnMaxNorms( A, maxAbsValsA );
        ColumnMaxNorms( B, maxAbsValsB );
        for( Int j=0; j<n; ++j )
            maxAbsValsA.Set
            ( j, 0, Max(maxAbsValsA.Get(j,0),maxAbsValsB.Get(j,0)) );
        ColumnMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        ColumnMinAbsNonzero( B, maxAbsValsA, minAbsValsB );
        for( Int j=0; j<n; ++j )
            minAbsValsA.Set
            ( j, 0, Min(minAbsValsA.Get(j,0),minAbsValsB.Get(j,0)) );
        for( Int j=0; j<n; ++j )
        {
            const Real maxAbs = maxAbsValsA.Get(j,0);
            const Real minAbs = minAbsValsA.Get(j,0);
            const Real propScale = Sqrt(minAbs*maxAbs);
            const Real scale = Max(propScale,sqrtDamp*maxAbs);
            colScale.Set(j,0,scale);
        }
        EntrywiseMap( colScale, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically rescale the rows
        // ------------------------------ 
        RowMaxNorms( A, maxAbsValsA );
        RowMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        for( Int i=0; i<mA; ++i )
        {
            const Real maxAbs = maxAbsValsA.Get(i,0);
            const Real minAbs = minAbsValsA.Get(i,0);
            const Real propScale = Sqrt(minAbs*maxAbs);
            const Real scale = Max(propScale,sqrtDamp*maxAbs);
            rowScaleA.Set( i, 0, scale );
        }
        EntrywiseMap( rowScaleA, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        RowMaxNorms( B, maxAbsValsB );
        ConeAllReduce( maxAbsValsB, orders, firstInds, mpi::MAX );
        ColumnMinAbsNonzero( B, maxAbsValsB, minAbsValsB );
        ConeAllReduce( minAbsValsB, orders, firstInds, mpi::MIN );
        for( Int i=0; i<mB; ++i )
        {
            const Real maxAbs = maxAbsValsB.Get(i,0);
            const Real minAbs = minAbsValsB.Get(i,0);
            const Real propScale = Sqrt(minAbs*maxAbs);
            const Real scale = Max(propScale,sqrtDamp*maxAbs);
            rowScaleB.Set( i, 0, scale );
        }
        EntrywiseMap( rowScaleB, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleB, dRowB );
        DiagonalSolve( LEFT, NORMAL, rowScaleB, B );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbsA = MaxAbs( A );
        auto newMaxAbsB = MaxAbs( B );
        const Real newMaxAbsVal = Max(newMaxAbsA.value,newMaxAbsB.value);
        const Real newMinAbsValA = MinAbsNonzero( A, newMaxAbsVal );
        const Real newMinAbsValB = MinAbsNonzero( B, newMaxAbsVal );
        const Real newMinAbsVal = Min(newMinAbsValA,newMinAbsValB);
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            Output("New ratio is ",newMaxAbsVal,"/",newMinAbsVal,"=",newRatio);
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }
    SetIndent( indent );

    // Scale each row so that its maximum entry is 1 or 0
    F* valBufA = A.ValueBuffer();
    for( Int i=0; i<mA; ++i )
    {
        const Int offset = A.RowOffset(i);
        const Int numConnect = A.NumConnections(i);

        // Compute the maximum value in this row
        Real maxRowAbs = 0;
        for( Int e=offset; e<offset+numConnect; ++e )
            maxRowAbs = Max(maxRowAbs,Abs(A.Value(e)));

        if( maxRowAbs > Real(0) )
        {
            dRowA.Set(i,0, maxRowAbs*dRowA.Get(i,0) );
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufA[e] /= maxRowAbs;
        }
    }
    // Scale the rows of B such that the maximum entry in each cone is 1 or 0
    RowMaxNorms( B, maxAbsValsB );
    ConeAllReduce( maxAbsValsB, orders, firstInds, mpi::MAX );
    for( Int i=0; i<mB; ++i )
    {
        const Real maxAbs = maxAbsValsB.Get(i,0);
        if( maxAbs == Real(0) )
            maxAbsValsB.Set(i,0,1);
    }
    DiagonalScale( LEFT, NORMAL, maxAbsValsB, dRowB );
    DiagonalSolve( LEFT, NORMAL, maxAbsValsB, B );
}

template<typename F>
void ConeGeomEquil
(       DistSparseMatrix<F>& A, 
        DistSparseMatrix<F>& B,
        DistMultiVec<Base<F>>& dRowA, 
        DistMultiVec<Base<F>>& dRowB, 
        DistMultiVec<Base<F>>& dCol, 
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff,
  bool progress )
{
    DEBUG_ONLY(CSE cse("ConeGeomEquil"))
    typedef Base<F> Real;
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);
    dRowA.SetComm( comm );
    dRowB.SetComm( comm );
    dCol.SetComm( comm );
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbsA = MaxAbs( A );
    auto maxAbsB = MaxAbs( B );
    const Real maxAbsVal = Max(maxAbsA.value,maxAbsB.value);
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsValA = MinAbsNonzero( A, maxAbsVal );
    const Real minAbsValB = MinAbsNonzero( B, maxAbsVal );
    const Real minAbsVal = Min(minAbsValA,minAbsValB);
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && commRank == 0 )
        Output("Original ratio is ",maxAbsVal,"/",minAbsVal,"=",ratio);

    const Real sqrtDamp = Sqrt(damp);
    DistMultiVec<Real> maxAbsValsA(comm), maxAbsValsB(comm),
                       minAbsValsA(comm), minAbsValsB(comm), scales(comm);
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        scales.Resize( n, 1 );
        ColumnMaxNorms( A, maxAbsValsA );
        ColumnMaxNorms( B, maxAbsValsB );
        const Int localWidth = maxAbsValsA.LocalHeight();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            maxAbsValsA.SetLocal
            ( jLoc, 0, Max(maxAbsValsA.GetLocal(jLoc,0),
                           maxAbsValsB.GetLocal(jLoc,0)) );
        ColumnMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        ColumnMinAbsNonzero( B, maxAbsValsA, minAbsValsB );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            minAbsValsA.SetLocal
            ( jLoc, 0, Min(minAbsValsA.GetLocal(jLoc,0),
                           minAbsValsB.GetLocal(jLoc,0)) );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Real maxAbs = maxAbsValsA.GetLocal(jLoc,0);
            const Real minAbs = minAbsValsA.GetLocal(jLoc,0);
            const Real propScale = Sqrt(minAbs*maxAbs);
            const Real scale = Max(propScale,sqrtDamp*maxAbs);
            scales.SetLocal( jLoc, 0, scale );
        }
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dCol );
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( RIGHT, NORMAL, scales, B );

        // Geometrically rescale the rows
        // ------------------------------ 
        scales.Resize( mA, 1 );
        RowMaxNorms( A, maxAbsValsA );
        ColumnMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        const Int localHeightA = maxAbsValsA.LocalHeight();
        for( Int iLoc=0; iLoc<localHeightA; ++iLoc )
        {
            const Real maxAbs = maxAbsValsA.GetLocal(iLoc,0);
            const Real minAbs = minAbsValsA.GetLocal(iLoc,0);
            const Real propScale = Sqrt(minAbs*maxAbs);
            const Real scale = Max(propScale,sqrtDamp*maxAbs);
            scales.SetLocal( iLoc, 0, scale );
        }
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowA );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        scales.Resize( mB, 1 );
        RowMaxNorms( B, maxAbsValsB );
        ConeAllReduce( maxAbsValsB, orders, firstInds, mpi::MAX, cutoff );
        ColumnMinAbsNonzero( B, maxAbsValsB, minAbsValsB );
        ConeAllReduce( minAbsValsB, orders, firstInds, mpi::MIN, cutoff );
        const Int localHeightB = maxAbsValsB.LocalHeight();
        for( Int iLoc=0; iLoc<localHeightB; ++iLoc )
        {
            const Real maxAbs = maxAbsValsB.GetLocal(iLoc,0);
            const Real minAbs = minAbsValsB.GetLocal(iLoc,0);
            const Real propScale = Sqrt(minAbs*maxAbs);
            const Real scale = Max(propScale,sqrtDamp*maxAbs);
            scales.SetLocal( iLoc, 0, scale );
        }
        EntrywiseMap( scales, function<Real(Real)>(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowB );
        DiagonalSolve( LEFT, NORMAL, scales, B );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbsA = MaxAbs( A );
        auto newMaxAbsB = MaxAbs( B );
        const Real newMaxAbsVal = Max(newMaxAbsA.value,newMaxAbsB.value);
        const Real newMinAbsValA = MinAbsNonzero( A, newMaxAbsVal );
        const Real newMinAbsValB = MinAbsNonzero( B, newMaxAbsVal );
        const Real newMinAbsVal = Min(newMinAbsValA,newMinAbsValB);
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && commRank == 0 )
            Output("New ratio is ",newMaxAbsVal,"/",newMinAbsVal,"=",newRatio);
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }
    SetIndent( indent );

    // Scale each row of A so that its maximum entry is 1 or 0
    F* valBufA = A.ValueBuffer();
    const Int localHeightA = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeightA; ++iLoc )
    {
        const Int offset = A.RowOffset(iLoc);
        const Int numConnect = A.NumConnections(iLoc);

        // Compute the maximum value in this row
        Real maxRowAbs = 0;
        for( Int e=offset; e<offset+numConnect; ++e )
            maxRowAbs = Max(maxRowAbs,Abs(A.Value(e)));

        if( maxRowAbs > Real(0) )
        {
            dRowA.SetLocal(iLoc,0, maxRowAbs*dRowA.GetLocal(iLoc,0) );
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufA[e] /= maxRowAbs;
        }
    }
    // Scale the rows of B such that the maximum entry in each cone is 1 or 0
    RowMaxNorms( B, maxAbsValsB );
    ConeAllReduce( maxAbsValsB, orders, firstInds, mpi::MAX, cutoff );
    const Int localHeightB = B.LocalHeight();
    for( Int iLoc=0; iLoc<localHeightB; ++iLoc )
    {
        const Real maxAbs = maxAbsValsB.GetLocal(iLoc,0);
        if( maxAbs == Real(0) )
            maxAbsValsB.SetLocal(iLoc,0,1);
    }
    DiagonalScale( LEFT, NORMAL, maxAbsValsB, dRowB );
    DiagonalSolve( LEFT, NORMAL, maxAbsValsB, B );
}

#define PROTO(F) \
  template void ConeGeomEquil \
  (       Matrix<F>& A, \
          Matrix<F>& B, \
          Matrix<Base<F>>& dRowA, \
          Matrix<Base<F>>& dRowB, \
          Matrix<Base<F>>& dCol, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    bool progress ); \
  template void ConeGeomEquil \
  (       AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<F>& B, \
          AbstractDistMatrix<Base<F>>& dRowA, \
          AbstractDistMatrix<Base<F>>& dRowB, \
          AbstractDistMatrix<Base<F>>& dCol, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Int cutoff, bool progress ); \
  template void ConeGeomEquil \
  (       SparseMatrix<F>& A, \
          SparseMatrix<F>& B, \
          Matrix<Base<F>>& dRowA, \
          Matrix<Base<F>>& dRowB, \
          Matrix<Base<F>>& dCol, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    bool progress ); \
  template void ConeGeomEquil \
  (       DistSparseMatrix<F>& A, \
          DistSparseMatrix<F>& B, \
          DistMultiVec<Base<F>>& dRowA, \
          DistMultiVec<Base<F>>& dRowB, \
          DistMultiVec<Base<F>>& dCol, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff, bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
