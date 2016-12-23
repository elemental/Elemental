/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
// TODO(poulson): Consider making MinAbsNonzero global
#include "../../../lapack_like/equilibrate/Util.hpp"

// The following are extensions to cones of the routine
// 'StackedGeomEquil'

namespace El {

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Real DampScaling( const Real& alpha )
{
    const Real tol = Pow(limits::Epsilon<Real>(),Real(0.33));
    if( alpha == Real(0) )
        return 1;
    else
        return Max(alpha,tol);
}

namespace cone {

template<typename Field>
void GeomEquil
(       Matrix<Field>& A,
        Matrix<Field>& B,
        Matrix<Base<Field>>& dRowA,
        Matrix<Base<Field>>& dRowB,
        Matrix<Base<Field>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6;
    const Real relTol = Real(9)/Real(10);

    // TODO(poulson): Incorporate damping
    //const Real damp = Real(1)/Real(1000);
    //const Real sqrtDamp = Sqrt(damp);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbsA = MaxAbsLoc( A );
    auto maxAbsB = MaxAbsLoc( B );
    const Real maxAbsVal = Max(maxAbsA.value,maxAbsB.value);
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsValA = MinAbsNonzero( A, maxAbsVal );
    const Real minAbsValB = MinAbsNonzero( B, maxAbsVal );
    const Real minAbsVal = Min(minAbsValA,minAbsValB);
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        Output("Original ratio is ",maxAbsVal,"/",minAbsVal,"=",ratio);

    Matrix<Real> rowScaleA, rowScaleB, colScale, colScaleB;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        StackedGeometricColumnScaling( A, B, colScale );
        EntrywiseMap( colScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically equilibrate the rows
        // ----------------------------------
        GeometricRowScaling( A, rowScaleA );
        EntrywiseMap( rowScaleA, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        // TODO(poulson): Get rid of GeometricRowScaling since we need
        // intrusive change
        GeometricRowScaling( B, rowScaleB );
        cone::AllReduce( rowScaleB, orders, firstInds, mpi::MAX );
        EntrywiseMap( rowScaleB, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleB, dRowB );
        DiagonalSolve( LEFT, NORMAL, rowScaleB, B );

        auto newMaxAbsA = MaxAbsLoc( A );
        auto newMaxAbsB = MaxAbsLoc( B );
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

    // Scale each column so that its maximum entry is 1 or 0
    // =====================================================
    ColumnMaxNorms( A, colScale );
    ColumnMaxNorms( B, colScaleB );
          Real* colScaleBuf = colScale.Buffer();
    const Real* colScaleBBuf = colScaleB.LockedBuffer();
    for( Int j=0; j<n; ++j )
    {
        Real maxScale = Max(colScaleBuf[j],colScaleBBuf[j]);
        colScaleBuf[j] = DampScaling<Real>(maxScale);
    }
    DiagonalScale( LEFT, NORMAL, colScale, dCol );
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
    DiagonalSolve( RIGHT, NORMAL, colScale, B );
}

// TODO(poulson): Use lower-level access
template<typename Field>
void GeomEquil
(       AbstractDistMatrix<Field>& APre,
        AbstractDistMatrix<Field>& BPre,
        AbstractDistMatrix<Base<Field>>& dRowAPre,
        AbstractDistMatrix<Base<Field>>& dRowBPre,
        AbstractDistMatrix<Base<Field>>& dColPre,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<Field,Field,MC,MR>
      AProx( APre, control ),
      BProx( BPre, control );
    DistMatrixWriteProxy<Real,Real,MC,STAR>
      dRowAProx( dRowAPre, control ),
      dRowBProx( dRowBPre, control );
    DistMatrixWriteProxy<Real,Real,MR,STAR>
      dColProx( dColPre, control );
    auto& A = AProx.Get();
    auto& B = BProx.Get();
    auto& dRowA = dRowAProx.Get();
    auto& dRowB = dRowBProx.Get();
    auto& dCol = dColProx.Get();

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6;
    const Real relTol = Real(9)/Real(10);

    // TODO(poulson): Incorporate damping
    //const Real damp = Real(1)/Real(1000);
    //const Real sqrtDamp = Sqrt(damp);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbsA = MaxAbsLoc( A );
    auto maxAbsB = MaxAbsLoc( B );
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
        EntrywiseMap( colScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically equilibrate the rows
        // ----------------------------------
        GeometricRowScaling( A, rowScaleA );
        EntrywiseMap( rowScaleA, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        // TODO(poulson): Get rid of GeometricRowScaling since we need
        // intrusive change
        GeometricRowScaling( B, rowScaleB );
        cone::AllReduce( rowScaleB, orders, firstInds, mpi::MAX, cutoff );
        EntrywiseMap( rowScaleB, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleB, dRowB );
        DiagonalSolve( LEFT, NORMAL, rowScaleB, B );

        auto newMaxAbsA = MaxAbsLoc( A );
        auto newMaxAbsB = MaxAbsLoc( B );
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
          Real* colScaleBuf = colScale.Buffer();
    const Real* colScaleBBuf = colScaleB.LockedBuffer();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real maxScale = Max(colScaleBuf[jLoc],colScaleBBuf[jLoc]);
        colScaleBuf[jLoc] = DampScaling<Real>(maxScale);
    }
    DiagonalScale( LEFT, NORMAL, colScale, dCol );
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
    DiagonalSolve( RIGHT, NORMAL, colScale, B );
}

template<typename Field>
void GeomEquil
(       SparseMatrix<Field>& A,
        SparseMatrix<Field>& B,
        Matrix<Base<Field>>& dRowA,
        Matrix<Base<Field>>& dRowB,
        Matrix<Base<Field>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6;
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbsA = MaxAbsLoc( A );
    auto maxAbsB = MaxAbsLoc( B );
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
    Real* colScaleBuf = colScale.Buffer();
    Real* rowScaleABuf = rowScaleA.Buffer();
    Real* rowScaleBBuf = rowScaleB.Buffer();
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        ColumnMaxNorms( A, maxAbsValsA );
        ColumnMaxNorms( B, maxAbsValsB );
        {
                  Real* maxAbsValsABuf = maxAbsValsA.Buffer();
            const Real* maxAbsValsBBuf = maxAbsValsB.LockedBuffer();
            for( Int j=0; j<n; ++j )
                maxAbsValsABuf[j] = Max(maxAbsValsABuf[j],maxAbsValsBBuf[j]);
        }

        ColumnMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        ColumnMinAbsNonzero( B, maxAbsValsA, minAbsValsB );
        {
                  Real* minAbsValsABuf = minAbsValsA.Buffer();
            const Real* minAbsValsBBuf = minAbsValsB.LockedBuffer();
            for( Int j=0; j<n; ++j )
                minAbsValsABuf[j] = Min(minAbsValsABuf[j],minAbsValsBBuf[j]);
        }

        {
            const Real* maxAbsValsABuf = maxAbsValsA.Buffer();
            const Real* minAbsValsABuf = minAbsValsA.Buffer();
            for( Int j=0; j<n; ++j )
            {
                const Real maxAbs = maxAbsValsABuf[j];
                const Real minAbs = minAbsValsABuf[j];
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                colScaleBuf[j] = scale;
            }
        }
        EntrywiseMap( colScale, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically rescale the rows
        // ------------------------------
        RowMaxNorms( A, maxAbsValsA );
        RowMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        {
            const Real* maxAbsValsABuf = maxAbsValsA.LockedBuffer();
            const Real* minAbsValsABuf = minAbsValsA.LockedBuffer();
            for( Int i=0; i<mA; ++i )
            {
                const Real maxAbs = maxAbsValsABuf[i];
                const Real minAbs = minAbsValsABuf[i];
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                rowScaleABuf[i] = scale;
            }
        }
        EntrywiseMap( rowScaleA, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        RowMaxNorms( B, maxAbsValsB );
        cone::AllReduce( maxAbsValsB, orders, firstInds, mpi::MAX );
        ColumnMinAbsNonzero( B, maxAbsValsB, minAbsValsB );
        cone::AllReduce( minAbsValsB, orders, firstInds, mpi::MIN );
        {
            const Real* maxAbsValsBBuf = maxAbsValsB.LockedBuffer();
            const Real* minAbsValsBBuf = minAbsValsB.LockedBuffer();
            for( Int i=0; i<mB; ++i )
            {
                const Real maxAbs = maxAbsValsBBuf[i];
                const Real minAbs = minAbsValsBBuf[i];
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                rowScaleBBuf[i] = scale;
            }
        }
        EntrywiseMap( rowScaleB, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, rowScaleB, dRowB );
        DiagonalSolve( LEFT, NORMAL, rowScaleB, B );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbsA = MaxAbsLoc( A );
        auto newMaxAbsB = MaxAbsLoc( B );
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
    Field* valBufA = A.ValueBuffer();
    Real* dRowABuf = dRowA.Buffer();
    for( Int i=0; i<mA; ++i )
    {
        const Int offset = A.RowOffset(i);
        const Int numConnect = A.NumConnections(i);

        // Compute the maximum value in this row
        Real maxRowAbs = 0;
        for( Int e=offset; e<offset+numConnect; ++e )
            maxRowAbs = Max(maxRowAbs,Abs(valBufA[e]));

        if( maxRowAbs > Real(0) )
        {
            dRowABuf[i] = maxRowAbs*dRowABuf[i];
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufA[e] /= maxRowAbs;
        }
    }
    // Scale the rows of B such that the maximum entry in each cone is 1 or 0
    RowMaxNorms( B, maxAbsValsB );
    cone::AllReduce( maxAbsValsB, orders, firstInds, mpi::MAX );
    {
        Real* maxAbsValsBBuf = maxAbsValsB.Buffer();
        for( Int i=0; i<mB; ++i )
        {
            const Real maxAbs = maxAbsValsBBuf[i];
            if( maxAbs == Real(0) )
                maxAbsValsBBuf[i] = 1;
        }
    }
    DiagonalScale( LEFT, NORMAL, maxAbsValsB, dRowB );
    DiagonalSolve( LEFT, NORMAL, maxAbsValsB, B );
}

template<typename Field>
void GeomEquil
(       DistSparseMatrix<Field>& A,
        DistSparseMatrix<Field>& B,
        DistMultiVec<Base<Field>>& dRowA,
        DistMultiVec<Base<Field>>& dRowB,
        DistMultiVec<Base<Field>>& dCol,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const int commRank = grid.Rank();

    dRowA.SetGrid( grid );
    dRowB.SetGrid( grid );
    dCol.SetGrid( grid );
    Ones( dRowA, mA, 1 );
    Ones( dRowB, mB, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6;
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbsA = MaxAbsLoc( A );
    auto maxAbsB = MaxAbsLoc( B );
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
    DistMultiVec<Real> maxAbsValsA(grid), maxAbsValsB(grid),
                       minAbsValsA(grid), minAbsValsB(grid), scales(grid);
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        scales.Resize( n, 1 );
        ColumnMaxNorms( A, maxAbsValsA );
        ColumnMaxNorms( B, maxAbsValsB );
        const Int localWidth = maxAbsValsA.LocalHeight();
        {
            auto maxAbsValsABuf = maxAbsValsA.Matrix().Buffer();
            auto maxAbsValsBBuf = maxAbsValsB.LockedMatrix().LockedBuffer();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                maxAbsValsABuf[jLoc] =
                  Max(maxAbsValsABuf[jLoc],maxAbsValsBBuf[jLoc]);
        }
        ColumnMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        ColumnMinAbsNonzero( B, maxAbsValsA, minAbsValsB );
        {
            auto minAbsValsABuf = minAbsValsA.Matrix().Buffer();
            auto minAbsValsBBuf = minAbsValsB.LockedMatrix().LockedBuffer();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
                minAbsValsABuf[jLoc] =
                  Min(minAbsValsABuf[jLoc],minAbsValsBBuf[jLoc]);

            auto scalesBuf = scales.Matrix().Buffer();
            auto maxAbsValsABuf = maxAbsValsA.LockedMatrix().LockedBuffer();
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Real maxAbs = maxAbsValsABuf[jLoc];
                const Real minAbs = minAbsValsABuf[jLoc];
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scalesBuf[jLoc] = scale;
            }
        }
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dCol );
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( RIGHT, NORMAL, scales, B );

        // Geometrically rescale the rows
        // ------------------------------
        scales.Resize( mA, 1 );
        RowMaxNorms( A, maxAbsValsA );
        ColumnMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        const Int localHeightA = maxAbsValsA.LocalHeight();
        {
            auto minAbsValsABuf = minAbsValsA.LockedMatrix().LockedBuffer();
            auto maxAbsValsABuf = maxAbsValsA.LockedMatrix().LockedBuffer();
            auto scalesBuf = scales.Matrix().Buffer();
            for( Int iLoc=0; iLoc<localHeightA; ++iLoc )
            {
                const Real maxAbs = maxAbsValsABuf[iLoc];
                const Real minAbs = minAbsValsABuf[iLoc];
                const Real propScale = Sqrt(minAbs*maxAbs);
                scalesBuf[iLoc] = Max(propScale,sqrtDamp*maxAbs);
            }
        }
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowA );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        scales.Resize( mB, 1 );
        RowMaxNorms( B, maxAbsValsB );
        cone::AllReduce( maxAbsValsB, orders, firstInds, mpi::MAX, cutoff );
        ColumnMinAbsNonzero( B, maxAbsValsB, minAbsValsB );
        cone::AllReduce( minAbsValsB, orders, firstInds, mpi::MIN, cutoff );
        const Int localHeightB = maxAbsValsB.LocalHeight();
        {
            auto minAbsValsBBuf = minAbsValsB.LockedMatrix().LockedBuffer();
            auto maxAbsValsBBuf = maxAbsValsB.LockedMatrix().LockedBuffer();
            auto scalesBuf = scales.Matrix().Buffer();
            for( Int iLoc=0; iLoc<localHeightB; ++iLoc )
            {
                const Real maxAbs = maxAbsValsBBuf[iLoc];
                const Real minAbs = minAbsValsBBuf[iLoc];
                const Real propScale = Sqrt(minAbs*maxAbs);
                scalesBuf[iLoc] = Max(propScale,sqrtDamp*maxAbs);
            }
        }
        EntrywiseMap( scales, MakeFunction(DampScaling<Real>) );
        DiagonalScale( LEFT, NORMAL, scales, dRowB );
        DiagonalSolve( LEFT, NORMAL, scales, B );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbsA = MaxAbsLoc( A );
        auto newMaxAbsB = MaxAbsLoc( B );
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
    Field* valBufA = A.ValueBuffer();
    const Int localHeightA = A.LocalHeight();
    Real* dRowABuf = dRowA.Matrix().Buffer();
    for( Int iLoc=0; iLoc<localHeightA; ++iLoc )
    {
        const Int offset = A.RowOffset(iLoc);
        const Int numConnect = A.NumConnections(iLoc);

        // Compute the maximum value in this row
        Real maxRowAbs = 0;
        for( Int e=offset; e<offset+numConnect; ++e )
            maxRowAbs = Max(maxRowAbs,Abs(valBufA[e]));

        if( maxRowAbs > Real(0) )
        {
            dRowABuf[iLoc] *= maxRowAbs;
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufA[e] /= maxRowAbs;
        }
    }
    // Scale the rows of B such that the maximum entry in each cone is 1 or 0
    RowMaxNorms( B, maxAbsValsB );
    cone::AllReduce( maxAbsValsB, orders, firstInds, mpi::MAX, cutoff );
    const Int localHeightB = B.LocalHeight();
    auto maxAbsValsBBuf = maxAbsValsB.Matrix().Buffer();
    for( Int iLoc=0; iLoc<localHeightB; ++iLoc )
    {
        const Real maxAbs = maxAbsValsBBuf[iLoc];
        if( maxAbs == Real(0) )
            maxAbsValsBBuf[iLoc] = 1;
    }
    DiagonalScale( LEFT, NORMAL, maxAbsValsB, dRowB );
    DiagonalSolve( LEFT, NORMAL, maxAbsValsB, B );
}

#define PROTO(Field) \
  template void GeomEquil \
  (       Matrix<Field>& A, \
          Matrix<Field>& B, \
          Matrix<Base<Field>>& dRowA, \
          Matrix<Base<Field>>& dRowB, \
          Matrix<Base<Field>>& dCol, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    bool progress ); \
  template void GeomEquil \
  (       AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& B, \
          AbstractDistMatrix<Base<Field>>& dRowA, \
          AbstractDistMatrix<Base<Field>>& dRowB, \
          AbstractDistMatrix<Base<Field>>& dCol, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    Int cutoff, bool progress ); \
  template void GeomEquil \
  (       SparseMatrix<Field>& A, \
          SparseMatrix<Field>& B, \
          Matrix<Base<Field>>& dRowA, \
          Matrix<Base<Field>>& dRowB, \
          Matrix<Base<Field>>& dCol, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    bool progress ); \
  template void GeomEquil \
  (       DistSparseMatrix<Field>& A, \
          DistSparseMatrix<Field>& B, \
          DistMultiVec<Base<Field>>& dRowA, \
          DistMultiVec<Base<Field>>& dRowB, \
          DistMultiVec<Base<Field>>& dCol, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    Int cutoff, bool progress );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace cone
} // namespace El
