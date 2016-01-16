/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "./Util.hpp"

// The following routines are adaptations of the approach uses by 
// Saunders et al. (originally recommended by Joseph Fourer) for iteratively 
// rescaling the columns and rows by their approximate geometric means in order
// to better scale the original problem. After this iteration is finished,
// the columns or rows are rescaled so that their maximum entry has magnitude 
// one (unless the row/column is identically zero).
//
// The implementation of Saunders et al. is commonly referred to as either
// gmscale.m or gmscal.m. 

namespace El {

// TODO: Make this consistent with ConeGeomEquil

template<typename F>
void GeomEquil
( Matrix<F>& A,
  Matrix<Base<F>>& dRow,
  Matrix<Base<F>>& dCol,
  bool progress )
{
    DEBUG_ONLY(CSE cse("GeomEquil"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbsLoc( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        Output("Original ratio is ",maxAbsVal,"/",minAbsVal,"=",ratio);

    const Real sqrtDamp = Sqrt(damp);
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        for( Int j=0; j<n; ++j )
        {
            auto aCol = A( ALL, IR(j) );
            auto maxColAbs = VectorMaxAbsLoc( aCol );
            const Real maxColAbsVal = maxColAbs.value;
            if( maxColAbsVal > Real(0) )
            {
                const Real minColAbsVal = MinAbsNonzero( aCol, maxColAbsVal );
                const Real propScale = Sqrt(minColAbsVal*maxColAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxColAbsVal);
                aCol *= 1/scale;
                dCol.Set( j, 0, scale*dCol.Get(j,0) );
            }
        }

        // Geometrically equilibrate the rows
        for( Int i=0; i<m; ++i )
        {
            auto aRow = A( IR(i), ALL );
            auto maxRowAbs = VectorMaxAbsLoc( aRow );
            const Real maxRowAbsVal = maxRowAbs.value;
            if( maxRowAbsVal > Real(0) )
            {
                const Real minRowAbsVal = MinAbsNonzero( aRow, maxRowAbsVal );
                const Real propScale = Sqrt(minRowAbsVal*maxRowAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbsVal);
                aRow *= 1/scale;
                dRow.Set( i, 0, scale*dRow.Get(i,0) );
            }
        }

        auto newMaxAbs = MaxAbsLoc( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            Output("New ratio is ",newMaxAbsVal,"/",newMinAbsVal,"=",newRatio);
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }
    SetIndent( indent );

    // Scale each column so that its maximum entry is 1 or 0
    for( Int j=0; j<n; ++j )
    {
        auto aCol = A( ALL, IR(j) );
        auto maxColAbs = VectorMaxAbsLoc( aCol );
        const Real maxColAbsVal = maxColAbs.value;
        if( maxColAbsVal > Real(0) )
        {
            aCol *= 1/maxColAbsVal;
            dCol.Set( j, 0, maxColAbsVal*dCol.Get(j,0) );
        }
    }
}

template<typename F>
void StackedGeomEquil
( Matrix<F>& A,
  Matrix<F>& B, 
  Matrix<Base<F>>& dRowA,
  Matrix<Base<F>>& dRowB, 
  Matrix<Base<F>>& dCol,
  bool progress )
{
    DEBUG_ONLY(CSE cse("StackedGeomEquil"))
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
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        for( Int j=0; j<n; ++j )
        {
            auto aCol = A( ALL, IR(j) );
            auto bCol = B( ALL, IR(j) );
            auto maxColAbsA = VectorMaxAbsLoc( aCol );
            auto maxColAbsB = VectorMaxAbsLoc( bCol );
            const Real maxColAbsVal = Max(maxColAbsA.value,maxColAbsB.value);
            if( maxColAbsVal > Real(0) )
            {
                const Real minColAbsAVal = MinAbsNonzero( aCol, maxColAbsVal );
                const Real minColAbsBVal = MinAbsNonzero( bCol, maxColAbsVal );
                const Real minColAbsVal = Min(minColAbsAVal,minColAbsBVal);
                const Real propScale = Sqrt(minColAbsVal*maxColAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxColAbsVal);
                aCol *= 1/scale;
                bCol *= 1/scale;
                dCol.Set( j, 0, scale*dCol.Get(j,0) );
            }
        }

        // Geometrically equilibrate the rows
        for( Int i=0; i<mA; ++i )
        {
            auto aRow = A( IR(i), ALL );
            auto maxRowAbs = VectorMaxAbsLoc( aRow );
            const Real maxRowAbsVal = maxRowAbs.value;
            if( maxRowAbsVal > Real(0) )
            {
                const Real minRowAbsVal = MinAbsNonzero( aRow, maxRowAbsVal );
                const Real propScale = Sqrt(minRowAbsVal*maxRowAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbsVal);
                aRow *= 1/scale;
                dRowA.Set( i, 0, scale*dRowA.Get(i,0) );
            }
        }
        for( Int i=0; i<mB; ++i )
        {
            auto bRow = B( IR(i), ALL );
            auto maxRowAbs = VectorMaxAbsLoc( bRow );
            const Real maxRowAbsVal = maxRowAbs.value;
            if( maxRowAbsVal > Real(0) )
            {
                const Real minRowAbsVal = MinAbsNonzero( bRow, maxRowAbsVal );
                const Real propScale = Sqrt(minRowAbsVal*maxRowAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbsVal);
                bRow *= 1/scale;
                dRowB.Set( i, 0, scale*dRowB.Get(i,0) );
            }
        }

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
    for( Int j=0; j<n; ++j )
    {
        auto aCol = A( ALL, IR(j) );
        auto bCol = B( ALL, IR(j) );
        auto maxColAbsA = VectorMaxAbsLoc( aCol );
        auto maxColAbsB = VectorMaxAbsLoc( bCol );
        const Real maxColAbsVal = Max(maxColAbsA.value,maxColAbsB.value);
        if( maxColAbsVal > Real(0) )
        {
            aCol *= 1/maxColAbsVal;
            bCol *= 1/maxColAbsVal;
            dCol.Set( j, 0, maxColAbsVal*dCol.Get(j,0) );
        }
    }
}

template<typename F>
void GeomEquil
( ElementalMatrix<F>& APre, 
  ElementalMatrix<Base<F>>& dRowPre,
  ElementalMatrix<Base<F>>& dColPre,
  bool progress )
{
    DEBUG_ONLY(CSE cse("GeomEquil"))
    typedef Base<F> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre, control );
    DistMatrixWriteProxy<Real,Real,MC,STAR> dRowProx( dRowPre, control );
    DistMatrixWriteProxy<Real,Real,MR,STAR> dColProx( dColPre, control );
    auto& A = AProx.Get();
    auto& dRow = dRowProx.Get();
    auto& dCol = dColProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6; 
    const Real relTol = Real(9)/Real(10);

    // TODO: Incorporate damping
    //const Real damp = Real(1)/Real(1000);
    //const Real sqrtDamp = Sqrt(damp);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbsLoc( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && A.Grid().Rank() == 0 )
        Output("Original ratio is ",maxAbsVal,"/",minAbsVal,"=",ratio);

    DistMatrix<Real,MC,STAR> rowScale(A.Grid());
    DistMatrix<Real,MR,STAR> colScale(A.Grid());
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        // TODO: Remove GeometricColumnScaling
        GeometricColumnScaling( A, colScale ); 
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            if( colScale.GetLocal(jLoc,0) == Real(0) )
                colScale.SetLocal(jLoc,0,Real(1));
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );

        // Geometrically equilibrate the rows
        // ----------------------------------
        // TODO: Remove GeometricRowScaling
        GeometricRowScaling( A, rowScale );
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            if( rowScale.GetLocal(iLoc,0) == Real(0) )
                rowScale.SetLocal(iLoc,0,Real(1));
        DiagonalScale( LEFT, NORMAL, rowScale, dRow );
        DiagonalSolve( LEFT, NORMAL, rowScale, A );

        auto newMaxAbs = MaxAbsLoc( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && A.Grid().Rank() == 0 )
            Output("New ratio is ",newMaxAbsVal,"/",newMinAbsVal,"=",newRatio);
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }
    SetIndent( indent );

    // Scale each column so that its maximum entry is 1 or 0
    ColumnMaxNorms( A, colScale );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        if( colScale.GetLocal(jLoc,0) == Real(0) ) 
            colScale.SetLocal(jLoc,0,Real(1));
    DiagonalScale( LEFT, NORMAL, colScale, dCol );
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
}

template<typename F>
void StackedGeomEquil
( ElementalMatrix<F>& APre, 
  ElementalMatrix<F>& BPre,
  ElementalMatrix<Base<F>>& dRowAPre, 
  ElementalMatrix<Base<F>>& dRowBPre,
  ElementalMatrix<Base<F>>& dColPre,
  bool progress )
{
    DEBUG_ONLY(CSE cse("StackedGeomEquil"))
    typedef Base<F> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre, control );
    DistMatrixReadWriteProxy<F,F,MC,MR> BProx( BPre, control );
    DistMatrixWriteProxy<Real,Real,MC,STAR> dRowAProx( dRowAPre, control );
    DistMatrixWriteProxy<Real,Real,MC,STAR> dRowBProx( dRowBPre, control );
    DistMatrixWriteProxy<Real,Real,MR,STAR> dColProx( dColPre, control );
    auto& A = AProx.Get();
    auto& B = BProx.Get();
    auto& dRowA = dRowAProx.Get();
    auto& dRowB = dRowBProx.Get();
    auto& dCol = dColProx.Get();

    const Int mA = A.Height();
    const Int mB = B.Height();
    const Int n = A.Width();
    const Int mLocalA = A.LocalHeight();
    const Int mLocalB = B.LocalHeight();
    const Int nLocal = A.LocalWidth();
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
        // TODO: Remove StackedGeometricColumnScaling
        StackedGeometricColumnScaling( A, B, colScale ); 
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            if( colScale.GetLocal(jLoc,0) == Real(0) )
                colScale.SetLocal(jLoc,0,Real(1));
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically equilibrate the rows
        // ----------------------------------
        // TODO: Remove GeometricRowScaling
        GeometricRowScaling( A, rowScaleA );
        for( Int iLoc=0; iLoc<mLocalA; ++iLoc )
            if( rowScaleA.GetLocal(iLoc,0) == Real(0) )
                rowScaleA.SetLocal(iLoc,0,Real(1));
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        // TODO: Remove GeometricRowScaling
        GeometricRowScaling( B, rowScaleB );
        for( Int iLoc=0; iLoc<mLocalB; ++iLoc )
            if( rowScaleB.GetLocal(iLoc,0) == Real(0) )
                rowScaleB.SetLocal(iLoc,0,Real(1));
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
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real maxScale = 
            Max(colScale.GetLocal(jLoc,0),colScaleB.GetLocal(jLoc,0));
        if( maxScale == Real(0) )
            maxScale = 1; 
        colScale.SetLocal(jLoc,0,maxScale);
    }
    DiagonalScale( LEFT, NORMAL, colScale, dCol );
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
    DiagonalSolve( RIGHT, NORMAL, colScale, B );
}

template<typename F>
void GeomEquil
( SparseMatrix<F>& A,
  Matrix<Base<F>>& dRow,
  Matrix<Base<F>>& dCol,
  bool progress )
{
    DEBUG_ONLY(CSE cse("GeomEquil"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbsLoc( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        Output("Original ratio is ",maxAbsVal,"/",minAbsVal,"=",ratio);

    const Real sqrtDamp = Sqrt(damp);
    Matrix<Real> rowScale(m,1), colScale(n,1), maxAbsVals, minAbsVals;
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        ColumnMaxNorms( A, maxAbsVals );
        ColumnMinAbsNonzero( A, maxAbsVals, minAbsVals );
        for( Int j=0; j<n; ++j )
        {
            const Real maxAbs = maxAbsVals.Get(j,0);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsVals.Get(j,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                colScale.Set(j,0,scale);
            }
            else
            {
                colScale.Set(j,0,1);
            }
        }
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );

        // Geometrically rescale the rows
        // ------------------------------ 
        RowMaxNorms( A, maxAbsVals );
        RowMinAbsNonzero( A, maxAbsVals, minAbsVals );
        for( Int i=0; i<m; ++i )
        {
            const Real maxAbs = maxAbsVals.Get(i,0);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsVals.Get(i,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                rowScale.Set(i,0,scale);
            } 
            else
                rowScale.Set(i,0,Real(1));
        }
        DiagonalScale( LEFT, NORMAL, rowScale, dRow );
        DiagonalSolve( LEFT, NORMAL, rowScale, A );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbsLoc( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            Output("New ratio is ",newMaxAbsVal,"/",newMinAbsVal,"=",newRatio);
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }
    SetIndent( indent );

    // Scale each row so that its maximum entry is 1 or 0
    F* valBuf = A.ValueBuffer();
    for( Int i=0; i<m; ++i )
    {
        const Int offset = A.RowOffset(i);
        const Int numConnect = A.NumConnections(i);

        // Compute the maximum value in this row
        Real maxRowAbs = 0;
        for( Int e=offset; e<offset+numConnect; ++e )
            maxRowAbs = Max(maxRowAbs,Abs(A.Value(e)));

        if( maxRowAbs > Real(0) )
        {
            dRow.Set(i,0, maxRowAbs*dRow.Get(i,0) );
            for( Int e=offset; e<offset+numConnect; ++e )
                valBuf[e] /= maxRowAbs;
        }
    }
}

template<typename F>
void StackedGeomEquil
( SparseMatrix<F>& A,
  SparseMatrix<F>& B,
  Matrix<Base<F>>& dRowA,
  Matrix<Base<F>>& dRowB,
  Matrix<Base<F>>& dCol,
  bool progress )
{
    DEBUG_ONLY(CSE cse("StackedGeomEquil"))
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
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsA.Get(j,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                colScale.Set(j,0,scale);
            }
            else
                colScale.Set(j,0,Real(1));
        }
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
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsA.Get(i,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                rowScaleA.Set(i,0,scale);
            } 
            else
                rowScaleA.Set(i,0,Real(1));
        }
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        RowMinAbsNonzero( B, maxAbsValsB, minAbsValsB );
        RowMaxNorms( B, maxAbsValsB );
        for( Int i=0; i<mB; ++i )
        {
            const Real maxAbs = maxAbsValsB.Get(i,0);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsB.Get(i,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                rowScaleB.Set(i,0,scale);
            } 
            else
                rowScaleB.Set(i,0,Real(1));
        }
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
    F* valBufB = B.ValueBuffer(); 
    for( Int i=0; i<mB; ++i )
    {
        const Int offset = B.RowOffset(i);
        const Int numConnect = B.NumConnections(i);

        // Compute the maximum value in this row
        Real maxRowAbs = 0;
        for( Int e=offset; e<offset+numConnect; ++e )
            maxRowAbs = Max(maxRowAbs,Abs(B.Value(e)));

        if( maxRowAbs > Real(0) )
        {
            dRowB.Set(i,0, maxRowAbs*dRowB.Get(i,0) );
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufB[e] /= maxRowAbs;
        }
    }
}

template<typename F>
void GeomEquil
( DistSparseMatrix<F>& A, 
  DistMultiVec<Base<F>>& dRow,
  DistMultiVec<Base<F>>& dCol, 
  bool progress )
{
    DEBUG_ONLY(CSE cse("GeomEquil"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);
    dRow.SetComm( comm );
    dCol.SetComm( comm );
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbsLoc( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && commRank == 0 )
        Output("Original ratio is ",maxAbsVal,"/",minAbsVal,"=",ratio);

    const Real sqrtDamp = Sqrt(damp);
    DistMultiVec<Real> maxAbsVals(comm), minAbsVals(comm), scales(comm);
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        ColumnMaxNorms( A, maxAbsVals );
        ColumnMinAbsNonzero( A, maxAbsVals, minAbsVals );
        scales.Resize( n, 1 );
        const Int localWidth = maxAbsVals.LocalHeight(); 
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Real maxAbs = maxAbsVals.GetLocal(jLoc,0);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsVals.GetLocal(jLoc,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scales.SetLocal( jLoc, 0, scale );
            }
            else
                scales.SetLocal( jLoc, 0, 1 ); 
        }
        DiagonalScale( LEFT, NORMAL, scales, dCol );
        DiagonalSolve( RIGHT, NORMAL, scales, A );

        // Geometrically rescale the rows
        // ------------------------------ 
        RowMaxNorms( A, maxAbsVals );
        RowMinAbsNonzero( A, maxAbsVals, minAbsVals );
        scales.Resize( m, 1 ); 
        const Int localHeight = maxAbsVals.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Real maxAbs = maxAbsVals.GetLocal(iLoc,0);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsVals.GetLocal(iLoc,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scales.SetLocal( iLoc, 0, scale );
            }
            else
                scales.SetLocal( iLoc, 0, 1 );
        }
        DiagonalScale( LEFT, NORMAL, scales, dRow );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbsLoc( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && commRank == 0 )
            Output("New ratio is ",newMaxAbsVal,"/",newMinAbsVal,"=",newRatio);
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }
    SetIndent( indent );

    // Scale each row so that its maximum entry is 1 or 0
    F* valBuf = A.ValueBuffer();
    const Int localHeight = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int offset = A.RowOffset(iLoc);
        const Int numConnect = A.NumConnections(iLoc);

        // Compute the maximum value in this row
        Real maxRowAbs = 0;
        for( Int e=offset; e<offset+numConnect; ++e )
            maxRowAbs = Max(maxRowAbs,Abs(A.Value(e)));

        if( maxRowAbs > Real(0) )
        {
            dRow.SetLocal(iLoc,0, maxRowAbs*dRow.GetLocal(iLoc,0) );
            for( Int e=offset; e<offset+numConnect; ++e )
                valBuf[e] /= maxRowAbs;
        }
    }
}

template<typename F>
void StackedGeomEquil
( DistSparseMatrix<F>& A,
  DistSparseMatrix<F>& B,
  DistMultiVec<Base<F>>& dRowA, 
  DistMultiVec<Base<F>>& dRowB, 
  DistMultiVec<Base<F>>& dCol, 
  bool progress )
{
    DEBUG_ONLY(CSE cse("StackedGeomEquil"))
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
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsA.GetLocal(jLoc,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scales.SetLocal( jLoc, 0, scale );
            }
            else
                scales.SetLocal( jLoc, 0, 1 );
        }
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
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsA.GetLocal(iLoc,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scales.SetLocal( iLoc, 0, scale );
            }
            else
                scales.SetLocal( iLoc, 0, 1 );
        }
        DiagonalScale( LEFT, NORMAL, scales, dRowA );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        scales.Resize( mB, 1 );
        RowMaxNorms( B, maxAbsValsB );
        ColumnMinAbsNonzero( B, maxAbsValsB, minAbsValsB );
        const Int localHeightB = maxAbsValsB.LocalHeight();
        for( Int iLoc=0; iLoc<localHeightB; ++iLoc )
        {
            const Real maxAbs = maxAbsValsB.GetLocal(iLoc,0);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsB.GetLocal(iLoc,0);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scales.SetLocal( iLoc, 0, scale );
            }
            else
                scales.SetLocal( iLoc, 0, 1 );
        }
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

    // Scale each row so that its maximum entry is 1 or 0
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
    F* valBufB = B.ValueBuffer();
    const Int localHeightB = B.LocalHeight();
    for( Int iLoc=0; iLoc<localHeightB; ++iLoc )
    {
        const Int offset = B.RowOffset(iLoc);
        const Int numConnect = B.NumConnections(iLoc);

        // Compute the maximum value in this row
        Real maxRowAbs = 0;
        for( Int e=offset; e<offset+numConnect; ++e )
            maxRowAbs = Max(maxRowAbs,Abs(B.Value(e)));

        if( maxRowAbs > Real(0) )
        {
            dRowB.SetLocal(iLoc,0, maxRowAbs*dRowB.GetLocal(iLoc,0) );
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufB[e] /= maxRowAbs;
        }
    }
}

#define PROTO(F) \
  template void GeomEquil \
  ( Matrix<F>& A, \
    Matrix<Base<F>>& dRow, \
    Matrix<Base<F>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& dRow, \
    ElementalMatrix<Base<F>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( SparseMatrix<F>& A, \
    Matrix<Base<F>>& dRow, \
    Matrix<Base<F>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( DistSparseMatrix<F>& A, \
    DistMultiVec<Base<F>>& dRow, \
    DistMultiVec<Base<F>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( Matrix<F>& A, \
    Matrix<F>& B, \
    Matrix<Base<F>>& dRowA, \
    Matrix<Base<F>>& dRowB, \
    Matrix<Base<F>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& B, \
    ElementalMatrix<Base<F>>& dRowA, \
    ElementalMatrix<Base<F>>& dRowB, \
    ElementalMatrix<Base<F>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( SparseMatrix<F>& A, \
    SparseMatrix<F>& B, \
    Matrix<Base<F>>& dRowA, \
    Matrix<Base<F>>& dRowB, \
    Matrix<Base<F>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( DistSparseMatrix<F>& A, \
    DistSparseMatrix<F>& B, \
    DistMultiVec<Base<F>>& dRowA, \
    DistMultiVec<Base<F>>& dRowB, \
    DistMultiVec<Base<F>>& dCol, \
    bool progress );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
