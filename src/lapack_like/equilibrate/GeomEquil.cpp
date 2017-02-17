/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
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

// TODO(poulson): Make this consistent with ConeGeomEquil

template<typename Field>
void GeomEquil
( Matrix<Field>& A,
  Matrix<Base<Field>>& dRow,
  Matrix<Base<Field>>& dCol,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
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
                dCol(j) *= scale;
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
                dRow(i) *= scale;
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
            dCol(j) *= maxColAbsVal;
        }
    }
}

template<typename Field>
void StackedGeomEquil
( Matrix<Field>& A,
  Matrix<Field>& B,
  Matrix<Base<Field>>& dRowA,
  Matrix<Base<Field>>& dRowB,
  Matrix<Base<Field>>& dCol,
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
                dCol(j) *= scale;
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
                dRowA(i) *= scale;
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
                dRowB(i) *= scale;
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
            dCol(j) *= maxColAbsVal;
        }
    }
}

template<typename Field>
void GeomEquil
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Base<Field>>& dRowPre,
  AbstractDistMatrix<Base<Field>>& dColPre,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre, control );
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

    // TODO(poulson): Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 6;
    const Real relTol = Real(9)/Real(10);

    // TODO(poulson): Incorporate damping
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
    auto& colScaleLoc = colScale.Matrix();
    auto& rowScaleLoc = rowScale.Matrix();
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        // TODO(poulson): Remove GeometricColumnScaling
        GeometricColumnScaling( A, colScale );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            if( colScaleLoc(jLoc) == Real(0) )
                colScaleLoc(jLoc) = Real(1);
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );

        // Geometrically equilibrate the rows
        // ----------------------------------
        // TODO(poulson): Remove GeometricRowScaling
        GeometricRowScaling( A, rowScale );
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            if( rowScaleLoc(iLoc) == Real(0) )
                rowScaleLoc(iLoc) = Real(1);
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
        if( colScaleLoc(jLoc) == Real(0) )
            colScaleLoc(jLoc) = Real(1);
    DiagonalScale( LEFT, NORMAL, colScale, dCol );
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
}

template<typename Field>
void StackedGeomEquil
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& BPre,
  AbstractDistMatrix<Base<Field>>& dRowAPre,
  AbstractDistMatrix<Base<Field>>& dRowBPre,
  AbstractDistMatrix<Base<Field>>& dColPre,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre, control );
    DistMatrixReadWriteProxy<Field,Field,MC,MR> BProx( BPre, control );
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
    auto& rowScaleALoc = rowScaleA.Matrix();
    auto& rowScaleBLoc = rowScaleB.Matrix();
    auto& colScaleLoc = colScale.Matrix();
    auto& colScaleBLoc = colScaleB.Matrix();
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        // TODO(poulson): Remove StackedGeometricColumnScaling
        StackedGeometricColumnScaling( A, B, colScale );
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            if( colScaleLoc(jLoc) == Real(0) )
                colScaleLoc(jLoc) = Real(1);
        DiagonalScale( LEFT, NORMAL, colScale, dCol );
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically equilibrate the rows
        // ----------------------------------
        // TODO(poulson): Remove GeometricRowScaling
        GeometricRowScaling( A, rowScaleA );
        for( Int iLoc=0; iLoc<mLocalA; ++iLoc )
            if( rowScaleALoc(iLoc) == Real(0) )
                rowScaleALoc(iLoc) = Real(1);
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        // TODO(poulson): Remove GeometricRowScaling
        GeometricRowScaling( B, rowScaleB );
        for( Int iLoc=0; iLoc<mLocalB; ++iLoc )
            if( rowScaleBLoc(iLoc) == Real(0) )
                rowScaleBLoc(iLoc) = Real(1);
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
        Real maxScale = Max(colScaleLoc(jLoc),colScaleBLoc(jLoc));
        if( maxScale == Real(0) )
            maxScale = 1;
        colScaleLoc(jLoc) = maxScale;
    }
    DiagonalScale( LEFT, NORMAL, colScale, dCol );
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
    DiagonalSolve( RIGHT, NORMAL, colScale, B );
}

template<typename Field>
void GeomEquil
( SparseMatrix<Field>& A,
  Matrix<Base<Field>>& dRow,
  Matrix<Base<Field>>& dCol,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
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
            const Real maxAbs = maxAbsVals(j);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsVals(j);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                colScale(j) = scale;
            }
            else
            {
                colScale(j) = 1;
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
            const Real maxAbs = maxAbsVals(i);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsVals(i);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                rowScale(i) = scale;
            }
            else
                rowScale(i) = Real(1);
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
    Field* valBuf = A.ValueBuffer();
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
            dRow(i) *= maxRowAbs;
            for( Int e=offset; e<offset+numConnect; ++e )
                valBuf[e] /= maxRowAbs;
        }
    }
}

template<typename Field>
void StackedGeomEquil
( SparseMatrix<Field>& A,
  SparseMatrix<Field>& B,
  Matrix<Base<Field>>& dRowA,
  Matrix<Base<Field>>& dRowB,
  Matrix<Base<Field>>& dCol,
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
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        ColumnMaxNorms( A, maxAbsValsA );
        ColumnMaxNorms( B, maxAbsValsB );
        for( Int j=0; j<n; ++j )
            maxAbsValsA(j) = Max(maxAbsValsA(j),maxAbsValsB(j));
        ColumnMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        ColumnMinAbsNonzero( B, maxAbsValsA, minAbsValsB );
        for( Int j=0; j<n; ++j )
            minAbsValsA(j) = Min(minAbsValsA(j),minAbsValsB(j));
        for( Int j=0; j<n; ++j )
        {
            const Real maxAbs = maxAbsValsA(j);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsA(j);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                colScale(j) = scale;
            }
            else
                colScale(j) = Real(1);
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
            const Real maxAbs = maxAbsValsA(i);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsA(i);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                rowScaleA(i) = scale;
            }
            else
                rowScaleA(i) = Real(1);
        }
        DiagonalScale( LEFT, NORMAL, rowScaleA, dRowA );
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );

        RowMinAbsNonzero( B, maxAbsValsB, minAbsValsB );
        RowMaxNorms( B, maxAbsValsB );
        for( Int i=0; i<mB; ++i )
        {
            const Real maxAbs = maxAbsValsB(i);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsB(i);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                rowScaleB(i) = scale;
            }
            else
                rowScaleB(i) = Real(1);
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
    Field* valBufA = A.ValueBuffer();
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
            dRowA(i) *= maxRowAbs;
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufA[e] /= maxRowAbs;
        }
    }
    Field* valBufB = B.ValueBuffer();
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
            dRowB(i) *= maxRowAbs;
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufB[e] /= maxRowAbs;
        }
    }
}

template<typename Field>
void GeomEquil
( DistSparseMatrix<Field>& A,
  DistMultiVec<Base<Field>>& dRow,
  DistMultiVec<Base<Field>>& dCol,
  bool progress )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const int commRank = grid.Rank();
    dRow.SetGrid(grid);
    dCol.SetGrid(grid);
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO(poulson): Expose these as control parameters
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
    DistMultiVec<Real> maxAbsVals(grid), minAbsVals(grid), scales(grid);
    auto& scalesLoc = scales.Matrix();
    const Int indent = PushIndent();
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        ColumnMaxNorms( A, maxAbsVals );
        ColumnMinAbsNonzero( A, maxAbsVals, minAbsVals );
        auto& maxAbsValsLoc = maxAbsVals.Matrix();
        auto& minAbsValsLoc = minAbsVals.Matrix();
        scales.Resize( n, 1 );
        const Int localWidth = maxAbsVals.LocalHeight();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Real maxAbs = maxAbsValsLoc(jLoc);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsLoc(jLoc);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scalesLoc(jLoc) = scale;
            }
            else
                scalesLoc(jLoc) = 1;
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
            const Real maxAbs = maxAbsValsLoc(iLoc);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsLoc(iLoc);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scalesLoc(iLoc) = scale;
            }
            else
                scalesLoc(iLoc) = 1;
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
    Field* valBuf = A.ValueBuffer();
    const Int localHeight = A.LocalHeight();
    auto& dRowLoc = dRow.Matrix();
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
            dRowLoc(iLoc) *= maxRowAbs;
            for( Int e=offset; e<offset+numConnect; ++e )
                valBuf[e] /= maxRowAbs;
        }
    }
}

template<typename Field>
void StackedGeomEquil
( DistSparseMatrix<Field>& A,
  DistSparseMatrix<Field>& B,
  DistMultiVec<Base<Field>>& dRowA,
  DistMultiVec<Base<Field>>& dRowB,
  DistMultiVec<Base<Field>>& dCol,
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
    auto& maxAbsValsALoc = maxAbsValsA.Matrix();
    auto& maxAbsValsBLoc = maxAbsValsB.Matrix();
    auto& minAbsValsALoc = minAbsValsA.Matrix();
    auto& minAbsValsBLoc = minAbsValsB.Matrix();
    auto& scalesLoc = scales.Matrix();
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
            maxAbsValsALoc(jLoc) =
              Max(maxAbsValsALoc(jLoc),maxAbsValsBLoc(jLoc));
        ColumnMinAbsNonzero( A, maxAbsValsA, minAbsValsA );
        ColumnMinAbsNonzero( B, maxAbsValsA, minAbsValsB );
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            minAbsValsALoc(jLoc) =
              Min(minAbsValsALoc(jLoc),minAbsValsBLoc(jLoc));
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Real maxAbs = maxAbsValsALoc(jLoc);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsALoc(jLoc);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scalesLoc(jLoc) = scale;
            }
            else
                scalesLoc(jLoc) = 1;
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
            const Real maxAbs = maxAbsValsALoc(iLoc);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsALoc(iLoc);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scalesLoc(iLoc) = scale;
            }
            else
                scalesLoc(iLoc) = 1;
        }
        DiagonalScale( LEFT, NORMAL, scales, dRowA );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        scales.Resize( mB, 1 );
        RowMaxNorms( B, maxAbsValsB );
        ColumnMinAbsNonzero( B, maxAbsValsB, minAbsValsB );
        const Int localHeightB = maxAbsValsB.LocalHeight();
        for( Int iLoc=0; iLoc<localHeightB; ++iLoc )
        {
            const Real maxAbs = maxAbsValsBLoc(iLoc);
            if( maxAbs > Real(0) )
            {
                const Real minAbs = minAbsValsBLoc(iLoc);
                const Real propScale = Sqrt(minAbs*maxAbs);
                const Real scale = Max(propScale,sqrtDamp*maxAbs);
                scalesLoc(iLoc) = scale;
            }
            else
                scalesLoc(iLoc) = 1;
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
    Field* valBufA = A.ValueBuffer();
    const Int localHeightA = A.LocalHeight();
    auto& dRowALoc = dRowA.Matrix();
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
            dRowALoc(iLoc) *= maxRowAbs;
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufA[e] /= maxRowAbs;
        }
    }
    Field* valBufB = B.ValueBuffer();
    const Int localHeightB = B.LocalHeight();
    auto& dRowBLoc = dRowB.Matrix();
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
            dRowBLoc(iLoc) *= maxRowAbs;
            for( Int e=offset; e<offset+numConnect; ++e )
                valBufB[e] /= maxRowAbs;
        }
    }
}

#define PROTO(Field) \
  template void GeomEquil \
  ( Matrix<Field>& A, \
    Matrix<Base<Field>>& dRow, \
    Matrix<Base<Field>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Base<Field>>& dRow, \
    AbstractDistMatrix<Base<Field>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( SparseMatrix<Field>& A, \
    Matrix<Base<Field>>& dRow, \
    Matrix<Base<Field>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( DistSparseMatrix<Field>& A, \
    DistMultiVec<Base<Field>>& dRow, \
    DistMultiVec<Base<Field>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( Matrix<Field>& A, \
    Matrix<Field>& B, \
    Matrix<Base<Field>>& dRowA, \
    Matrix<Base<Field>>& dRowB, \
    Matrix<Base<Field>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Field>& B, \
    AbstractDistMatrix<Base<Field>>& dRowA, \
    AbstractDistMatrix<Base<Field>>& dRowB, \
    AbstractDistMatrix<Base<Field>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( SparseMatrix<Field>& A, \
    SparseMatrix<Field>& B, \
    Matrix<Base<Field>>& dRowA, \
    Matrix<Base<Field>>& dRowB, \
    Matrix<Base<Field>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( DistSparseMatrix<Field>& A, \
    DistSparseMatrix<Field>& B, \
    DistMultiVec<Base<Field>>& dRowA, \
    DistMultiVec<Base<Field>>& dRowB, \
    DistMultiVec<Base<Field>>& dCol, \
    bool progress );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
