/*
   Copyright (c) 2009-2015, Jack Poulson
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

template<typename F>
void GeomEquil
( Matrix<F>& A, Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol, bool progress )
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
    auto maxAbs = MaxAbs( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        for( Int j=0; j<n; ++j )
        {
            auto aCol = A( ALL, IR(j) );
            auto maxColAbs = VectorMaxAbs( aCol );
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
            auto maxRowAbs = VectorMaxAbs( aRow );
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

        auto newMaxAbs = MaxAbs( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "    New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each column so that its maximum entry is 1 or 0
    for( Int j=0; j<n; ++j )
    {
        auto aCol = A( ALL, IR(j) );
        auto maxColAbs = VectorMaxAbs( aCol );
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
( Matrix<F>& A, Matrix<F>& B, 
  Matrix<Base<F>>& dRowA, Matrix<Base<F>>& dRowB, 
  Matrix<Base<F>>& dCol, bool progress )
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
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        for( Int j=0; j<n; ++j )
        {
            auto aCol = A( ALL, IR(j) );
            auto bCol = B( ALL, IR(j) );
            auto maxColAbsA = VectorMaxAbs( aCol );
            auto maxColAbsB = VectorMaxAbs( bCol );
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
            auto maxRowAbs = VectorMaxAbs( aRow );
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
            auto maxRowAbs = VectorMaxAbs( bRow );
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

        auto newMaxAbsA = MaxAbs( A );
        auto newMaxAbsB = MaxAbs( B );
        const Real newMaxAbsVal = Max(newMaxAbsA.value,newMaxAbsB.value);
        const Real newMinAbsValA = MinAbsNonzero( A, newMaxAbsVal );
        const Real newMinAbsValB = MinAbsNonzero( B, newMaxAbsVal );
        const Real newMinAbsVal = Min(newMinAbsValA,newMinAbsValB);
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "    New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each column so that its maximum entry is 1 or 0
    for( Int j=0; j<n; ++j )
    {
        auto aCol = A( ALL, IR(j) );
        auto bCol = B( ALL, IR(j) );
        auto maxColAbsA = VectorMaxAbs( aCol );
        auto maxColAbsB = VectorMaxAbs( bCol );
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
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<Base<F>>& dRowPre, AbstractDistMatrix<Base<F>>& dColPre,
  bool progress )
{
    DEBUG_ONLY(CSE cse("GeomEquil"))
    typedef Base<F> Real;

    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    auto APtr    = ReadWriteProxy<F,MC,MR>(&APre,control);
    auto dRowPtr = WriteProxy<Real,MC,STAR>(&dRowPre,control); 
    auto dColPtr = WriteProxy<Real,MR,STAR>(&dColPre,control);
    auto& A = *APtr;
    auto& dRow = *dRowPtr;
    auto& dCol = *dColPtr;

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
    auto maxAbs = MaxAbs( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && A.Grid().Rank() == 0 )
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    DistMatrix<Real,MC,STAR> rowScale(A.Grid());
    DistMatrix<Real,MR,STAR> colScale(A.Grid());
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        GeometricColumnScaling( A, colScale ); 
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            if( colScale.GetLocal(jLoc,0) == Real(0) )
                colScale.SetLocal(jLoc,0,Real(1));
            dCol.SetLocal
            ( jLoc, 0, colScale.GetLocal(jLoc,0)*dCol.GetLocal(jLoc,0) );
        }
        DiagonalSolve( RIGHT, NORMAL, colScale, A );

        // Geometrically equilibrate the rows
        // ----------------------------------
        GeometricRowScaling( A, rowScale );
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
        {
            if( rowScale.GetLocal(iLoc,0) == Real(0) )
                rowScale.SetLocal(iLoc,0,Real(1));
            dRow.SetLocal
            ( iLoc, 0, rowScale.GetLocal(iLoc,0)*dRow.GetLocal(iLoc,0) );
        }
        DiagonalSolve( LEFT, NORMAL, rowScale, A );

        auto newMaxAbs = MaxAbs( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && A.Grid().Rank() == 0 )
            cout << "    New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each column so that its maximum entry is 1 or 0
    MaxEntryColumnScaling( A, colScale );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        if( colScale.GetLocal(jLoc,0) == Real(0) ) 
            colScale.SetLocal(jLoc,0,Real(1));
        dCol.SetLocal
        ( jLoc, 0, colScale.GetLocal(jLoc,0)*dCol.GetLocal(jLoc,0) );
    }
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
}

template<typename F>
void StackedGeomEquil
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<F>& BPre,
  AbstractDistMatrix<Base<F>>& dRowAPre, 
  AbstractDistMatrix<Base<F>>& dRowBPre,
  AbstractDistMatrix<Base<F>>& dColPre,
  bool progress )
{
    DEBUG_ONLY(CSE cse("GeomEquil"))
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
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    DistMatrix<Real,MC,STAR> rowScaleA(A.Grid()),
                             rowScaleB(A.Grid());
    DistMatrix<Real,MR,STAR> colScale(A.Grid());
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        StackedGeometricColumnScaling( A, B, colScale ); 
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            if( colScale.GetLocal(jLoc,0) == Real(0) )
                colScale.SetLocal(jLoc,0,Real(1));
            dCol.SetLocal
            ( jLoc, 0, colScale.GetLocal(jLoc,0)*dCol.GetLocal(jLoc,0) );
        }
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically equilibrate the rows
        // ----------------------------------
        GeometricRowScaling( A, rowScaleA );
        GeometricRowScaling( B, rowScaleB );
        for( Int iLoc=0; iLoc<mLocalA; ++iLoc )
        {
            if( rowScaleA.GetLocal(iLoc,0) == Real(0) )
                rowScaleA.SetLocal(iLoc,0,Real(1));
            dRowA.SetLocal
            ( iLoc, 0, rowScaleA.GetLocal(iLoc,0)*dRowA.GetLocal(iLoc,0) );
        }
        for( Int iLoc=0; iLoc<mLocalB; ++iLoc )
        {
            if( rowScaleB.GetLocal(iLoc,0) == Real(0) )
                rowScaleB.SetLocal(iLoc,0,Real(1));
            dRowB.SetLocal
            ( iLoc, 0, rowScaleB.GetLocal(iLoc,0)*dRowB.GetLocal(iLoc,0) );
        }
        DiagonalSolve( LEFT, NORMAL, rowScaleA, A );
        DiagonalSolve( LEFT, NORMAL, rowScaleB, B );

        auto newMaxAbsA = MaxAbs( A );
        auto newMaxAbsB = MaxAbs( B );
        const Real newMaxAbsVal = Max(newMaxAbsA.value,newMaxAbsB.value);
        const Real newMinAbsValA = MinAbsNonzero( A, newMaxAbsVal );
        const Real newMinAbsValB = MinAbsNonzero( B, newMaxAbsVal );
        const Real newMinAbsVal = Min(newMinAbsValA,newMinAbsValB);
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && A.Grid().Rank() == 0 )
            cout << "    New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each column so that its maximum entry is 1 or 0
    MaxEntryColumnScaling( A, colScale );
    MaxEntryColumnScaling( B, colScale );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        if( colScale.GetLocal(jLoc,0) == Real(0) ) 
            colScale.SetLocal(jLoc,0,Real(1));
        dCol.SetLocal
        ( jLoc, 0, colScale.GetLocal(jLoc,0)*dCol.GetLocal(jLoc,0) );
    }
    DiagonalSolve( RIGHT, NORMAL, colScale, A );
    DiagonalSolve( RIGHT, NORMAL, colScale, B );
}

template<typename F>
void GeomEquil
( SparseMatrix<F>& A, Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol,
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
    auto maxAbs = MaxAbs( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    SparseMatrix<F> ATrans;
    Transpose( A, ATrans );

    F* valBuf = A.ValueBuffer();
    F* transValBuf = ATrans.ValueBuffer();

    const Real sqrtDamp = Sqrt(damp);
    Matrix<Real> rowScale(m,1), colScale(n,1);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        for( Int j=0; j<n; ++j )
        {
            const Int offset = ATrans.EntryOffset(j);
            const Int numConnect = ATrans.NumConnections(j);

            // Compute the maximum value in this column
            Real maxColAbs = 0;
            for( Int e=offset; e<offset+numConnect; ++e )
                maxColAbs = Max(maxColAbs,Abs(ATrans.Value(e)));

            if( maxColAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this column
                Real minColAbs = maxColAbs;
                for( Int e=offset; e<offset+numConnect; ++e )
                {
                    const Real absVal = Abs(ATrans.Value(e));
                    if( absVal > 0 )
                        minColAbs = Min(minColAbs,absVal);  
                }

                const Real propScale = Sqrt(minColAbs*maxColAbs);
                const Real scale = Max(propScale,sqrtDamp*maxColAbs);
                colScale.Set(j,0,scale);
                dCol.Set( j, 0, scale*dCol.Get(j,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    transValBuf[e] /= scale;
            }
            else
                colScale.Set(j,0,Real(1));
        }
        DiagonalSolve( RIGHT, NORMAL, colScale, A );

        // Geometrically rescale the rows
        // ------------------------------ 
        for( Int i=0; i<m; ++i )
        {
            const Int offset = A.EntryOffset(i);
            const Int numConnect = A.NumConnections(i);

            // Compute the maximum value in this row
            Real maxRowAbs = 0;
            for( Int e=offset; e<offset+numConnect; ++e )
                maxRowAbs = Max(maxRowAbs,Abs(A.Value(e)));

            if( maxRowAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this row
                Real minRowAbs = maxRowAbs;
                for( Int e=offset; e<offset+numConnect; ++e )
                {
                    const Real absVal = Abs(A.Value(e));
                    if( absVal > 0 )
                        minRowAbs = Min(minRowAbs,absVal);  
                }

                const Real propScale = Sqrt(minRowAbs*maxRowAbs);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbs);
                rowScale.Set(i,0,scale);
                dRow.Set( i, 0, scale*dRow.Get(i,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    valBuf[e] /= scale;
            } 
            else
                rowScale.Set(i,0,Real(1));
        }
        DiagonalSolve( RIGHT, NORMAL, rowScale, ATrans );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbs( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "    New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each row so that its maximum entry is 1 or 0
    for( Int i=0; i<m; ++i )
    {
        const Int offset = A.EntryOffset(i);
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
( SparseMatrix<F>& A, SparseMatrix<F>& B,
  Matrix<Base<F>>& dRowA, Matrix<Base<F>>& dRowB,
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
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    SparseMatrix<F> ATrans, BTrans;
    Transpose( A, ATrans );
    Transpose( B, BTrans );

    F* valBufA = A.ValueBuffer();
    F* valBufB = B.ValueBuffer(); 
    F* transValBufA = ATrans.ValueBuffer();
    F* transValBufB = BTrans.ValueBuffer();

    const Real sqrtDamp = Sqrt(damp);
    Matrix<Real> rowScaleA(mA,1), rowScaleB(mB,1), colScale(n,1);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        for( Int j=0; j<n; ++j )
        {
            const Int offsetA = ATrans.EntryOffset(j);
            const Int offsetB = BTrans.EntryOffset(j);
            const Int numConnectA = ATrans.NumConnections(j);
            const Int numConnectB = BTrans.NumConnections(j);

            // Compute the maximum value in this column
            Real maxColAbs = 0;
            for( Int e=offsetA; e<offsetA+numConnectA; ++e )
                maxColAbs = Max(maxColAbs,Abs(ATrans.Value(e)));
            for( Int e=offsetB; e<offsetB+numConnectB; ++e )
                maxColAbs = Max(maxColAbs,Abs(BTrans.Value(e)));

            if( maxColAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this column
                Real minColAbs = maxColAbs;
                for( Int e=offsetA; e<offsetA+numConnectA; ++e )
                {
                    const Real absVal = Abs(ATrans.Value(e));
                    if( absVal > 0 )
                        minColAbs = Min(minColAbs,absVal);  
                }
                for( Int e=offsetB; e<offsetB+numConnectB; ++e )
                {
                    const Real absVal = Abs(BTrans.Value(e));
                    if( absVal > 0 )
                        minColAbs = Min(minColAbs,absVal);  
                }

                const Real propScale = Sqrt(minColAbs*maxColAbs);
                const Real scale = Max(propScale,sqrtDamp*maxColAbs);
                colScale.Set(j,0,scale);
                dCol.Set( j, 0, scale*dCol.Get(j,0) );
                for( Int e=offsetA; e<offsetA+numConnectA; ++e )
                    transValBufA[e] /= scale;
                for( Int e=offsetB; e<offsetB+numConnectB; ++e )
                    transValBufB[e] /= scale;
            }
            else
                colScale.Set(j,0,Real(1));
        }
        DiagonalSolve( RIGHT, NORMAL, colScale, A );
        DiagonalSolve( RIGHT, NORMAL, colScale, B );

        // Geometrically rescale the rows
        // ------------------------------ 
        for( Int i=0; i<mA; ++i )
        {
            const Int offset = A.EntryOffset(i);
            const Int numConnect = A.NumConnections(i);

            // Compute the maximum value in this row
            Real maxRowAbs = 0;
            for( Int e=offset; e<offset+numConnect; ++e )
                maxRowAbs = Max(maxRowAbs,Abs(A.Value(e)));

            if( maxRowAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this row
                Real minRowAbs = maxRowAbs;
                for( Int e=offset; e<offset+numConnect; ++e )
                {
                    const Real absVal = Abs(A.Value(e));
                    if( absVal > 0 )
                        minRowAbs = Min(minRowAbs,absVal);  
                }

                const Real propScale = Sqrt(minRowAbs*maxRowAbs);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbs);
                rowScaleA.Set(i,0,scale);
                dRowA.Set( i, 0, scale*dRowA.Get(i,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    valBufA[e] /= scale;
            } 
            else
                rowScaleA.Set(i,0,Real(1));
        }
        for( Int i=0; i<mB; ++i )
        {
            const Int offset = B.EntryOffset(i);
            const Int numConnect = B.NumConnections(i);

            // Compute the maximum value in this row
            Real maxRowAbs = 0;
            for( Int e=offset; e<offset+numConnect; ++e )
                maxRowAbs = Max(maxRowAbs,Abs(B.Value(e)));

            if( maxRowAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this row
                Real minRowAbs = maxRowAbs;
                for( Int e=offset; e<offset+numConnect; ++e )
                {
                    const Real absVal = Abs(B.Value(e));
                    if( absVal > 0 )
                        minRowAbs = Min(minRowAbs,absVal);  
                }

                const Real propScale = Sqrt(minRowAbs*maxRowAbs);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbs);
                rowScaleB.Set(i,0,scale);
                dRowB.Set( i, 0, scale*dRowB.Get(i,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    valBufB[e] /= scale;
            } 
            else
                rowScaleB.Set(i,0,Real(1));
        }
        DiagonalSolve( RIGHT, NORMAL, rowScaleA, ATrans );
        DiagonalSolve( RIGHT, NORMAL, rowScaleB, BTrans );

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
            cout << "    New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each row so that its maximum entry is 1 or 0
    for( Int i=0; i<mA; ++i )
    {
        const Int offset = A.EntryOffset(i);
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
    for( Int i=0; i<mB; ++i )
    {
        const Int offset = B.EntryOffset(i);
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
  DistMultiVec<Base<F>>& dRow, DistMultiVec<Base<F>>& dCol, 
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
    auto maxAbs = MaxAbs( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && commRank == 0 )
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    DistSparseMatrix<F> ATrans(comm);

    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        Transpose( A, ATrans );
        F* transValBuf = ATrans.ValueBuffer();
        const Int localWidth = ATrans.LocalHeight();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int offset = ATrans.EntryOffset(jLoc);
            const Int numConnect = ATrans.NumConnections(jLoc);

            // Compute the maximum value in this column
            Real maxColAbs = 0;
            for( Int e=offset; e<offset+numConnect; ++e )
                maxColAbs = Max(maxColAbs,Abs(ATrans.Value(e)));

            if( maxColAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this column
                Real minColAbs = maxColAbs;
                for( Int e=offset; e<offset+numConnect; ++e )
                {
                    const Real absVal = Abs(ATrans.Value(e));
                    if( absVal > 0 )
                        minColAbs = Min(minColAbs,absVal);  
                }

                const Real propScale = Sqrt(minColAbs*maxColAbs);
                const Real scale = Max(propScale,sqrtDamp*maxColAbs);
                dCol.SetLocal( jLoc, 0, scale*dCol.GetLocal(jLoc,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    transValBuf[e] /= scale;
            }
        }

        // Geometrically rescale the rows
        // ------------------------------ 
        Transpose( ATrans, A );
        F* valBuf = A.ValueBuffer();
        const Int localHeight = A.LocalHeight();
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int offset = A.EntryOffset(iLoc);
            const Int numConnect = A.NumConnections(iLoc);

            // Compute the maximum value in this row
            Real maxRowAbs = 0;
            for( Int e=offset; e<offset+numConnect; ++e )
                maxRowAbs = Max(maxRowAbs,Abs(A.Value(e)));

            if( maxRowAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this row
                Real minRowAbs = maxRowAbs;
                for( Int e=offset; e<offset+numConnect; ++e )
                {
                    const Real absVal = Abs(A.Value(e));
                    if( absVal > 0 )
                        minRowAbs = Min(minRowAbs,absVal);  
                }

                const Real propScale = Sqrt(minRowAbs*maxRowAbs);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbs);
                dRow.SetLocal( iLoc, 0, scale*dRow.GetLocal(iLoc,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    valBuf[e] /= scale;
            }
        }

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbs( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && commRank == 0 )
            cout << "    New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each row so that its maximum entry is 1 or 0
    F* valBuf = A.ValueBuffer();
    const Int localHeight = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int offset = A.EntryOffset(iLoc);
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
( DistSparseMatrix<F>& A, DistSparseMatrix<F>& B,
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
        cout << "    Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    DistSparseMatrix<F> ATrans(comm), BTrans(comm);

    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        Transpose( A, ATrans );
        Transpose( B, BTrans );
        F* transValBufA = ATrans.ValueBuffer();
        F* transValBufB = BTrans.ValueBuffer();
        const Int localWidth = ATrans.LocalHeight();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
        {
            const Int offsetA = ATrans.EntryOffset(jLoc);
            const Int offsetB = BTrans.EntryOffset(jLoc);
            const Int numConnectA = ATrans.NumConnections(jLoc);
            const Int numConnectB = BTrans.NumConnections(jLoc);

            // Compute the maximum value in this column
            Real maxColAbs = 0;
            for( Int e=offsetA; e<offsetA+numConnectA; ++e )
                maxColAbs = Max(maxColAbs,Abs(ATrans.Value(e)));
            for( Int e=offsetB; e<offsetB+numConnectB; ++e )
                maxColAbs = Max(maxColAbs,Abs(BTrans.Value(e)));

            if( maxColAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this column
                Real minColAbs = maxColAbs;
                for( Int e=offsetA; e<offsetA+numConnectA; ++e )
                {
                    const Real absVal = Abs(ATrans.Value(e));
                    if( absVal > 0 )
                        minColAbs = Min(minColAbs,absVal);  
                }
                for( Int e=offsetB; e<offsetB+numConnectB; ++e )
                {
                    const Real absVal = Abs(BTrans.Value(e));
                    if( absVal > 0 )
                        minColAbs = Min(minColAbs,absVal);  
                }

                const Real propScale = Sqrt(minColAbs*maxColAbs);
                const Real scale = Max(propScale,sqrtDamp*maxColAbs);
                dCol.SetLocal( jLoc, 0, scale*dCol.GetLocal(jLoc,0) );
                for( Int e=offsetA; e<offsetA+numConnectA; ++e )
                    transValBufA[e] /= scale;
                for( Int e=offsetB; e<offsetB+numConnectB; ++e )
                    transValBufB[e] /= scale;
            }
        }

        // Geometrically rescale the rows
        // ------------------------------ 
        Transpose( ATrans, A );
        Transpose( BTrans, B );
        F* valBufA = A.ValueBuffer();
        F* valBufB = B.ValueBuffer();
        const Int localHeightA = A.LocalHeight();
        for( Int iLoc=0; iLoc<localHeightA; ++iLoc )
        {
            const Int offset = A.EntryOffset(iLoc);
            const Int numConnect = A.NumConnections(iLoc);

            // Compute the maximum value in this row
            Real maxRowAbs = 0;
            for( Int e=offset; e<offset+numConnect; ++e )
                maxRowAbs = Max(maxRowAbs,Abs(A.Value(e)));

            if( maxRowAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this row
                Real minRowAbs = maxRowAbs;
                for( Int e=offset; e<offset+numConnect; ++e )
                {
                    const Real absVal = Abs(A.Value(e));
                    if( absVal > 0 )
                        minRowAbs = Min(minRowAbs,absVal);  
                }

                const Real propScale = Sqrt(minRowAbs*maxRowAbs);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbs);
                dRowA.SetLocal( iLoc, 0, scale*dRowA.GetLocal(iLoc,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    valBufA[e] /= scale;
            }
        }
        const Int localHeightB = B.LocalHeight();
        for( Int iLoc=0; iLoc<localHeightB; ++iLoc )
        {
            const Int offset = B.EntryOffset(iLoc);
            const Int numConnect = B.NumConnections(iLoc);

            // Compute the maximum value in this row
            Real maxRowAbs = 0;
            for( Int e=offset; e<offset+numConnect; ++e )
                maxRowAbs = Max(maxRowAbs,Abs(B.Value(e)));

            if( maxRowAbs > Real(0) )
            {
                // Compute the minimum nonzero value in this row
                Real minRowAbs = maxRowAbs;
                for( Int e=offset; e<offset+numConnect; ++e )
                {
                    const Real absVal = Abs(B.Value(e));
                    if( absVal > 0 )
                        minRowAbs = Min(minRowAbs,absVal);  
                }

                const Real propScale = Sqrt(minRowAbs*maxRowAbs);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbs);
                dRowB.SetLocal( iLoc, 0, scale*dRowB.GetLocal(iLoc,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    valBufB[e] /= scale;
            }
        }

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
            cout << "    New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each row so that its maximum entry is 1 or 0
    F* valBufA = A.ValueBuffer();
    const Int localHeightA = A.LocalHeight();
    for( Int iLoc=0; iLoc<localHeightA; ++iLoc )
    {
        const Int offset = A.EntryOffset(iLoc);
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
        const Int offset = B.EntryOffset(iLoc);
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
  ( Matrix<F>& A, Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Base<F>>& dRow, AbstractDistMatrix<Base<F>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( SparseMatrix<F>& A, Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol, \
    bool progress ); \
  template void GeomEquil \
  ( DistSparseMatrix<F>& A, \
    DistMultiVec<Base<F>>& dRow, DistMultiVec<Base<F>>& dCol, \
    bool progress ); \
  template void StackedGeomEquil \
  ( Matrix<F>& A, Matrix<F>& B, \
    Matrix<Base<F>>& dRowA, Matrix<Base<F>>& dRowB, \
    Matrix<Base<F>>& dCol, bool progress ); \
  template void StackedGeomEquil \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, \
    AbstractDistMatrix<Base<F>>& dRowA, \
    AbstractDistMatrix<Base<F>>& dRowB, \
    AbstractDistMatrix<Base<F>>& dCol, bool progress ); \
  template void StackedGeomEquil \
  ( SparseMatrix<F>& A, SparseMatrix<F>& B, \
    Matrix<Base<F>>& dRowA, Matrix<Base<F>>& dRowB, \
    Matrix<Base<F>>& dCol, bool progress ); \
  template void StackedGeomEquil \
  ( DistSparseMatrix<F>& A, DistSparseMatrix<F>& B, \
    DistMultiVec<Base<F>>& dRowA, DistMultiVec<Base<F>>& dRowB, \
    DistMultiVec<Base<F>>& dCol, bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
