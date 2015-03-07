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
void SymmetricGeomEquil( Matrix<F>& A, Matrix<Base<F>>& d, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricGeomEquil"))
    // TODO: Ensure A is symmetric
    typedef Base<F> Real;
    const Int n = A.Height();
    Ones( d, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10; 
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
        cout << "Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    Matrix<Real> scales;
    Zeros( scales, n, 1 );
    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        for( Int j=0; j<n; ++j )
        {
            auto aCol = A( IR(0,n), IR(j,j+1) );
            auto maxColAbs = VectorMaxAbs( aCol );
            const Real maxColAbsVal = maxColAbs.value;
            if( maxColAbsVal > Real(0) )
            {
                const Real minColAbsVal = MinAbsNonzero( aCol, maxColAbsVal );
                const Real propScale = Sqrt(minColAbsVal*maxColAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxColAbsVal);
                scales.Set( j, 0, scale );
                d.Set( j, 0, scale*d.Get(j,0) );
            }
            else
                scales.Set( j, 0, Real(1) );
        }
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        auto newMaxAbs = MaxAbs( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Normalize the maximum absolute values towards one
    for( Int normIter=0; normIter<3; ++normIter )
    {
        for( Int j=0; j<n; ++j )
        {
            Real maxAbs = 1;
            for( Int i=0; i<n; ++i )
                maxAbs = Max( Abs(A.Get(i,j)), maxAbs );
            const Real scale = Sqrt(maxAbs);
            scales.Set( j, 0, scale );
            d.Set( j, 0, scale*d.Get(j,0) );
        }
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );
    }
    auto newMaxAbs = MaxAbs( A );
    const Real newMaxAbsVal = newMaxAbs.value;
    const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
    const Real newRatio = newMaxAbsVal / newMinAbsVal;
    if( progress )
        cout << "Final ratio is " << newMaxAbsVal << "/" 
             << newMinAbsVal << "=" << newRatio << endl;
}

template<typename F>
void SymmetricGeomEquil
( AbstractDistMatrix<F>& APre, 
  AbstractDistMatrix<Base<F>>& dPre, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricGeomEquil"))
    typedef Base<F> Real;

    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    auto APtr    = ReadWriteProxy<F,MC,MR>(&APre,control);
    auto dPtr = WriteProxy<Real,MC,STAR>(&dPre,control); 
    auto& A = *APtr;
    auto& d = *dPtr;

    const Int n = A.Height();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    Ones( d, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbs( A );
    const Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    const Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && A.Grid().Rank() == 0 )
        cout << "Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    const Real sqrtDamp = Sqrt(damp);
    DistMatrix<Real,MR,STAR> scales(A.Grid());
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        // -------------------------------------
        GeometricColumnScaling( A, scales ); 
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            if( scales.GetLocal(jLoc,0) == Real(0) )
                scales.SetLocal(jLoc,0,Real(1));
            d.SetLocal
            ( jLoc, 0, scales.GetLocal(jLoc,0)*d.GetLocal(jLoc,0) );
        }
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        auto newMaxAbs = MaxAbs( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && A.Grid().Rank() == 0 )
            cout << "New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Normalize the maximum absolute values towards one
    scales.AlignWith( A );
    Zeros( scales, n, 1 );
    for( Int normIter=0; normIter<3; ++normIter )
    {
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            Real maxAbs = 1;
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
                maxAbs = Max( Abs(A.GetLocal(iLoc,jLoc)), maxAbs );
            const Real scale = Sqrt(maxAbs);
            scales.SetLocal( jLoc, 0, scale );
        }
        mpi::AllReduce( scales.Buffer(), nLocal, mpi::MAX, A.ColComm() );
        DiagonalScale( LEFT, NORMAL, scales, d );
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );
    }
    auto newMaxAbs = MaxAbs( A );
    const Real newMaxAbsVal = newMaxAbs.value;
    const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
    const Real newRatio = newMaxAbsVal / newMinAbsVal;
    if( progress && A.Grid().Rank() == 0 ) 
        cout << "Final ratio is " << newMaxAbsVal << "/"
             << newMinAbsVal << "=" << newRatio << endl;
}

template<typename F>
void SymmetricGeomEquil
( SparseMatrix<F>& A, Matrix<Base<F>>& d, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricGeomEquil"))
    typedef Base<F> Real;
    const Int n = A.Height();
    const Int numEntries = A.NumEntries();
    Ones( d, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10; 
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
        cout << "Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    SparseMatrix<F> ATrans;
    const Real sqrtDamp = Sqrt(damp);
    Matrix<Real> scales(n,1);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        ATrans = A; // NOTE: Same as Transpose( A, ATrans ) by assumption
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
                scales.Set(j,0,scale);
                d.Set( j, 0, scale*d.Get(j,0) );
            }
            else
                scales.Set(j,0,Real(1));
        }
        DiagonalSolve( LEFT, NORMAL, scales, ATrans );
        Transpose( ATrans, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbs( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Normalize the maximum absolute values towards one
    for( Int normIter=0; normIter<3; ++normIter )
    {
        for( Int i=0; i<n; ++i )
        {
            Real maxAbs = 1;
            const Int numConnections = A.NumConnections(i);
            const Int entryOff = A.EntryOffset(i);
            for( Int k=0; k<numConnections; ++k )
                maxAbs = Max( Abs(A.Value(entryOff+k)), maxAbs );
            const Real scale = Sqrt(maxAbs);
            scales.Set( i, 0, scale );
            d.Set( i, 0, scale*d.Get(i,0) );
        }
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );
    }
    auto newMaxAbs = MaxAbs( A );
    const Real newMaxAbsVal = newMaxAbs.value;
    const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
    const Real newRatio = newMaxAbsVal / newMinAbsVal;
    if( progress )
        cout << "Final ratio is " << newMaxAbsVal << "/" 
             << newMinAbsVal << "=" << newRatio << endl;
}

template<typename F>
void SymmetricGeomEquil
( DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& d, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("SymmetricGeomEquil"))
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank(comm);
    d.SetComm( comm );
    Ones( d, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10; 
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
        cout << "Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    DistSparseMatrix<F> ATrans(comm);

    DistMultiVec<Real> scales(comm);
    Zeros( scales, n, 1 );

    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically rescale the columns
        // ---------------------------------
        ATrans = A; // NOTE: Same as Transpose( A, ATrans ) by assumption
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
                scales.SetLocal( jLoc, 0, scale );
                d.SetLocal( jLoc, 0, scale*d.GetLocal(jLoc,0) );
            }
            else
                scales.SetLocal( jLoc, 0, Real(1) );
        }
        DiagonalSolve( LEFT, NORMAL, scales, ATrans );
        Transpose( ATrans, A ); 
        DiagonalSolve( LEFT, NORMAL, scales, A );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbs( A );
        const Real newMaxAbsVal = newMaxAbs.value;
        const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        const Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && commRank == 0 )
            cout << "New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Normalize the maximum absolute values towards one
    const Int localHeight = A.LocalHeight();
    for( Int normIter=0; normIter<3; ++normIter )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            Real maxAbs = 1;
            const Int numConnections = A.NumConnections(iLoc);
            const Int entryOff = A.EntryOffset(iLoc);
            for( Int k=0; k<numConnections; ++k )
                maxAbs = Max( Abs(A.Value(entryOff+k)), maxAbs );
            const Real scale = Sqrt(maxAbs);
            scales.SetLocal( iLoc, 0, scale );
            d.SetLocal( iLoc, 0, scale*d.GetLocal(iLoc,0) );
        }
        DiagonalSolve( RIGHT, NORMAL, scales, A );
        DiagonalSolve( LEFT, NORMAL, scales, A );
    }
    auto newMaxAbs = MaxAbs( A );
    const Real newMaxAbsVal = newMaxAbs.value;
    const Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
    const Real newRatio = newMaxAbsVal / newMinAbsVal;
    if( progress && commRank == 0 )
        cout << "Final ratio is " << newMaxAbsVal << "/" 
             << newMinAbsVal << "=" << newRatio << endl;
}

#define PROTO(F) \
  template void SymmetricGeomEquil \
  ( Matrix<F>& A, Matrix<Base<F>>& d, bool progress ); \
  template void SymmetricGeomEquil \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Base<F>>& d, bool progress ); \
  template void SymmetricGeomEquil \
  ( SparseMatrix<F>& A, Matrix<Base<F>>& d, bool progress ); \
  template void SymmetricGeomEquil \
  ( DistSparseMatrix<F>& A, DistMultiVec<Base<F>>& d, bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
