/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

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
inline Base<F> MinAbsNonzero( const Matrix<F>& A, Base<F> upperBound )
{
    DEBUG_ONLY(CallStackEntry cse("MinAbsNonzero"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Real minAbs = upperBound;
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Real alphaAbs = Abs(A.Get(i,j));
            if( alphaAbs > Real(0) && alphaAbs < minAbs )
                minAbs = alphaAbs;
        }
    }
    return minAbs;
}

template<typename F>
inline Base<F> MinAbsNonzero( const SparseMatrix<F>& A, Base<F> upperBound )
{
    DEBUG_ONLY(CallStackEntry cse("MinAbsNonzero"))
    typedef Base<F> Real;
    const Int numEntries = A.NumEntries();
    Real minAbs = upperBound;
    for( Int e=0; e<numEntries; ++e )
    {
        const Real absVal = Abs(A.Value(e));
        if( absVal > Real(0) && absVal < minAbs )
            minAbs = absVal;
    }
    return minAbs;
}

template<typename F>
inline Base<F> MinAbsNonzero
( const AbstractDistMatrix<F>& A, Base<F> upperBound )
{
    DEBUG_ONLY(CallStackEntry cse("MinAbsNonzero"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Real minAbs;
    if( A.Participating() )
    {
        const Real minLocAbs = MinAbsNonzero( A.LockedMatrix(), upperBound );
        minAbs = mpi::AllReduce( minLocAbs, mpi::MAX, A.DistComm() ); 
    }
    mpi::Broadcast( minAbs, A.Root(), A.CrossComm() );
    return minAbs;
}

template<typename F>
inline Base<F> MinAbsNonzero( const DistSparseMatrix<F>& A, Base<F> upperBound )
{
    DEBUG_ONLY(CallStackEntry cse("MinAbsNonzero"))
    typedef Base<F> Real;
    const Int numEntries = A.NumLocalEntries();
    Real minLocAbs = upperBound;
    for( Int e=0; e<numEntries; ++e )
    {
        const Real absVal = Abs(A.Value(e));
        if( absVal > Real(0) && absVal < minLocAbs )
            minLocAbs = absVal;
    }
    return mpi::AllReduce( minLocAbs, mpi::MAX, A.Comm() );
}

template<typename F,Dist U,Dist V>
inline void MaxEntryColumnScaling
( const DistMatrix<F,      U,V   >& A, 
        DistMatrix<Base<F>,V,STAR>& scaling )
{
    DEBUG_ONLY(CallStackEntry cse("MaxEntryColumnScaling"))
    typedef Base<F> Real;
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    scaling.AlignWith( A );
    scaling.Resize( A.Width(), 1 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real maxAbs = 0;
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            maxAbs = Max(maxAbs,Abs(A.GetLocal(iLoc,jLoc)));
        scaling.SetLocal( jLoc, 0, maxAbs );
    }
    mpi::AllReduce( scaling.Buffer(), nLocal, mpi::MAX, A.ColComm() );
}

template<typename F,Dist U,Dist V>
inline void GeometricColumnScaling
( const DistMatrix<F,      U,V   >& A, 
        DistMatrix<Base<F>,V,STAR>& geomScaling )
{
    DEBUG_ONLY(CallStackEntry cse("GeometricColumnScaling"))
    typedef Base<F> Real;

    DistMatrix<Real,V,STAR> maxScaling(A.Grid());
    MaxEntryColumnScaling( A, maxScaling );

    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    geomScaling.AlignWith( maxScaling );
    geomScaling.Resize( A.Width(), 1 );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        Real minAbs = maxScaling.GetLocal(jLoc,0);
        for( Int iLoc=0; iLoc<mLocal; ++iLoc )
        {
            const Real absVal = Abs(A.GetLocal(iLoc,jLoc));  
            if( absVal > 0 && absVal < minAbs )
                minAbs = Min(minAbs,absVal);
        }
        geomScaling.SetLocal( jLoc, 0, minAbs );
    }
    mpi::AllReduce( geomScaling.Buffer(), nLocal, mpi::MIN, A.ColComm() );

    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        const Real maxAbs = maxScaling.GetLocal(iLoc,0);
        const Real minAbs = geomScaling.GetLocal(iLoc,0);
        geomScaling.SetLocal(iLoc,0,Sqrt(minAbs*maxAbs));
    }
}

template<typename F,Dist U,Dist V>
inline void MaxEntryRowScaling
( const DistMatrix<F,      U,V   >& A, 
        DistMatrix<Base<F>,U,STAR>& scaling )
{
    DEBUG_ONLY(CallStackEntry cse("MaxEntryRowScaling"))
    typedef Base<F> Real;
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    scaling.AlignWith( A );

    scaling.Resize( A.Height(), 1 );
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        Real maxAbs = 0;
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
            maxAbs = Max(maxAbs,Abs(A.GetLocal(iLoc,jLoc)));
        scaling.SetLocal( iLoc, 0, maxAbs );
    }
    mpi::AllReduce( scaling.Buffer(), mLocal, mpi::MAX, A.RowComm() );
}

template<typename F,Dist U,Dist V>
inline void GeometricRowScaling
( const DistMatrix<F,      U,V   >& A, 
        DistMatrix<Base<F>,U,STAR>& geomScaling )
{
    DEBUG_ONLY(CallStackEntry cse("GeometricRowScaling"))
    typedef Base<F> Real;

    DistMatrix<Real,U,STAR> maxScaling(A.Grid());
    MaxEntryRowScaling( A, maxScaling );

    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    geomScaling.AlignWith( maxScaling );
    geomScaling.Resize( A.Height(), 1 );
    for( Int iLoc=0; iLoc<mLocal; ++iLoc )
    {
        Real minAbs = maxScaling.GetLocal(iLoc,0);
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Real absVal = Abs(A.GetLocal(iLoc,jLoc));  
            if( absVal > 0 && absVal < minAbs )
                minAbs = Min(minAbs,absVal);
        }
        geomScaling.SetLocal( iLoc, 0, minAbs );
    }
    mpi::AllReduce( geomScaling.Buffer(), mLocal, mpi::MIN, A.RowComm() );

    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Real maxAbs = maxScaling.GetLocal(jLoc,0);
        const Real minAbs = geomScaling.GetLocal(jLoc,0);
        geomScaling.SetLocal(jLoc,0,Sqrt(minAbs*maxAbs));
    }
}

template<typename F>
void GeomEquil
( Matrix<F>& A, Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("GeomEquil"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbs( A );
    Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        cout << "Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    const Real sqrtDamp = Sqrt(damp);
    for( Int iter=0; iter<maxIter; ++iter )
    {
        // Geometrically equilibrate the columns
        for( Int j=0; j<n; ++j )
        {
            auto aCol = A( IR(0,m), IR(j,j+1) );
            auto maxColAbs = VectorMaxAbs( aCol );
            const Real maxColAbsVal = maxColAbs.value;
            const Real minColAbsVal = MinAbsNonzero( aCol, maxColAbsVal );
            if( maxColAbsVal > Real(0) )
            {
                const Real propScale = Sqrt(minColAbsVal*maxColAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxColAbsVal);
                Scale( Real(1)/scale, aCol );         
                dCol.Set( j, 0, scale*dCol.Get(j,0) );
            }
        }

        // Geometrically equilibrate the rows
        for( Int i=0; i<m; ++i )
        {
            auto aRow = A( IR(i,i+1), IR(0,n) );
            auto maxRowAbs = VectorMaxAbs( aRow );
            const Real maxRowAbsVal = maxRowAbs.value;
            const Real minRowAbsVal = MinAbsNonzero( aRow, maxRowAbsVal );
            if( maxRowAbsVal > Real(0) )
            {
                const Real propScale = Sqrt(minRowAbsVal*maxRowAbsVal);
                const Real scale = Max(propScale,sqrtDamp*maxRowAbsVal);
                Scale( Real(1)/scale, aRow );         
                dRow.Set( i, 0, scale*dRow.Get(i,0) );
            }
        }

        auto newMaxAbs = MaxAbs( A );
        Real newMaxAbsVal = newMaxAbs.value;
        Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "New ratio is " << newMaxAbsVal << "/" 
                 << newMinAbsVal << "=" << newRatio << endl;
        if( iter >= minIter && newRatio >= ratio*relTol )
            break;
        ratio = newRatio;
    }

    // Scale each column so that its maximum entry is 1 or 0
    for( Int j=0; j<n; ++j )
    {
        auto aCol = A( IR(0,m), IR(j,j+1) );
        auto maxColAbs = VectorMaxAbs( aCol );
        const Real maxColAbsVal = maxColAbs.value;
        if( maxColAbsVal > Real(0) )
        {
            Scale( Real(1)/maxColAbsVal, aCol );         
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
    DEBUG_ONLY(CallStackEntry cse("GeomEquil"))
    typedef Base<F> Real;

    auto APtr    = ReadWriteProxy<F,MC,MR>(&APre);     auto& A    = *APtr;
    auto dRowPtr = WriteProxy<Real,MC,STAR>(&dRowPre); auto& dRow = *dRowPtr;
    auto dColPtr = WriteProxy<Real,MR,STAR>(&dColPre); auto& dCol = *dColPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbs( A );
    Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && A.Grid().Rank() == 0 )
        cout << "Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
             << ratio << endl;

    const Real sqrtDamp = Sqrt(damp);
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
        Real newMaxAbsVal = newMaxAbs.value;
        Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && A.Grid().Rank() == 0 )
            cout << "New ratio is " << newMaxAbsVal << "/" 
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
void GeomEquil
( SparseMatrix<F>& A, Matrix<Base<F>>& dRow, Matrix<Base<F>>& dCol,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("GeomEquil"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    Ones( dRow, m, 1 );
    Ones( dCol, n, 1 );

    // TODO: Expose these as control parameters
    const Int minIter = 3;
    const Int maxIter = 10; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbs( A );
    Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress )
        cout << "Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
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
                rowScale.Set(j,0,scale);
                dCol.Set( j, 0, scale*dCol.Get(j,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    transValBuf[e] /= scale;
            }
            else
                rowScale.Set(j,0,Real(1));
        }
        DiagonalSolve( RIGHT, NORMAL, rowScale, A );

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
                colScale.Set(i,0,scale);
                dRow.Set( i, 0, scale*dRow.Get(i,0) );
                for( Int e=offset; e<offset+numConnect; ++e )
                    valBuf[e] /= scale;
            } 
            else
                colScale.Set(i,0,Real(1));
        }
        DiagonalSolve( RIGHT, NORMAL, colScale, ATrans );

        // Determine whether we are done or not
        // ------------------------------------
        auto newMaxAbs = MaxAbs( A );
        Real newMaxAbsVal = newMaxAbs.value;
        Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress )
            cout << "New ratio is " << newMaxAbsVal << "/" 
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
void GeomEquil
( DistSparseMatrix<F>& A, 
  DistMultiVec<Base<F>>& dRow, DistMultiVec<Base<F>>& dCol, 
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("GeomEquil"))
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
    const Int maxIter = 10; 
    const Real damp = Real(1)/Real(1000);
    const Real relTol = Real(9)/Real(10);

    // Compute the original ratio of the maximum to minimum nonzero
    auto maxAbs = MaxAbs( A );
    Real maxAbsVal = maxAbs.value;
    if( maxAbsVal == Real(0) )
        return;
    Real minAbsVal = MinAbsNonzero( A, maxAbsVal );
    Real ratio = maxAbsVal / minAbsVal;
    if( progress && commRank == 0 )
        cout << "Original ratio is " << maxAbsVal << "/" << minAbsVal << "="
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
        Real newMaxAbsVal = newMaxAbs.value;
        Real newMinAbsVal = MinAbsNonzero( A, newMaxAbsVal );
        Real newRatio = newMaxAbsVal / newMinAbsVal;
        if( progress && commRank == 0 )
            cout << "New ratio is " << newMaxAbsVal << "/" 
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
    bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
