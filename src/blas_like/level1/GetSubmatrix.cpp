/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Contiguous
// ==========
template<typename T>
void GetSubmatrix
( const Matrix<T>& A, Range<Int> I, Range<Int> J, 
        Matrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    auto ASubView = A(I,J);
    ASub = ASubView;
}

template<typename T>
Matrix<T> GetSubmatrix
( const Matrix<T>& A, Range<Int> I, Range<Int> J )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    Matrix<T> ASub;
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A, Range<Int> I, Range<Int> J, 
        AbstractDistMatrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    unique_ptr<AbstractDistMatrix<T>> 
      ASubView( A.Construct(A.Grid(),A.Root()) );
    LockedView( *ASubView, A, I, J );
    Copy( *ASubView, ASub );
}

template<typename T>
DistMatrix<T> GetSubmatrix
( const AbstractDistMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    DistMatrix<T> ASub( A.Grid() );
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J,
        SparseMatrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    if( I.end == END ) I.end = A.Height();
    if( J.end == END ) J.end = A.Width();
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    Zeros( ASub, mSub, nSub );

    // Reserve the number of nonzeros that live within the submatrix
    Int numNonzerosSub = 0;
    for( Int i=I.beg; i<I.end; ++i )
    {
        const Int rowOff = A.RowOffset(i);
        const Int numConn = A.NumConnections(i);
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = A.Col(e);
            if( j >= J.beg && j < J.end )
                ++numNonzerosSub;
        }
    }
    ASub.Reserve( numNonzerosSub );

    // Insert the nonzeros
    for( Int i=I.beg; i<I.end; ++i ) 
    {
        const Int rowOff = A.RowOffset(i);
        const Int numConn = A.NumConnections(i);
        for( Int e=rowOff; e<rowOff+numConn; ++e )
        {
            const Int j = A.Col(e);
            if( j >= J.beg && j < J.end )
                ASub.QueueUpdate( i-I.beg, j-J.beg, A.Value(e) );
        }
    }
    ASub.ProcessQueues();
}

template<typename T>
SparseMatrix<T> GetSubmatrix
( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    SparseMatrix<T> ASub;
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J,
        DistSparseMatrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    if( I.end == END ) I.end = A.Height();
    if( J.end == END ) J.end = A.Width();

    ASub.SetComm( A.Comm() );
    Zeros( ASub, I.end-I.beg, J.end-J.beg );

    // Count the number of updates
    // ===========================
    Int numUpdates = 0;
    for( Int e=0; e<A.NumLocalEntries(); ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
            ++numUpdates;
    }

    // Queue and process the updates
    // =============================
    ASub.Reserve( numUpdates, numUpdates );
    for( Int e=0; e<A.NumLocalEntries(); ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        if( i >= I.end )
            break;
        else if( i >= I.beg && j >= J.beg && j < J.end )
            ASub.QueueUpdate( i-I.beg, j-J.beg, A.Value(e), false );
    }
    ASub.ProcessQueues();
}

template<typename T>
DistSparseMatrix<T> GetSubmatrix
( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J )
{
    DistSparseMatrix<T> ASub(A.Comm());
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J,
        DistMultiVec<T>& ASub )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    if( I.end == END ) I.end = A.Height();
    if( J.end == END ) J.end = A.Width();
    const Int mSub = I.end-I.beg;
    const Int nSub = J.end-J.beg;
    const Int localHeight = A.LocalHeight();

    ASub.SetComm( A.Comm() );
    Zeros( ASub, mSub, nSub );
    
    // If no communication is necessary, take the easy and fast approach
    if( mSub == A.Height() )
    {
        for( Int j=J.beg; j<J.end; ++j )
            for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                ASub.SetLocal( iLoc, j-J.beg, A.GetLocal(iLoc,j) );
        return;
    }

    // Compute the number of updates
    // =============================
    Int numUpdates = 0;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
            numUpdates += nSub;
    }

    // Queue and process the updates
    // =============================
    ASub.Reserve( numUpdates );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = A.GlobalRow(iLoc);
        if( i >= I.end )
            break;
        else if( i >= I.beg )
        {
            for( Int j=J.beg; j<J.end; ++j )
                ASub.QueueUpdate( i-I.beg, j-J.beg, A.GetLocal(iLoc,j) );
        }
    }
    ASub.ProcessQueues();
}

template<typename T>
DistMultiVec<T> GetSubmatrix
( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J )
{
    DistMultiVec<T> ASub(A.Comm());
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

// Noncontiguous
// =============
template<typename T>
void GetSubmatrix
( const Matrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
        Matrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    const Int m = I.size();
    const Int n = J.size();

    Zeros( ASub, m, n );
    for( Int jSub=0; jSub<n; ++jSub )
    {
        const Int j = J[jSub];
        for( Int iSub=0; iSub<m; ++iSub )
        {
            const Int i = I[iSub];
            ASub.Set( iSub, jSub, A.Get(i,j) );
        }
    }
}

template<typename T>
Matrix<T> GetSubmatrix
( const Matrix<T>& A, const vector<Int>& I, const vector<Int>& J )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    Matrix<T> ASub;
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

template<typename T>
void GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J, 
        AbstractDistMatrix<T>& ASub )
{
    DEBUG_ONLY(CSE cse("GetSubmatrix"))
    const Grid& g = A.Grid();
    ASub.SetGrid( g ); 
    Zeros( ASub, I.size(), J.size() );

    // TODO: Intelligently pick the redundant rank to pack from?

    // Count the number of updates
    // ===========================
    Int numUpdates = 0;
    if( A.RedundantRank() == 0 )
        for( auto& i : I )
            if( A.IsLocalRow(i) )
                for( auto& j : J )
                    if( A.IsLocalCol(j) )
                        ++numUpdates;

    // Queue and process the updates
    // =============================
    ASub.Reserve( numUpdates );
    if( A.RedundantRank() == 0 )
    {
        for( size_t iSub=0; iSub<I.size(); ++iSub )
        {
            const Int i = I[iSub];
            if( A.IsLocalRow(i) )
            {
                const Int iLoc = A.LocalRow(i);
                for( size_t jSub=0; jSub<J.size(); ++jSub )
                {
                    const Int j = J[jSub];
                    if( A.IsLocalCol(j) )
                    {
                        const Int jLoc = A.LocalCol(j);
                        const T value = A.GetLocal(iLoc,jLoc);
                        ASub.QueueUpdate( iSub, jSub, value );
                    }
                }
            }
        }
    }
    ASub.ProcessQueues();
}

template<typename T>
DistMatrix<T> GetSubmatrix
( const AbstractDistMatrix<T>& A, 
  const vector<Int>& I, const vector<Int>& J )
{
    DistMatrix<T> ASub( A.Grid() );
    GetSubmatrix( A, I, J, ASub );
    return ASub;
}

#define PROTO(T) \
  /* Contiguous */ \
  template void GetSubmatrix \
  ( const Matrix<T>& A, Range<Int> I, Range<Int> J, \
          Matrix<T>& ASub ); \
  template Matrix<T> GetSubmatrix \
  ( const Matrix<T>& A, Range<Int> I, Range<Int> J ); \
  template void GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, Range<Int> I, Range<Int> J, \
          AbstractDistMatrix<T>& ASub ); \
  template DistMatrix<T> GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template void GetSubmatrix \
  ( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J, \
          SparseMatrix<T>& ASub ); \
  template SparseMatrix<T> GetSubmatrix \
  ( const SparseMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template void GetSubmatrix \
  ( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J, \
          DistSparseMatrix<T>& ASub ); \
  template DistSparseMatrix<T> GetSubmatrix \
  ( const DistSparseMatrix<T>& A, Range<Int> I, Range<Int> J ); \
  template void GetSubmatrix \
  ( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J, \
          DistMultiVec<T>& ASub ); \
  template DistMultiVec<T> GetSubmatrix \
  ( const DistMultiVec<T>& A, Range<Int> I, Range<Int> J ); \
  /* Noncontiguous */ \
  template void GetSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
          Matrix<T>& ASub ); \
  template Matrix<T> GetSubmatrix \
  ( const Matrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J ); \
  template void GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J, \
          AbstractDistMatrix<T>& ASub ); \
  template DistMatrix<T> GetSubmatrix \
  ( const AbstractDistMatrix<T>& A, \
    const vector<Int>& I, const vector<Int>& J ); \


#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
