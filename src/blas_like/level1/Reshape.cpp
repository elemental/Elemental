/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Reshape
( Int mNew, Int nNew, 
  const Matrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Reshape"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m*n != mNew*nNew )
        LogicError
        ("Reshape from ",m," x ",n," to ",mNew," x ",nNew,
         " did not preserve the total number of entries");

    Zeros( B, mNew, nNew );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<m; ++i )
        {
            const Int iNew = (i+j*m) % mNew;
            const Int jNew = (i+j*m) / mNew;
            B.Set( iNew, jNew, A.Get(i,j) );
        }
    }
}

template<typename T>
Matrix<T> Reshape( Int mNew, Int nNew, const Matrix<T>& A )
{
    DEBUG_ONLY(CSE cse("Reshape"))
    Matrix<T> B;
    Reshape( mNew, nNew, A, B );
    return B;
}

// TODO: Merge with implementation of GetSubmatrix via a function which maps
//       the coordinates in A to the coordinates in B
template<typename T>
void Reshape
( Int mNew, Int nNew, 
  const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Reshape"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const Grid& g = A.Grid();
    if( m*n != mNew*nNew )
        LogicError
        ("Reshape from ",m," x ",n," to ",mNew," x ",nNew,
         " did not preserve the total number of entries");

    B.SetGrid( g ); 
    Zeros( B, mNew, nNew );
    mpi::Comm comm = g.ViewingComm();
    const int commSize = mpi::Size( comm );  

    // TODO: Intelligently pick the redundant rank to pack from?

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);    
    if( A.RedundantRank() == 0 )
    {
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                const Int iNew = (i+j*m) % mNew;
                const Int jNew = (i+j*m) / mNew;
                const int owner = 
                  g.VCToViewing(
                    g.CoordsToVC
                    ( B.ColDist(), B.RowDist(), 
                      B.Owner(iNew,jNew), B.Root() ) 
                  );
                ++sendCounts[owner];
            }
        }
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    vector<Entry<T>> sendBuf(totalSend);
    if( A.RedundantRank() == 0 )
    {
        auto offs = sendOffs;
        for( Int jLoc=0; jLoc<nLocal; ++jLoc )
        {
            const Int j = A.GlobalCol(jLoc);
            for( Int iLoc=0; iLoc<mLocal; ++iLoc )
            {
                const Int i = A.GlobalRow(iLoc);
                const Int iNew = (i+j*m) % mNew;
                const Int jNew = (i+j*m) / mNew;
                const int owner = 
                  g.VCToViewing(
                    g.CoordsToVC
                    ( B.ColDist(), B.RowDist(), 
                      B.Owner(iNew,jNew), B.Root() ) 
                  );
                const T value = A.GetLocal(iLoc,jLoc);
                sendBuf[offs[owner]++] = 
                  Entry<T>{ Int(iNew), Int(jNew), value };
            }
        }
    }

    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll(sendBuf,sendCounts,sendOffs,comm);
    for( auto& entry : recvBuf )
        B.Set( entry );

    // Broadcast from the assigned member of the redundant rank team
    // =============================================================
    Broadcast( B, B.RedundantComm(), 0 );
}

template<typename T>
DistMatrix<T> Reshape( Int mNew, Int nNew, const AbstractDistMatrix<T>& A )
{
    DistMatrix<T> B( A.Grid() );
    Reshape( mNew, nNew, A, B );
    return B;
}

template<typename T>
void Reshape
( Int mNew, Int nNew, const SparseMatrix<T>& A, SparseMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Reshape"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    if( m*n != mNew*nNew )
        LogicError
        ("Reshape from ",m," x ",n," to ",mNew," x ",nNew,
         " did not preserve the total number of entries");

    Zeros( B, mNew, nNew );
    B.Reserve( numEntries );

    // Insert the nonzeros
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Int iNew = (i+j*m) % mNew;
        const Int jNew = (i+j*m) / mNew;
        B.QueueUpdate( iNew, jNew, A.Value(e) );
    }
    B.ProcessQueues();
}

template<typename T>
SparseMatrix<T> Reshape( Int mNew, Int nNew, const SparseMatrix<T>& A )
{
    SparseMatrix<T> B;
    Reshape( mNew, nNew, A, B );
    return B;
}

template<typename T>
void Reshape
( Int mNew, Int nNew, 
  const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B )
{
    DEBUG_ONLY(CSE cse("Reshape"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumLocalEntries();
    if( m*n != mNew*nNew )
        LogicError
        ("Reshape from ",m," x ",n," to ",mNew," x ",nNew,
         " did not preserve the total number of entries");

    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    B.SetComm( comm );
    Zeros( B, mNew, nNew );

    // Compute the metadata
    // ====================
    vector<int> sendCounts(commSize,0);
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const Int iNew = (i+j*m) % mNew;
        ++sendCounts[ B.RowOwner(iNew) ];
    }

    // Pack the data
    // =============
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    auto offs = sendOffs;
    vector<Entry<T>> sendBuf(totalSend);
    for( Int e=0; e<numEntries; ++e )
    {
        const Int i = A.Row(e);
        const Int j = A.Col(e);
        const T value = A.Value(e);
        const Int iNew = (i+j*m) % mNew;
        const Int jNew = (i+j*m) / mNew;
        const int owner = B.RowOwner( iNew );
        sendBuf[offs[owner]++] = Entry<T>{ iNew, jNew, value };
    }
    
    // Exchange and unpack the data
    // ============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    B.Reserve( recvBuf.size() );
    for( auto& entry : recvBuf )
        B.QueueUpdate( entry );
    B.ProcessQueues();
}

template<typename T>
DistSparseMatrix<T> Reshape( Int mNew, Int nNew, const DistSparseMatrix<T>& A )
{
    DistSparseMatrix<T> B(A.Comm());
    Reshape( mNew, nNew, A, B );
    return B;
}

#define PROTO(T) \
  template void Reshape \
  ( Int mNew, Int nNew, const Matrix<T>& A, Matrix<T>& B ); \
  template Matrix<T> Reshape( Int mNew, Int nNew, const Matrix<T>& A ); \
  template void Reshape \
  ( Int mNew, Int nNew, \
    const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template DistMatrix<T> Reshape \
  ( Int mNew, Int nNew, const AbstractDistMatrix<T>& A ); \
  template void Reshape \
  ( Int mNew, Int nNew, const SparseMatrix<T>& A, SparseMatrix<T>& B ); \
  template SparseMatrix<T> Reshape \
  ( Int mNew, Int nNew, const SparseMatrix<T>& A ); \
  template void Reshape \
  ( Int mNew, Int nNew, \
    const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B ); \
  template DistSparseMatrix<T> Reshape \
  ( Int mNew, Int nNew, const DistSparseMatrix<T>& A );

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
