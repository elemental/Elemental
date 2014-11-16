/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const Int m = A.Height();
    const Int n = A.Width();
    B.Resize( n, m );
    if( conjugate )
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                B.Set(j,i,Conj(A.Get(i,j)));
    }
    else
    {
        for( Int j=0; j<n; ++j )
            for( Int i=0; i<m; ++i )
                B.Set(j,i,A.Get(i,j));
    }
}

template<typename T>
void Transpose
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const DistData ADistData = A.DistData();
    const DistData BDistData = B.DistData();
    if( ADistData.colDist == BDistData.rowDist &&
        ADistData.rowDist == BDistData.colDist &&
        ((ADistData.colAlign==BDistData.rowAlign) || !B.RowConstrained()) &&
        ((ADistData.rowAlign==BDistData.colAlign) || !B.ColConstrained()) )
    {
        B.Align( A.RowAlign(), A.ColAlign() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        std::unique_ptr<AbstractDistMatrix<T>> 
            C( B.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignRowsWith( B.DistData() );
        C->AlignColsWith( B.DistData() );
        Copy( A, *C );
        B.Resize( A.Width(), A.Height() );
        Transpose( C->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void Transpose
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B, 
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const BlockDistData ADistData = A.DistData();
    const BlockDistData BDistData = B.DistData();
    if( ADistData.colDist == BDistData.rowDist &&
        ADistData.rowDist == BDistData.colDist &&
        ((ADistData.colAlign    == BDistData.rowAlign && 
          ADistData.blockHeight == BDistData.blockWidth &&
          ADistData.colCut      == BDistData.rowCut) || !B.RowConstrained()) &&
        ((ADistData.rowAlign   == BDistData.colAlign && 
          ADistData.blockWidth == BDistData.blockHeight &&
          ADistData.rowCut     == BDistData.colCut) || !B.ColConstrained()))
    {
        B.Align
        ( A.BlockWidth(), A.BlockHeight(), 
          A.RowAlign(), A.ColAlign(), A.RowCut(), A.ColCut() );
        B.Resize( A.Width(), A.Height() );
        Transpose( A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        std::unique_ptr<AbstractBlockDistMatrix<T>> 
            C( B.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignRowsWith( B.DistData() );
        C->AlignColsWith( B.DistData() );
        Copy( A, *C );
        B.Resize( A.Width(), A.Height() );
        Transpose( C->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T>
void Transpose
( const SparseMatrix<T>& A, SparseMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    Zeros( B, n, m );
    B.Reserve( numEntries );
    for( Int k=0; k<numEntries; ++k )
    {
        const Int i = A.Row(k);
        const Int j = A.Col(k);
        const T alpha = A.Value(k);
        if( conjugate )
            B.QueueUpdate( j, i, Conj(alpha) );
        else
            B.QueueUpdate( j, i, alpha );
    }
    B.MakeConsistent();
}

template<typename T>
void Transpose
( const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Transpose"))
    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const Int commSize = mpi::Size( comm );

    B.SetComm( comm );
    Zeros( B, n, m );

    // Compute the number of entries of A to send to each process
    // ==========================================================
    std::vector<int> sendCounts(commSize,0);
    for( Int k=0; k<A.NumLocalEntries(); ++k )
        ++sendCounts[ B.RowOwner(A.Col(k)) ];
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );

    // Convert the send/recv counts into offsets and total sizes
    // =========================================================
    Int totalSend=0, totalRecv=0;
    std::vector<int> sendOffsets(commSize), recvOffsets(commSize);
    for( Int q=0; q<commSize; ++q )
    {
        sendOffsets[q] = totalSend;
        recvOffsets[q] = totalRecv;
        totalSend += sendCounts[q];
        totalRecv += recvCounts[q];
    }

    // Pack the triplets
    // =================
    std::vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
    std::vector<T> vSendBuf(totalSend);
    std::vector<int> offsets = sendOffsets;
    for( Int k=0; k<A.NumLocalEntries(); ++k )
    {
        const Int j = A.Col(k);
        const Int owner = B.RowOwner(j);
        const Int s = offsets[owner];
        sSendBuf[s] = j;
        tSendBuf[s] = A.Row(k);
        if( conjugate )
            vSendBuf[s] = Conj(A.Value(k));
        else
            vSendBuf[s] = A.Value(k);
        ++offsets[owner];
    }

    // Exchange and unpack the triplets
    // ================================
    // TODO: Switch to a mechanism which directly unpacks into the 
    //       class's local storage?
    std::vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
    std::vector<T> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      sRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( tSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      tRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendCounts.data(), sendOffsets.data(),
      vRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    B.Reserve( totalRecv );
    for( Int k=0; k<totalRecv; ++k )
        B.QueueLocalUpdate
        ( sRecvBuf[k]-B.FirstLocalRow(), tRecvBuf[k], vRecvBuf[k] );
    B.MakeConsistent();
}

#define PROTO(T) \
  template void Transpose( const Matrix<T>& A, Matrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const AbstractDistMatrix<T>& A, \
          AbstractDistMatrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const AbstractBlockDistMatrix<T>& A, \
          AbstractBlockDistMatrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const SparseMatrix<T>& A, \
          SparseMatrix<T>& B, bool conjugate ); \
  template void Transpose \
  ( const DistSparseMatrix<T>& A, \
          DistSparseMatrix<T>& B, bool conjugate );

#include "El/macros/Instantiate.h"

} // namespace El
