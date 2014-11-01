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
    B.Empty();
    B.Resize( n, m );
    B.Reserve( A.NumEntries() );
    for( Int i=0; i<m; ++i )
    {
        const Int offset = A.EntryOffset( i ); 
        const Int rowSize = A.NumConnections( i );
        for( Int k=0; k<rowSize; ++k )
        {
            const Int j = A.Col(offset+k);
            const T AVal = A.Value(offset+k);
            if( conjugate )
                B.QueueUpdate( j, i, Conj(AVal) );
            else
                B.QueueUpdate( j, i, AVal );
        }
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

    B.Empty();
    B.SetComm( A.Comm() );
    B.Resize( n, m );

    // Compute the number of entries of A to send to each process
    // ==========================================================
    std::vector<int> sendCounts(commSize,0);
    const Int firstLocalRow = A.FirstLocalRow();
    const Int localHeight = A.LocalHeight();
    const Int blocksizeB = B.Blocksize();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int offset = A.EntryOffset( iLoc );
        const Int rowSize = A.NumConnections( iLoc );
        for( Int k=0; k<rowSize; ++k )
        {
            const Int j = A.Col(offset+k);
            const Int owner = RowToProcess( j, blocksizeB, commSize );
            ++sendCounts[owner];
        }
    }

    // Communicate to determine the number we receive from each process
    // ================================================================
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
    std::vector<Int> sourceBuf(totalSend), targetBuf(totalSend);
    std::vector<T> valueBuf(totalSend);
    std::vector<int> offsets = sendOffsets;
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstLocalRow + iLoc;
        const Int rowSize = A.NumConnections( iLoc );
        const Int offset = A.EntryOffset( iLoc );
        for( Int k=0; k<rowSize; ++k )
        {
            const Int j = A.Col(offset+k);
            const Int owner = RowToProcess( j, blocksizeB, commSize );
            const Int s = offsets[owner];
            sourceBuf[s] = j;
            targetBuf[s] = i;
            valueBuf[s] = A.Value(offset+k); 
            ++offsets[owner];
        }
    }

    // Exchange and unpack the triplets
    // ================================
    // TODO: Switch to a mechanism which directly unpacks into the 
    //       class's local storage?
    std::vector<Int> sourceRecvBuf(totalRecv), targetRecvBuf(totalRecv);
    std::vector<T> valueRecvBuf(totalRecv);
    mpi::AllToAll
    ( sourceBuf.data(), sendCounts.data(), sendOffsets.data(),
      sourceRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( targetBuf.data(), sendCounts.data(), sendOffsets.data(),
      targetRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( valueBuf.data(), sendCounts.data(), sendOffsets.data(),
      valueRecvBuf.data(), recvCounts.data(), recvOffsets.data(), comm );
    B.Reserve( totalRecv );
    for( Int k=0; k<totalRecv; ++k )
        B.QueueLocalUpdate
        ( sourceRecvBuf[k]-B.FirstLocalRow(), 
          targetRecvBuf[k], valueRecvBuf[k] );
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
