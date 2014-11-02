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
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    const Int n = A.Width();
    if( A.Height() != n )
        LogicError("Cannot make non-square matrix symmetric");

    if( conjugate )
        A.MakeDiagonalReal();

    T* ABuf = A.Buffer();
    const Int ldim = A.LDim();
    if( uplo == LOWER )
    {
        for( Int j=0; j<n; ++j )
        {
            for( Int i=j+1; i<n; ++i )
            {
                if( conjugate )
                    ABuf[j+i*ldim] = Conj(ABuf[i+j*ldim]); 
                else
                    ABuf[j+i*ldim] = ABuf[i+j*ldim];
            }    
        }
    }
    else
    {
        for( Int j=0; j<n; ++j )
        {
            for( Int i=0; i<j; ++i )
            {
                if( conjugate )
                    ABuf[j+i*ldim] = Conj(ABuf[i+j*ldim]);
                else
                    ABuf[j+i*ldim] = ABuf[i+j*ldim];
            }
        }
    }
}

template<typename T>
void MakeSymmetric
( UpperOrLower uplo, AbstractDistMatrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    MakeTrapezoidal( uplo, A );
    if( conjugate )
        A.MakeDiagonalReal();

    std::unique_ptr<AbstractDistMatrix<T>> 
      ATrans( A.Construct(A.Grid(),A.Root()) );
    Transpose( A, *ATrans, conjugate );
    if( uplo == LOWER )
        AxpyTrapezoid( UPPER, T(1), *ATrans, A, 1 );
    else
        AxpyTrapezoid( LOWER, T(1), *ATrans, A, -1 );
}

template<typename T>
void MakeSymmetric( UpperOrLower uplo, SparseMatrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    MakeTrapezoidal( uplo, A );

    const Int numEntries = A.NumEntries();
    const Int* sBuf = A.LockedSourceBuffer();
    const Int* tBuf = A.LockedTargetBuffer();
    T* vBuf = A.ValueBuffer();
    if( conjugate && IsComplex<T>::val )
        for( Int k=0; k<numEntries; ++k )
            if( sBuf[k] == tBuf[k] ) 
                vBuf[k] = RealPart(vBuf[k]);

    if( uplo == LOWER )
    {
        // Transpose the strictly lower triangle onto the upper triangle
        for( Int k=0; k<numEntries; ++k ) 
        {
            if( sBuf[k] > tBuf[k] )
            {
                if( conjugate )
                    A.QueueUpdate( tBuf[k], sBuf[k], Conj(vBuf[k]) );
                else
                    A.QueueUpdate( tBuf[k], sBuf[k], vBuf[k] );
            }
        }
    }
    else
    {
        // Transpose the strictly upper triangle onto the lower triangle
        for( Int k=0; k<numEntries; ++k ) 
        {
            if( sBuf[k] < tBuf[k] )
            {
                if( conjugate )
                    A.QueueUpdate( tBuf[k], sBuf[k], Conj(vBuf[k]) );
                else
                    A.QueueUpdate( tBuf[k], sBuf[k], vBuf[k] );
            }
        }
    }
    A.MakeConsistent();
}

template<typename T>
void MakeSymmetric( UpperOrLower uplo, DistSparseMatrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    MakeTrapezoidal( uplo, A );
    if( conjugate && IsComplex<T>::val )
    {
        const Int numLocalEntries = A.NumLocalEntries();
        T* vBuf = A.ValueBuffer();
        const Int* sBuf = A.LockedSourceBuffer();
        const Int* tBuf = A.LockedTargetBuffer();
        for( Int k=0; k<numLocalEntries; ++k )
            if( sBuf[k] == tBuf[k] )
                vBuf[k] = RealPart(vBuf[k]);
    }

    // Compute the number of entries to send to each process
    // =====================================================
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size(comm);
    std::vector<int> sendCounts(commSize,0);
    const Int firstLocalRow = A.FirstLocalRow();
    const Int localHeight = A.LocalHeight();
    const Int blocksize = A.Blocksize();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = firstLocalRow + iLoc;
        const Int offset = A.EntryOffset( iLoc );
        const Int rowSize = A.NumConnections( iLoc );
        for( Int k=0; k<rowSize; ++k )
        {
            const Int j = A.Col(offset+k);
            if( (uplo == LOWER && i < j) || (uplo == UPPER && i > j) )
            { 
                const Int owner = RowToProcess( j, blocksize, commSize );
                ++sendCounts[owner];
            }
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
            if( (uplo == LOWER && i > j) || (uplo == UPPER && i < j) )
            {
                const Int owner = RowToProcess( j, blocksize, commSize );
                const Int s = offsets[owner];
                sourceBuf[s] = j;
                targetBuf[s] = i;
                valueBuf[s] = A.Value(offset+k);
                ++offsets[owner];
            }
        }
    }

    // Exchange and unpack the triplets
    // ================================
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
    A.Reserve( A.NumLocalEntries()+totalRecv );
    for( Int k=0; k<totalRecv; ++k )
        A.QueueLocalUpdate
        ( sourceRecvBuf[k]-firstLocalRow,
          targetRecvBuf[k], valueRecvBuf[k] );
    A.MakeConsistent();
}

#define PROTO(T) \
  template void MakeSymmetric \
  ( UpperOrLower uplo, Matrix<T>& A, bool conjugate ); \
  template void MakeSymmetric \
  ( UpperOrLower uplo, AbstractDistMatrix<T>& A, bool conjugate ); \
  template void MakeSymmetric \
  ( UpperOrLower uplo, SparseMatrix<T>& A, bool conjugate ); \
  template void MakeSymmetric \
  ( UpperOrLower uplo, DistSparseMatrix<T>& A, bool conjugate );

#include "El/macros/Instantiate.h"

} // namespace El
