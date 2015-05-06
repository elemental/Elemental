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
void MakeSymmetric( UpperOrLower uplo, Matrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CSE cse("MakeSymmetric"))
    const Int n = A.Width();
    if( A.Height() != n )
        LogicError("Cannot make non-square matrix symmetric");

    if( conjugate )
        MakeDiagonalReal(A);

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
    DEBUG_ONLY(CSE cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    MakeTrapezoidal( uplo, A );
    if( conjugate )
        MakeDiagonalReal(A);

    unique_ptr<AbstractDistMatrix<T>> 
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
    DEBUG_ONLY(CSE cse("MakeSymmetric"))
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
    A.ProcessQueues();
}

template<typename T>
void MakeSymmetric( UpperOrLower uplo, DistSparseMatrix<T>& A, bool conjugate )
{
    DEBUG_ONLY(CSE cse("MakeSymmetric"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make non-square matrix symmetric");

    MakeTrapezoidal( uplo, A );
    const Int numLocalEntries = A.NumLocalEntries();
    T* vBuf = A.ValueBuffer();
    const Int* sBuf = A.LockedSourceBuffer();
    const Int* tBuf = A.LockedTargetBuffer();
    if( conjugate && IsComplex<T>::val )
    {
        for( Int k=0; k<numLocalEntries; ++k )
            if( sBuf[k] == tBuf[k] )
                vBuf[k] = RealPart(vBuf[k]);
    }

    // Compute the number of entries to send to each process
    // =====================================================
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size(comm);
    vector<int> sendCounts(commSize,0);
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = sBuf[k];
        const Int j = tBuf[k];
        if( (uplo == LOWER && i > j) || (uplo == UPPER && i < j) )
            ++sendCounts[ A.RowOwner(j) ];
    }

    // Pack the triplets
    // =================
    vector<int> sendOffs;
    const int totalSend = Scan( sendCounts, sendOffs );
    vector<Entry<T>> sendBuf(totalSend);
    auto offs = sendOffs;
    for( Int k=0; k<numLocalEntries; ++k )
    {
        const Int i = sBuf[k];
        const Int j = tBuf[k];
        if( (uplo == LOWER && i > j) || (uplo == UPPER && i < j) )
        {
            const int owner = A.RowOwner(j);
            const Int s = offs[owner]++;
            sendBuf[s].i = j;
            sendBuf[s].j = i;
            sendBuf[s].value = ( conjugate ? Conj(vBuf[k]) : vBuf[k] );
        }
    }

    // Exchange and unpack the triplets
    // ================================
    auto recvBuf = mpi::AllToAll( sendBuf, sendCounts, sendOffs, comm );
    A.Reserve( A.NumLocalEntries()+recvBuf.size() );
    for( auto& entry : recvBuf )
        A.QueueUpdate( entry );
    A.ProcessQueues();
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

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
