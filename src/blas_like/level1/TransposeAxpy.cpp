/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T,typename S>
void TransposeAxpy( S alphaS, const Matrix<T>& X, Matrix<T>& Y, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("TransposeAxpy"))
    const T alpha = T(alphaS);
    const Int mX = X.Height();
    const Int nX = X.Width();
    const Int nY = Y.Width();
    const Int ldX = X.LDim();
    const Int ldY = Y.LDim();
    const T* XBuf = X.LockedBuffer();
          T* YBuf = Y.Buffer();
    // If X and Y are vectors, we can allow one to be a column and the other
    // to be a row. Otherwise we force X and Y to be the same dimension.
    if( mX == 1 || nX == 1 )
    {
        const Int lengthX = ( nX==1 ? mX : nX );
        const Int incX = ( nX==1 ? 1  : ldX );
        const Int incY = ( nY==1 ? 1  : ldY );
        DEBUG_ONLY(
            const Int mY = Y.Height();
            const Int lengthY = ( nY==1 ? mY : nY );
            if( lengthX != lengthY )
                LogicError("Nonconformal TransposeAxpy");
        )
        if( conjugate )
            for( Int j=0; j<lengthX; ++j )
                YBuf[j*incY] += alpha*Conj(XBuf[incX]);
        else
            blas::Axpy( lengthX, alpha, XBuf, incX, YBuf, incY );
    }
    else
    {
        DEBUG_ONLY(
            const Int mY = Y.Height();
            if( mX != nY || nX != mY )
                LogicError("Nonconformal TransposeAxpy");
        )
        if( nX <= mX )
        {
            if( conjugate )
                for( Int j=0; j<nX; ++j )
                    for( Int i=0; i<mX; ++i )
                        YBuf[j+i*ldY] += alpha*Conj(XBuf[i+j*ldX]);
            else
                for( Int j=0; j<nX; ++j )
                    blas::Axpy( mX, alpha, &XBuf[j*ldX], 1, &YBuf[j], ldY );
        }
        else
        {
            if( conjugate )
                for( Int i=0; i<mX; ++i )
                    for( Int j=0; j<nX; ++j )
                        YBuf[j+i*ldY] += alpha*Conj(XBuf[i+j*ldX]);
            else
                for( Int i=0; i<mX; ++i )
                    blas::Axpy( nX, alpha, &XBuf[i], ldX, &YBuf[i*ldY], 1 );
        }
    }
}

template<typename T,typename S>
void TransposeAxpy
( S alphaS, const SparseMatrix<T>& X, SparseMatrix<T>& Y, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("TransposeAxpy"))
    if( X.Height() != Y.Width() || X.Width() != Y.Height() )
        LogicError("X and Y must have transposed dimensions");
    const T alpha = T(alphaS);
    const Int numEntries = X.NumEntries();
    Y.Reserve( Y.NumEntries()+numEntries );
    if( conjugate )
        for( Int k=0; k<numEntries; ++k ) 
            Y.QueueUpdate( X.Col(k), X.Row(k), alpha*Conj(X.Value(k)) );
    else
        for( Int k=0; k<numEntries; ++k ) 
            Y.QueueUpdate( X.Col(k), X.Row(k), alpha*X.Value(k) );
    Y.MakeConsistent();
}

template<typename T,typename S>
void TransposeAxpy
( S alphaS, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B,
  bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("TransposeAxpy");
        AssertSameGrids( A, B );
        if( A.Height() != B.Width() || A.Width() != B.Height() )
            LogicError("A and B must have transposed dimensions");
    )
    const T alpha = T(alphaS);

    const DistData ADistData = A.DistData();
    const DistData BDistData = B.DistData();
    if( ADistData.colDist == BDistData.rowDist &&
        ADistData.rowDist == BDistData.colDist &&
        ADistData.colAlign==BDistData.rowAlign &&
        ADistData.rowAlign==BDistData.colAlign )
    {
        TransposeAxpy( alpha, A.LockedMatrix(), B.Matrix(), conjugate );
    }
    else
    {
        std::unique_ptr<AbstractDistMatrix<T>>
            C( B.ConstructTranspose(A.Grid(),A.Root()) );
        C->AlignRowsWith( B.DistData() );
        C->AlignColsWith( B.DistData() );
        Copy( A, *C );
        TransposeAxpy( alpha, C->LockedMatrix(), B.Matrix(), conjugate );
    }
}

template<typename T,typename S>
void TransposeAxpy
( S alphaS, const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B,
  bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("TransposeAxpy"))
    if( A.Height() != B.Width() || A.Width() != B.Height() )
        LogicError("A and B must have transposed dimensions");
    if( A.Comm() != B.Comm() )
        LogicError("A and B must have the same communicator");

    // Compute the number of entries of A to send to each process
    // ==========================================================
    mpi::Comm comm = A.Comm();
    const Int commSize = mpi::Size( comm );
    std::vector<int> sendCounts(commSize,0);
    for( Int k=0; k<A.NumLocalEntries(); ++k )
        ++sendCounts[ B.RowOwner(A.Col(k)) ];
    std::vector<int> recvCounts(commSize);
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, comm );

    // Convert the send/recv counts into offsets and total sizes
    // =========================================================
    std::vector<int> sendOffsets, recvOffsets;
    const int totalSend = Scan( sendCounts, sendOffsets );
    const int totalRecv = Scan( recvCounts, recvOffsets );

    // Pack the triplets
    // =================
    std::vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
    std::vector<T> vSendBuf(totalSend);
    auto offsets = sendOffsets;
    for( Int k=0; k<A.NumLocalEntries(); ++k )
    {
        const Int j = A.Col(k);
        const int owner = B.RowOwner(j);
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
    B.Reserve( B.NumLocalEntries()+totalRecv );
    for( Int k=0; k<totalRecv; ++k )
        B.QueueLocalUpdate
        ( sRecvBuf[k]-B.FirstLocalRow(), tRecvBuf[k], vRecvBuf[k] );
    B.MakeConsistent();
}

#define PROTO_TYPES(T,S) \
  template void TransposeAxpy \
  ( S alpha, const Matrix<T>& A, Matrix<T>& B, bool conjugate ); \
  template void TransposeAxpy \
  ( S alpha, const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B, \
    bool conjugate ); \
  template void TransposeAxpy \
  ( S alpha, const SparseMatrix<T>& A, SparseMatrix<T>& B, bool conjugate ); \
  template void TransposeAxpy \
  ( S alpha, const DistSparseMatrix<T>& A, DistSparseMatrix<T>& B, \
    bool conjugate );

#define PROTO_INT(T) PROTO_TYPES(T,T)

#define PROTO_REAL(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,T)

#define PROTO_COMPLEX(T) \
  PROTO_TYPES(T,Int) \
  PROTO_TYPES(T,Base<T>) \
  PROTO_TYPES(T,T)

#include "El/macros/Instantiate.h"

} // namespace El
