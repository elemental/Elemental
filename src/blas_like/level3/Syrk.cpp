/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Syrk/LN.hpp"
#include "./Syrk/LT.hpp"
#include "./Syrk/UN.hpp"
#include "./Syrk/UT.hpp"

namespace El {

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, T beta, Matrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(
        CallStackEntry cse("Syrk");
        if( orientation == NORMAL )
        {
            if( A.Height() != C.Height() || A.Height() != C.Width() )
                LogicError("Nonconformal Syrk");
        }
        else
        {
            if( A.Width() != C.Height() || A.Width() != C.Width() )
                LogicError("Nonconformal Syrk");
        }
    )
    const char uploChar = UpperOrLowerToChar( uplo );
    const char transChar = OrientationToChar( orientation );
    const Int k = ( orientation == NORMAL ? A.Width() : A.Height() );
    if( conjugate )
    {
        blas::Herk
        ( uploChar, transChar, C.Height(), k,
          RealPart(alpha), A.LockedBuffer(), A.LDim(),
          RealPart(beta),  C.Buffer(),       C.LDim() );
    }
    else
    {
        blas::Syrk
        ( uploChar, transChar, C.Height(), k,
          alpha, A.LockedBuffer(), A.LDim(),
          beta,  C.Buffer(),       C.LDim() );
    }
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const Matrix<T>& A, Matrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, alpha, A, T(0), C, conjugate );
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, 
  T beta,        AbstractDistMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))
    if( uplo == LOWER && orientation == NORMAL )
        syrk::LN( alpha, A, beta, C, conjugate );
    else if( uplo == LOWER )
        syrk::LT( alpha, A, beta, C, conjugate );
    else if( orientation == NORMAL )
        syrk::UN( alpha, A, beta, C, conjugate );
    else
        syrk::UT( alpha, A, beta, C, conjugate );
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, 
                 AbstractDistMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, alpha, A, T(0), C, conjugate );
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const SparseMatrix<T>& A, 
  T beta, SparseMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))

    if( orientation == NORMAL )
    {
        SparseMatrix<T> B;
        Transpose( A, B, conjugate );
        const Orientation newOrient = ( conjugate ? ADJOINT : TRANSPOSE );
        Syrk( uplo, newOrient, alpha, B, beta, C, conjugate );
        return;
    }

    const Int m = A.Height();
    const Int n = A.Width();
    if( C.Height() != n || C.Width() != n )
        LogicError("C was of the incorrect size");

    ScaleTrapezoid( beta, uplo, C );

    // Compute an upper bound on the required capacity
    // ===============================================
    Int newCapacity = C.NumEntries();
    for( Int k=0; k<m; ++k )
    {
        const Int numConn = A.NumConnections(k);
        newCapacity += numConn*numConn;
    }
    C.Reserve( newCapacity );

    // Queue the updates
    // =================
    for( Int k=0; k<m; ++k )
    {
        const Int offset = A.EntryOffset(k);
        const Int numConn = A.NumConnections(k);
        for( Int iConn=0; iConn<numConn; ++iConn )
        {
            const Int i = A.Col(offset+iConn);
            const T A_ki = A.Value(offset+iConn);
            for( Int jConn=0; jConn<numConn; ++jConn )
            {
                const Int j = A.Col(offset+jConn);
                if( (uplo == LOWER && i >= j) || (uplo == UPPER && i <= j) )
                {
                    const T A_kj = A.Value(offset+jConn);
                    if( conjugate )
                        C.QueueUpdate( i, j, T(alpha)*Conj(A_ki)*A_kj ); 
                    else
                        C.QueueUpdate( i, j, T(alpha)*A_ki*A_kj );
                }
            }
        }
    }

    // Force the result to be consistent
    // =================================
    C.MakeConsistent();
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const SparseMatrix<T>& A, 
                 SparseMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( orientation == NORMAL )
        Zeros( C, m, m );
    else
        Zeros( C, n, n );
    Syrk( uplo, orientation, alpha, A, T(0), C, conjugate );
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistSparseMatrix<T>& A, 
  T beta, DistSparseMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))

    if( orientation == NORMAL )
    {
        DistSparseMatrix<T> B(A.Comm());
        Transpose( A, B, conjugate );
        const Orientation newOrient = ( conjugate ? ADJOINT : TRANSPOSE );
        Syrk( uplo, newOrient, alpha, B, beta, C, conjugate );
        return;
    }

    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const Int commSize = mpi::Size( comm );

    if( C.Height() != n || C.Width() != n )
        LogicError("C was of the incorrect size");
    if( C.Comm() != comm )
        LogicError("Communicators of A and C must match");

    ScaleTrapezoid( beta, uplo, C );

    // Count the number of entries that we will send to each process
    // =============================================================
    std::vector<int> sendSizes(commSize,0);
    const Int localHeightA = A.LocalHeight();
    for( Int kLoc=0; kLoc<localHeightA; ++kLoc )
    {
        const Int offset = A.EntryOffset(kLoc);
        const Int numConn = A.NumConnections(kLoc);
        for( Int iConn=0; iConn<numConn; ++iConn )
        {
            const Int i = A.Col(offset+iConn);
            const int owner = C.RowOwner(i);
            for( Int jConn=0; jConn<numConn; ++jConn )
            {
                const Int j = A.Col(offset+jConn);
                if( (uplo==LOWER && i>=j) || (uplo==UPPER && i<=j) )
                    ++sendSizes[owner];
            }
        }
    }
    std::vector<int> recvSizes(commSize);
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm );

    // Convert the send and recv counts to offsets and total sizes
    // ===========================================================
    std::vector<int> sendOffsets, recvOffsets;
    const int totalSend = Scan( sendSizes, sendOffsets );
    const int totalRecv = Scan( recvSizes, recvOffsets );

    // Pack the send buffers
    // ===================== 
    std::vector<Int> sourceBuf(totalSend), targetBuf(totalSend);
    std::vector<T> valueBuf(totalSend);
    auto offsets = sendOffsets;
    for( Int kLoc=0; kLoc<localHeightA; ++kLoc )
    {
        const Int offset = A.EntryOffset(kLoc);
        const Int numConn = A.NumConnections(kLoc);
        for( Int iConn=0; iConn<numConn; ++iConn ) 
        {
            const Int i = A.Col(offset+iConn);
            const T A_ki = A.Value(offset+iConn);
            const int owner = C.RowOwner(i);
            for( Int jConn=0; jConn<numConn; ++jConn )
            {
                const Int j = A.Col(offset+jConn);
                if( (uplo==LOWER && i>=j) || (uplo==UPPER && i<=j) )
                {
                    const T A_kj = A.Value(offset+jConn);
                    const Int s = offsets[owner]++;
                    sourceBuf[s] = i;
                    targetBuf[s] = j;
                    if( conjugate )
                        valueBuf[s] = T(alpha)*Conj(A_ki)*A_kj;
                    else
                        valueBuf[s] = T(alpha)*A_ki*A_kj;
                }
            }
        }
    }

    // Receive the updates
    // ===================
    const Int oldSize = C.NumLocalEntries();
    const Int newCapacity = oldSize + totalRecv;
    std::vector<Int> sourceRecvBuf(totalRecv), targetRecvBuf(totalRecv);
    std::vector<T> valueRecvBuf(totalRecv);
    mpi::AllToAll
    ( sourceBuf.data(), sendSizes.data(), sendOffsets.data(),
      sourceRecvBuf.data(), recvSizes.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( targetBuf.data(), sendSizes.data(), sendOffsets.data(),
      targetRecvBuf.data(), recvSizes.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( valueBuf.data(), sendSizes.data(), sendOffsets.data(),
      valueRecvBuf.data(), recvSizes.data(), recvOffsets.data(), comm );
    C.Reserve( newCapacity );
    for( Int k=0; k<totalRecv; ++k )
        C.QueueLocalUpdate
        ( sourceRecvBuf[k]-C.FirstLocalRow(), 
          targetRecvBuf[k], valueRecvBuf[k] );
    // NOTE: In order to avoid unnecessary workspace, we can directly manipulate
    //       the internals of the class
    /*
    C.distGraph_.consistent_ = false;
    C.distGraph_.sources_.resize( newCapacity );
    C.distGraph_.targets_.resize( newCapacity );
    C.vals_.resize( newCapacity );
    mpi::AllToAll
    ( sourceBuf.data(), sendSizes.data(), sendOffsets.data(),
      C.distGraph_.sources_.data()+oldSize, 
      recvSizes.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( targetBuf.data(), sendSizes.data(), sendOffsets.data(),
      C.distGraph_.targets_.data()+oldSize, 
      recvSizes.data(), recvOffsets.data(), comm );
    mpi::AllToAll
    ( valueBuf.data(), sendSizes.data(), sendOffsets.data(),
      C.vals_.data()+oldSize, 
      recvSizes.data(), recvOffsets.data(), comm );
    */

    // Make the distributed matrix consistent
    // ======================================
    C.MakeConsistent();
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistSparseMatrix<T>& A, 
                 DistSparseMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CallStackEntry cse("Syrk"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( orientation == NORMAL )
        Zeros( C, m, m );
    else
        Zeros( C, n, n );
    Syrk( uplo, orientation, alpha, A, T(0), C, conjugate );
}

#define PROTO(T) \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const Matrix<T>& A, T beta, Matrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const Matrix<T>& A, Matrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const AbstractDistMatrix<T>& A, \
    T beta, AbstractDistMatrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const AbstractDistMatrix<T>& A, \
                   AbstractDistMatrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const SparseMatrix<T>& A, \
    T beta,        SparseMatrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const SparseMatrix<T>& A, \
                   SparseMatrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const DistSparseMatrix<T>& A, \
    T beta,        DistSparseMatrix<T>& C, bool conjugate ); \
  template void Syrk \
  ( UpperOrLower uplo, Orientation orientation, \
    T alpha, const DistSparseMatrix<T>& A, \
                   DistSparseMatrix<T>& C, bool conjugate );

// blas::Syrk is not yet supported for Int
#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
