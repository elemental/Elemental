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
      CSE cse("Syrk");
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
    DEBUG_ONLY(CSE cse("Syrk"))
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
    DEBUG_ONLY(CSE cse("Syrk"))
    ScaleTrapezoid( beta, uplo, C );
    if( uplo == LOWER && orientation == NORMAL )
        syrk::LN( alpha, A, C, conjugate );
    else if( uplo == LOWER )
        syrk::LT( alpha, A, C, conjugate );
    else if( orientation == NORMAL )
        syrk::UN( alpha, A, C, conjugate );
    else
        syrk::UT( alpha, A, C, conjugate );
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const AbstractDistMatrix<T>& A, 
                 AbstractDistMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CSE cse("Syrk"))
    const Int n = ( orientation==NORMAL ? A.Height() : A.Width() );
    Zeros( C, n, n );
    Syrk( uplo, orientation, alpha, A, T(0), C, conjugate );
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const SparseMatrix<T>& A, 
  T beta,        SparseMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CSE cse("Syrk"))

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
        const Int offset = A.RowOffset(k);
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
    C.ProcessQueues();
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const SparseMatrix<T>& A, 
                 SparseMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CSE cse("Syrk"))
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
  T beta,        DistSparseMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CSE cse("Syrk"))

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
    vector<int> sendSizes(commSize,0);
    const Int localHeightA = A.LocalHeight();
    for( Int kLoc=0; kLoc<localHeightA; ++kLoc )
    {
        const Int offset = A.RowOffset(kLoc);
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

    // Pack the send buffers
    // ===================== 
    vector<int> sendOffs;
    const int totalSend = Scan( sendSizes, sendOffs );
    vector<Entry<T>> sendBuf(totalSend);
    auto offs = sendOffs;
    for( Int kLoc=0; kLoc<localHeightA; ++kLoc )
    {
        const Int offset = A.RowOffset(kLoc);
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
                    const Int s = offs[owner]++;
                    sendBuf[s].i = i;
                    sendBuf[s].j = j;
                    if( conjugate )
                        sendBuf[s].value = T(alpha)*Conj(A_ki)*A_kj;
                    else
                        sendBuf[s].value = T(alpha)*A_ki*A_kj;
                }
            }
        }
    }

    // Receive and apply the updates
    // =============================
    auto recvBuf = mpi::AllToAll( sendBuf, sendSizes, sendOffs, comm );
    C.Reserve( C.NumLocalEntries() + recvBuf.size() );
    for( auto& entry : recvBuf )
        C.QueueUpdate( entry );
    C.ProcessQueues();
}

template<typename T>
void Syrk
( UpperOrLower uplo, Orientation orientation,
  T alpha, const DistSparseMatrix<T>& A, 
                 DistSparseMatrix<T>& C, bool conjugate )
{
    DEBUG_ONLY(CSE cse("Syrk"))
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
