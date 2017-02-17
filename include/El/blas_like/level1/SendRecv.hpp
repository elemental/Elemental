/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_SENDRECV_HPP
#define EL_BLAS_SENDRECV_HPP

namespace El {

template<typename T>
void SendRecv
( const Matrix<T>& A, Matrix<T>& B, mpi::Comm comm, int sendRank, int recvRank )
{
    EL_DEBUG_CSE
    const Int heightA = A.Height();
    const Int heightB = B.Height();
    const Int widthA = A.Width();
    const Int widthB = B.Width();
    const Int sizeA = heightA*widthA;
    const Int sizeB = heightB*widthB;
    if( heightA == A.LDim() && heightB == B.LDim() )
    {
        mpi::SendRecv
        ( A.LockedBuffer(), sizeA, sendRank,
          B.Buffer(),       sizeB, recvRank, comm );
    }
    else if( heightA == A.LDim() )
    {
        vector<T> recvBuf;
        FastResize( recvBuf, sizeB );
        mpi::SendRecv
        ( A.LockedBuffer(), sizeA, sendRank,
          recvBuf.data(),   sizeB, recvRank, comm );
        copy::util::InterleaveMatrix
        ( heightB, widthB,
          recvBuf.data(), 1, heightB,
          B.Buffer(),     1, B.LDim() );
    }
    else
    {
        vector<T> sendBuf;
        FastResize( sendBuf, sizeA );
        copy::util::InterleaveMatrix
        ( heightA, widthA,
          A.LockedBuffer(), 1, A.LDim(),
          sendBuf.data(),   1, heightA );

        vector<T> recvBuf;
        FastResize( recvBuf, sizeB );
        mpi::SendRecv
        ( sendBuf.data(), sizeA, sendRank,
          recvBuf.data(), sizeB, recvRank, comm );
        copy::util::InterleaveMatrix
        ( heightB, widthB,
          recvBuf.data(), 1, heightB,
          B.Buffer(),     1, B.LDim() );
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void SendRecv \
  ( const Matrix<T>& A, Matrix<T>& B, mpi::Comm comm, \
    int sendRank, int recvRank );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_SENDRECV_HPP
