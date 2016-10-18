/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_BROADCAST_HPP
#define EL_BLAS_BROADCAST_HPP

namespace El {

template<typename T>
void Broadcast( Matrix<T>& A, mpi::Comm comm, int rank )
{
    DEBUG_CSE
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commSize == 1 )
        return;

    const Int height = A.Height();
    const Int width = A.Width();
    const Int size = height*width;
    if( height == A.LDim() )
    {
        mpi::Broadcast( A.Buffer(), size, rank, comm );
    }
    else
    {
        vector<T> buf;
        FastResize( buf, size );

        // Pack
        if( commRank == rank )
            copy::util::InterleaveMatrix
            ( height, width,
              A.LockedBuffer(), 1, A.LDim(),
              buf.data(),       1, height );

        mpi::Broadcast( buf.data(), size, rank, comm );

        // Unpack
        if( commRank != rank )
            copy::util::InterleaveMatrix
            ( height,        width,
              buf.data(), 1, height,
              A.Buffer(), 1, A.LDim() );
    }
}

template<typename T>
void Broadcast( AbstractDistMatrix<T>& A, mpi::Comm comm, int rank )
{
    DEBUG_CSE
    const int commSize = mpi::Size( comm );
    const int commRank = mpi::Rank( comm );
    if( commSize == 1 )
        return;
    if( !A.Participating() )
        return;

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int localSize = localHeight*localWidth;
    if( localHeight == A.LDim() )
    {
        mpi::Broadcast( A.Buffer(), localSize, rank, comm );
    }
    else
    {
        vector<T> buf;
        FastResize( buf, localSize );

        // Pack
        if( commRank == rank )
            copy::util::InterleaveMatrix
            ( localHeight, localWidth,
              A.LockedBuffer(), 1, A.LDim(),
              buf.data(),       1, localHeight );

        mpi::Broadcast( buf.data(), localSize, rank, comm );

        // Unpack
        if( commRank != rank )
            copy::util::InterleaveMatrix
            ( localHeight, localWidth,
              buf.data(), 1, localHeight,
              A.Buffer(), 1, A.LDim() );
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Broadcast \
  ( Matrix<T>& A, mpi::Comm comm, int rank ); \
  EL_EXTERN template void Broadcast \
  ( AbstractDistMatrix<T>& A, mpi::Comm comm, int rank );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_BROADCAST_HPP
