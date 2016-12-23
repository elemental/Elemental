/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_SEND_HPP
#define EL_BLAS_SEND_HPP

namespace El {

template<typename T>
void Send( const Matrix<T>& A, mpi::Comm comm, int destination )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    const Int size = height*width;
    if( height == A.LDim() )
    {
        mpi::Send( A.LockedBuffer(), size, destination, comm );
    }
    else
    {
        vector<T> buf;
        FastResize( buf, size );

        // Pack
        copy::util::InterleaveMatrix
        ( height, width,
          A.LockedBuffer(), 1, A.LDim(),
          buf.data(),       1, height );

        mpi::Send( buf.data(), size, destination, comm );
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Send( const Matrix<T>& A, mpi::Comm comm, int rank );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_SEND_HPP
