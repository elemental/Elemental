/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_RECV_HPP
#define EL_BLAS_RECV_HPP

namespace El {

// Recall that A must already be the correct size
template<typename T>
void Recv( Matrix<T>& A, mpi::Comm comm, int source )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    const Int size = height*width;
    if( height == A.LDim() )
    {
        mpi::Recv( A.Buffer(), size, source, comm );
    }
    else
    {
        vector<T> buf;
        FastResize( buf, size );
        mpi::Recv( buf.data(), size, source, comm );

        // Unpack
        copy::util::InterleaveMatrix
        ( height,        width,
          buf.data(), 1, height,
          A.Buffer(), 1, A.LDim() );
    }
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Recv( Matrix<T>& A, mpi::Comm comm, int source );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_RECV_HPP
