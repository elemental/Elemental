/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/matrices.hpp>

namespace El {

template<typename T>
void Zeros( Matrix<T>& A, Int m, Int n )
{
    EL_DEBUG_CSE
    A.Resize( m, n );
    Zero( A );
}

template<typename T>
void Zeros( AbstractDistMatrix<T>& A, Int m, Int n )
{
    EL_DEBUG_CSE
    A.Resize( m, n );
    Zero( A );
}

template<typename T>
void Zeros( SparseMatrix<T>& A, Int m, Int n )
{
    EL_DEBUG_CSE
    A.Resize( m, n );
    Zero( A );
}

template<typename T>
void Zeros( DistSparseMatrix<T>& A, Int m, Int n )
{
    EL_DEBUG_CSE
    A.Resize( m, n );
    Zero( A );
}

template<typename T>
void Zeros( DistMultiVec<T>& A, Int m, Int n )
{
    EL_DEBUG_CSE
    A.Resize( m, n );
    Zero( A );
}

#define PROTO(T) \
  template void Zeros( Matrix<T>& A, Int m, Int n ); \
  template void Zeros( AbstractDistMatrix<T>& A, Int m, Int n ); \
  template void Zeros( SparseMatrix<T>& A, Int m, Int n ); \
  template void Zeros( DistSparseMatrix<T>& A, Int m, Int n ); \
  template void Zeros( DistMultiVec<T>& A, Int m, Int n );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
