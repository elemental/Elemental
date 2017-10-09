/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_ROUND_HPP
#define EL_BLAS_ROUND_HPP

namespace El {

// TODO(poulson): Sparse matrix versions

template<typename T>
void Round( Matrix<T>& A )
{
    EL_DEBUG_CSE
    auto round = []( const T& alpha ) { return Round(alpha); };
    EntrywiseMap( A, MakeFunction(round) );
}

template<>
inline void Round( Matrix<Int>& /*A*/ )
{ }

#ifdef EL_HAVE_MPC
template<>
inline void Round( Matrix<BigInt>& /*A*/ )
{ }
#endif

template<typename T>
void Round( AbstractDistMatrix<T>& A )
{ Round( A.Matrix() ); }

template<typename T>
void Round( DistMultiVec<T>& A )
{ Round( A.Matrix() ); }

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void Round( Matrix<T>& A ); \
  EL_EXTERN template void Round( AbstractDistMatrix<T>& A ); \
  EL_EXTERN template void Round( DistMultiVec<T>& A );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_ROUND_HPP
