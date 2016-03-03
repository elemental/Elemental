/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_FULL_HPP
#define EL_BLAS_FULL_HPP

namespace El {

template<typename T>
Matrix<T> Full( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("Full"))
    Matrix<T> B;
    Copy( A, B );
    return B;
}

// NOTE: A DistSparseMatrix version does not exist since it is not yet clear
//       whether Elemental can currently handle creating a grid in such a 
//       routine without a memory leak

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template Matrix<T> Full( const SparseMatrix<T>& A );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_FULL_HPP
