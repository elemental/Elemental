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
void MinIJ( Matrix<T>& M, Int n )
{
    EL_DEBUG_CSE
    M.Resize( n, n );
    auto minIJFill = []( Int i, Int j ) { return T(Min(i+1,j+1)); };
    IndexDependentFill( M, function<T(Int,Int)>(minIJFill) );
}

template<typename T>
void MinIJ( AbstractDistMatrix<T>& M, Int n )
{
    EL_DEBUG_CSE
    M.Resize( n, n );
    auto minIJFill = []( Int i, Int j ) { return T(Min(i+1,j+1)); };
    IndexDependentFill( M, function<T(Int,Int)>(minIJFill) );
}

#define PROTO(T) \
  template void MinIJ( Matrix<T>& M, Int n ); \
  template void MinIJ( AbstractDistMatrix<T>& M, Int n );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
