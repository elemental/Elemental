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
void GCDMatrix( Matrix<T>& G, Int m, Int n )
{
    DEBUG_CSE
    G.Resize( m, n );
    auto gcdFill = []( Int i, Int j ) { return T(GCD(i+1,j+1)); };
    IndexDependentFill( G, function<T(Int,Int)>(gcdFill) );
}

template<typename T>
void GCDMatrix( AbstractDistMatrix<T>& G, Int m, Int n )
{
    DEBUG_CSE
    G.Resize( m, n );
    auto gcdFill = []( Int i, Int j ) { return T(GCD(i+1,j+1)); };
    IndexDependentFill( G, function<T(Int,Int)>(gcdFill) );
}

#define PROTO(T) \
  template void GCDMatrix( Matrix<T>& G, Int m, Int n ); \
  template void GCDMatrix( AbstractDistMatrix<T>& G, Int m, Int n );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
