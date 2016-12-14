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

template<typename F>
void Hilbert( Matrix<F>& A, Int n )
{
    EL_DEBUG_CSE
    A.Resize( n, n );
    auto hilbertFill = []( Int i, Int j ) { return F(1)/F(i+j+1); };
    IndexDependentFill( A, function<F(Int,Int)>(hilbertFill) );
}

template<typename F>
void Hilbert( AbstractDistMatrix<F>& A, Int n )
{
    EL_DEBUG_CSE
    A.Resize( n, n );
    auto hilbertFill = []( Int i, Int j ) { return F(1)/F(i+j+1); };
    IndexDependentFill( A, function<F(Int,Int)>(hilbertFill) );
}

#define PROTO(F) \
  template void Hilbert( Matrix<F>& A, Int n ); \
  template void Hilbert( AbstractDistMatrix<F>& A, Int n );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
