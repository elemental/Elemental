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
void KMS( Matrix<T>& K, Int n, T rho )
{
    DEBUG_CSE
    K.Resize( n, n );
    auto kmsFill = 
      [=]( Int i, Int j ) -> T
      { if( i < j ) { return Pow(rho,T(j-i));       } 
        else        { return Conj(Pow(rho,T(i-j))); } };
    IndexDependentFill( K, function<T(Int,Int)>(kmsFill) );
}

template<typename T>
void KMS( AbstractDistMatrix<T>& K, Int n, T rho )
{
    DEBUG_CSE
    K.Resize( n, n );
    auto kmsFill = 
      [=]( Int i, Int j ) -> T
      { if( i < j ) { return Pow(rho,T(j-i));       } 
        else        { return Conj(Pow(rho,T(i-j))); } };
    IndexDependentFill( K, function<T(Int,Int)>(kmsFill) );
}

#define PROTO(T) \
  template void KMS( Matrix<T>& K, Int n, T rho ); \
  template void KMS( AbstractDistMatrix<T>& K, Int n, T rho );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
