/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T> 
void KMS( Matrix<T>& K, Int n, T rho )
{
    DEBUG_ONLY(CallStackEntry cse("KMS"))
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
    DEBUG_ONLY(CallStackEntry cse("KMS"))
    K.Resize( n, n );
    auto kmsFill = 
      [=]( Int i, Int j ) -> T
      { if( i < j ) { return Pow(rho,T(j-i));       } 
        else        { return Conj(Pow(rho,T(i-j))); } };
    IndexDependentFill( K, function<T(Int,Int)>(kmsFill) );
}

template<typename T>
void KMS( AbstractBlockDistMatrix<T>& K, Int n, T rho )
{
    DEBUG_ONLY(CallStackEntry cse("KMS"))
    K.Resize( n, n );
    auto kmsFill = 
      [=]( Int i, Int j ) -> T
      { if( i < j ) { return Pow(rho,T(j-i));       } 
        else        { return Conj(Pow(rho,T(i-j))); } };
    IndexDependentFill( K, function<T(Int,Int)>(kmsFill) );
}

#define PROTO(T) \
  template void KMS( Matrix<T>& K, Int n, T rho ); \
  template void KMS( AbstractDistMatrix<T>& K, Int n, T rho ); \
  template void KMS( AbstractBlockDistMatrix<T>& K, Int n, T rho );

#include "El/macros/Instantiate.h"

} // namespace El
