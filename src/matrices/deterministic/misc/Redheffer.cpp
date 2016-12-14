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
void Redheffer( Matrix<T>& R, Int n )
{
    EL_DEBUG_CSE
    R.Resize( n, n );
    auto redhefferFill = 
      []( Int i, Int j ) -> T
      { if( j == 0 || ((j+1)%(i+1))==0 ) { return T(1); }
        else                             { return T(0); } };
    IndexDependentFill( R, function<T(Int,Int)>(redhefferFill) );
}

template<typename T>
void Redheffer( AbstractDistMatrix<T>& R, Int n )
{
    EL_DEBUG_CSE
    R.Resize( n, n );
    auto redhefferFill = 
      []( Int i, Int j ) -> T
      { if( j == 0 || ((j+1)%(i+1))==0 ) { return T(1); }
        else                             { return T(0); } };
    IndexDependentFill( R, function<T(Int,Int)>(redhefferFill) );
}

#define PROTO(T) \
  template void Redheffer( Matrix<T>& R, Int n ); \
  template void Redheffer( AbstractDistMatrix<T>& R, Int n );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
