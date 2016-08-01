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
void ThreeValued( Matrix<T>& A, Int m, Int n, double p )
{
    DEBUG_CSE
    A.Resize( m, n );
    auto tripleCoin = [=]() -> T
    { 
        const double alpha = SampleUniform<double>(0,1);
        if( alpha <= p/2 ) return T(-1);
        else if( alpha <= p ) return T(1);
        else return T(0);
    };
    EntrywiseFill( A, function<T()>(tripleCoin) );
}

template<typename T>
void ThreeValued( AbstractDistMatrix<T>& A, Int m, Int n, double p )
{
    DEBUG_CSE
    A.Resize( m, n );
    if( A.RedundantRank() == 0 )
        ThreeValued( A.Matrix(), A.LocalHeight(), A.LocalWidth(), p );
    Broadcast( A, A.RedundantComm(), 0 );
}

#define PROTO(T) \
  template void ThreeValued \
  ( Matrix<T>& A, Int m, Int n, double p ); \
  template void ThreeValued \
  ( AbstractDistMatrix<T>& A, Int m, Int n, double p );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
