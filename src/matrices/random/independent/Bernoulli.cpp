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
void Bernoulli( Matrix<T>& A, Int m, Int n, double p )
{ 
    EL_DEBUG_CSE
    if( p < 0. || p > 1. )
        LogicError
        ("Invalid choice of parameter p for Bernoulli distribution: ",p);
    A.Resize( m, n );
    const double q = 1-p;
    auto doubleCoin = [=]() -> T
    {
        const double alpha = SampleUniform<double>(0,1);
        if( alpha <= q ) return T(0); 
        else             return T(1);
    };
    EntrywiseFill( A, function<T()>(doubleCoin) );
}

template<typename T>
void Bernoulli( AbstractDistMatrix<T>& A, Int m, Int n, double p )
{
    EL_DEBUG_CSE
    if( p < 0. || p > 1. )
        LogicError
        ("Invalid choice of parameter p for Bernoulli distribution: ",p);
    A.Resize( m, n );
    const double q = 1-p;
    auto doubleCoin = [=]() -> T
    {
        const double alpha = SampleUniform<double>(0,1);
        if( alpha <= q ) return T(0); 
        else             return T(1);
    };
    EntrywiseFill( A, function<T()>(doubleCoin) );
}

#define PROTO(T) \
  template void Bernoulli( Matrix<T>& A, Int m, Int n, double p ); \
  template void Bernoulli( AbstractDistMatrix<T>& A, Int m, Int n, double p );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
