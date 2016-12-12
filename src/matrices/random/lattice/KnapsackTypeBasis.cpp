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

// See Subsection 3.4 of Nguyen and Stehle's "LLL on the Average"

template<typename T>
void KnapsackTypeBasis( Matrix<T>& A, Int n, Base<T> radius )
{
    EL_DEBUG_CSE
    A.Resize( n+1, n );
    auto AT = A( IR(0,n), IR(0,n) );
    auto aB = A( IR(n),   IR(0,n) );
    Identity( AT, n, n );
    Uniform( aB, 1, n, T(0), radius );
    Round( aB );
}

template<typename T>
void KnapsackTypeBasis( AbstractDistMatrix<T>& APre, Int n, Base<T> radius )
{
    EL_DEBUG_CSE
    DistMatrixWriteProxy<T,T,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    A.Resize( n+1, n );
    auto AT = A( IR(0,n), IR(0,n) );
    auto aB = A( IR(n),   IR(0,n) );
    Identity( AT, n, n );
    Uniform( aB, 1, n, T(0), radius );
    Round( aB );
}

#define PROTO(T) \
  template void KnapsackTypeBasis \
  ( Matrix<T>& A, Int n, Base<T> radius ); \
  template void KnapsackTypeBasis \
  ( AbstractDistMatrix<T>& A, Int n, Base<T> radius );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
