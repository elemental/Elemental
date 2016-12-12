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
void Gear( Matrix<T>& G, Int n, Int s, Int t )
{
    EL_DEBUG_CSE
    if( s == 0 || s > n || s < -n )
        LogicError("Invalid s value");
    if( t == 0 || t > n || t < -n )
        LogicError("Invalid t value");
    Zeros( G, n, n );
    FillDiagonal( G, T(1), -1 );
    FillDiagonal( G, T(1),  1 );
    G.Set( 0,   Abs(s)-1, Sgn(s) );
    G.Set( n-1, n-Abs(t), Sgn(t) );
}

template<typename T>
void Gear( AbstractDistMatrix<T>& G, Int n, Int s, Int t )
{
    EL_DEBUG_CSE
    if( s == 0 || s > n || s < -n )
        LogicError("Invalid s value");
    if( t == 0 || t > n || t < -n )
        LogicError("Invalid t value");
    Zeros( G, n, n );
    FillDiagonal( G, T(1), -1 );
    FillDiagonal( G, T(1),  1 );
    G.Set( 0,   Abs(s)-1, Sgn(s) );
    G.Set( n-1, n-Abs(t), Sgn(t) );
}

#define PROTO(T) \
  template void Gear( Matrix<T>& G, Int n, Int s, Int t ); \
  template void Gear( AbstractDistMatrix<T>& G, Int n, Int s, Int t );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
