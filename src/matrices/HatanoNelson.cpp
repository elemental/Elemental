/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Please see Section 36 of Trefethen and Embree's "Spectra and Pseudospectra"

template<typename F> 
void HatanoNelson
( Matrix<F>& A, Int n, F center, Base<F> radius, F g, bool periodic )
{
    DEBUG_ONLY(CallStackEntry cse("HatanoNelson"))
    if( n < 3 )
        LogicError("Hatano Nelson requires at least a 3x3 matrix");
    Zeros( A, n, n );
    auto d = A.GetDiagonal();
    MakeUniform( d, center, radius );
    A.SetDiagonal( d );
    SetDiagonal( A, Exp(g),   1 );
    SetDiagonal( A, Exp(-g), -1 );
    if( periodic )
    {
        A.Set( 0,   n-1, Exp(-g) );
        A.Set( n-1, 0,   Exp( g) );
    }
}

template<typename F>
void HatanoNelson
( AbstractDistMatrix<F>& APre, Int n, 
  F center, Base<F> radius, F g, bool periodic )
{
    DEBUG_ONLY(CallStackEntry cse("HatanoNelson"))
    if( n < 3 )
        LogicError("Hatano Nelson requires at least a 3x3 matrix");

    auto APtr = WriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    Zeros( A, n, n );
    auto d = A.GetDiagonal();
    MakeUniform( d, center, radius );
    A.SetDiagonal( d );
    SetDiagonal( A, Exp(g),   1 );
    SetDiagonal( A, Exp(-g), -1 );
    if( periodic )
    {
        A.Set( 0,   n-1, Exp(-g) );
        A.Set( n-1, 0,   Exp( g) );
    }
}

/*
template<typename F>
void HatanoNelson
( AbstractBlockDistMatrix<F>& APre, Int n, 
  F center, Base<F> radius, F g, bool periodic )
{
    DEBUG_ONLY(CallStackEntry cse("HatanoNelson"))
    if( n < 3 )
        LogicError("Hatano Nelson requires at least a 3x3 matrix");

    auto APtr = WriteProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    Zeros( A, n, n );
    auto d = A.GetDiagonal();
    MakeUniform( d, center, radius );
    A.SetDiagonal( d );
    SetDiagonal( A, Exp(g),   1 );
    SetDiagonal( A, Exp(-g), -1 );
    if( periodic )
    {
        A.Set( 0,   n-1, Exp(-g) );
        A.Set( n-1, 0,   Exp( g) );
    }
}
*/

#define PROTO(F) \
  template void HatanoNelson \
  ( Matrix<F>& A, Int n, F center, Base<F> radius, F g, bool periodic ); \
  template void HatanoNelson \
  ( AbstractDistMatrix<F>& A, Int n, \
    F center, Base<F> radius, F g, bool periodic );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
