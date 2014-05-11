/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HATANONELSON_HPP
#define ELEM_HATANONELSON_HPP

#include ELEM_SETDIAGONAL_INC
#include ELEM_UNIFORM_INC
#include ELEM_ZEROS_INC

namespace elem {

// Please see Section 36 of Trefethen and Embree's "Spectra and Pseudospectra"

template<typename F> 
inline void
HatanoNelson
( Matrix<F>& A, Int n, F center, Base<F> radius, F g, bool periodic=true )
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

template<typename F,Dist U,Dist V>
inline void
HatanoNelson
( DistMatrix<F,U,V>& A, Int n, F center, Base<F> radius, F g, 
  bool periodic=true )
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

template<typename F,Dist U,Dist V>
inline void
HatanoNelson
( BlockDistMatrix<F,U,V>& A, Int n, F center, Base<F> radius, F g, 
  bool periodic=true )
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

} // namespace elem

#endif // ifndef ELEM_HATANONELSON_HPP
