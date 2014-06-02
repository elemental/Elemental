/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_GEAR_HPP
#define EL_GEAR_HPP

#include "./Zeros.hpp"

namespace El {

template<typename T> 
inline void
Gear( Matrix<T>& G, Int n, Int s, Int t )
{
    DEBUG_ONLY(CallStackEntry cse("Gear"))
    if( s == 0 || s > n || s < -n )
        LogicError("Invalid s value");
    if( t == 0 || t > n || t < -n )
        LogicError("Invalid t value");
    
    Zeros( G, n, n );
    for( Int j=0; j<n-1; ++j ) 
    {
        G.Set( j, j+1, T(1) );
        G.Set( j+1, j, T(1) );
    } 

    if( s > 0 )
        G.Set( 0, s-1, T(1) );
    else
        G.Set( 0, (-s)-1, T(-1) );

    if( t > 0 )
        G.Set( n-1, n-t, T(1) );
    else
        G.Set( n-1, n+t, T(-1) );
}

template<typename T>
inline void
Gear( AbstractDistMatrix<T>& G, Int n, Int s, Int t )
{
    DEBUG_ONLY(CallStackEntry cse("Gear"))
    if( s == 0 || s > n || s < -n )
        LogicError("Invalid s value");
    if( t == 0 || t > n || t < -n )
        LogicError("Invalid t value");
    
    Zeros( G, n, n );
    for( Int j=0; j<n-1; ++j ) 
    {
        G.Set( j, j+1, T(1) );
        G.Set( j+1, j, T(1) );
    } 

    if( s > 0 )
        G.Set( 0, s-1, T(1) );
    else
        G.Set( 0, (-s)-1, T(-1) );

    if( t > 0 )
        G.Set( n-1, n-t, T(1) );
    else
        G.Set( n-1, n+t, T(-1) );
}

template<typename T>
inline void
Gear( AbstractBlockDistMatrix<T>& G, Int n, Int s, Int t )
{
    DEBUG_ONLY(CallStackEntry cse("Gear"))
    if( s == 0 || s > n || s < -n )
        LogicError("Invalid s value");
    if( t == 0 || t > n || t < -n )
        LogicError("Invalid t value");
    
    Zeros( G, n, n );
    for( Int j=0; j<n-1; ++j ) 
    {
        G.Set( j, j+1, T(1) );
        G.Set( j+1, j, T(1) );
    } 

    if( s > 0 )
        G.Set( 0, s-1, T(1) );
    else
        G.Set( 0, (-s)-1, T(-1) );

    if( t > 0 )
        G.Set( n-1, n-t, T(1) );
    else
        G.Set( n-1, n+t, T(-1) );
}

} // namespace El

#endif // ifndef EL_GEAR_HPP
