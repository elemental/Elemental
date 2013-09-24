/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_GEAR_HPP
#define ELEM_MATRICES_GEAR_HPP

#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T> 
inline void
Gear( Matrix<T>& G, Int n, Int s, Int t )
{
#ifndef RELEASE
    CallStackEntry cse("Gear");
#endif
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
inline Matrix<T>
Gear( Int n, Int s, Int t )
{
    Matrix<T> G;
    Gear( G, n, s, t );
    return G;
}

template<typename T,Distribution U,Distribution V>
inline void
Gear( DistMatrix<T,U,V>& G, Int n, Int s, Int t )
{
#ifndef RELEASE
    CallStackEntry cse("Gear");
#endif
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

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Gear( const Grid& g, Int n, Int s, Int t )
{
    DistMatrix<T,U,V> G(g);
    Gear( G, n, s, t );
    return G;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_GEAR_HPP
