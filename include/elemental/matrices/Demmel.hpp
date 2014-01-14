/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_DEMMEL_HPP
#define ELEM_MATRICES_DEMMEL_HPP

#include "elemental/matrices/Toeplitz.hpp"

namespace elem {

template<typename F> 
inline void
Demmel( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Demmel"))
    typedef Base<F> Real;
    if( n == 0 )
    {
        A.ResizeTo( 0, 0 );
        return;
    }
    else if( n == 1 )
    {
        A.ResizeTo( 1, 1 );
        A.Set( 0, 0, -Real(1) );
        return;
    }

    const Real B = Pow(10.,4./(n-1));

    const Int numDiags = 2*n-1;
    std::vector<F> a( numDiags, 0 );
    for( Int j=0; j<n-1; ++j )
        a[j] = -Pow(B,Real(n-1-j));
    a[n-1] = -1;
    for( Int j=n; j<numDiags; ++j )
        a[j] = 0;
    Toeplitz( A, n, n, a );
}

template<typename F> 
inline Matrix<F>
Demmel( Int n )
{
    Matrix<F> A;
    Demmel( A, n );
    return A;
}

template<typename F,Distribution U,Distribution V>
inline void
Demmel( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Demmel"))
    typedef Base<F> Real;
    if( n == 0 )
    {
        A.ResizeTo( 0, 0 );
        return;
    }
    else if( n == 1 )
    {
        A.ResizeTo( 1, 1 );
        A.Set( 0, 0, -Real(1) );
        return;
    }
    
    const Real B = Pow(10.,4./(n-1));

    const Int numDiags = 2*n-1;
    std::vector<F> a( numDiags, 0 );
    for( Int j=0; j<n-1; ++j )
        a[j] = -Pow(B,Real(n-1-j));
    a[n-1] = -1;
    for( Int j=n; j<numDiags; ++j )
        a[j] = 0;
    Toeplitz( A, n, n, a );
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Demmel( const Grid& g, Int n )
{
    DistMatrix<F,U,V> A(g);
    Demmel( A, n );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_DEMMEL_HPP
