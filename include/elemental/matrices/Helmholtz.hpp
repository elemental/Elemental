/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_HELMHOLTZ_HPP
#define ELEM_HELMHOLTZ_HPP

#include "./Zeros.hpp"

namespace elem {

// 1D Helmholtz
// ============

template<typename F> 
inline void
Helmholtz( Matrix<F>& H, Int n, F shift )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    typedef Base<F> R;
    Zeros( H, n, n );

    const R hInv = n+1; 
    const R hInvSquared = hInv*hInv;
    const F mainTerm = 2*hInvSquared - shift;
    for( Int i=0; i<n; ++i )
    {
        H.Set( i, i, mainTerm );
        if( i != 0 )
            H.Set( i, i-1, -hInvSquared );
        if( i != n-1 )
            H.Set( i, i+1, -hInvSquared );
    }
}

template<typename F,Dist U,Dist V>
inline void
Helmholtz( DistMatrix<F,U,V>& H, Int n, F shift )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    typedef Base<F> R;
    Zeros( H, n, n );

    const R hInv = n+1; 
    const R hInvSquared = hInv*hInv;
    const F mainTerm = 2*hInvSquared - shift;

    const Int localHeight = H.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);

        H.Set( i, i, mainTerm );
        if( i != 0 )
            H.Set( i, i-1, -hInvSquared );
        if( i != n-1 )
            H.Set( i, i+1, -hInvSquared );
    }
}

// 2D Helmholtz
// ============

template<typename F> 
inline void
Helmholtz( Matrix<F>& H, Int nx, Int ny, F shift )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    typedef Base<F> R;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const R hxInv = nx+1; 
    const R hyInv = ny+1;
    const R hxInvSquared = hxInv*hxInv;
    const R hyInvSquared = hyInv*hyInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared) - shift;
    for( Int i=0; i<n; ++i )
    {
        const Int x = i % nx;
        const Int y = i/nx;

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -hxInvSquared );
        if( x != nx-1 )
            H.Set( i, i+1, -hxInvSquared );
        if( y != 0 )
            H.Set( i, i-nx, -hyInvSquared );
        if( y != ny-1 )
            H.Set( i, i+nx, -hyInvSquared );
    }
}

template<typename F,Dist U,Dist V>
inline void
Helmholtz( DistMatrix<F,U,V>& H, Int nx, Int ny, F shift )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    typedef Base<F> R;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const R hxInv = nx+1; 
    const R hyInv = ny+1;
    const R hxInvSquared = hxInv*hxInv;
    const R hyInvSquared = hyInv*hyInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared) - shift;

    const Int localHeight = H.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);
        const Int x = i % nx;
        const Int y = i/nx;

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -hxInvSquared );
        if( x != nx-1 )
            H.Set( i, i+1, -hxInvSquared );
        if( y != 0 )
            H.Set( i, i-nx, -hyInvSquared );
        if( y != ny-1 )
            H.Set( i, i+nx, -hyInvSquared );
    }
}

// 3D Helmholtz
// ============

template<typename F> 
inline void
Helmholtz( Matrix<F>& H, Int nx, Int ny, Int nz, F shift )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    typedef Base<F> R;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const R hxInv = nx+1; 
    const R hyInv = ny+1;
    const R hzInv = nz+1;
    const R hxInvSquared = hxInv*hxInv;
    const R hyInvSquared = hyInv*hyInv;
    const R hzInvSquared = hzInv*hzInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared+hzInvSquared) - shift;
    for( Int i=0; i<n; ++i )
    {
        const Int x = i % nx;
        const Int y = (i/nx) % ny;
        const Int z = i/(nx*ny);

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -hxInvSquared );
        if( x != nx-1 )
            H.Set( i, i+1, -hxInvSquared );
        if( y != 0 )
            H.Set( i, i-nx, -hyInvSquared );
        if( y != ny-1 )
            H.Set( i, i+nx, -hyInvSquared );
        if( z != 0 )
            H.Set( i, i-nx*ny, -hzInvSquared );
        if( z != nz-1 )
            H.Set( i, i+nx*ny, -hzInvSquared );
    }
}

template<typename F,Dist U,Dist V>
inline void
Helmholtz( DistMatrix<F,U,V>& H, Int nx, Int ny, Int nz, F shift )
{
    DEBUG_ONLY(CallStackEntry cse("Helmholtz"))
    typedef Base<F> R;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const R hxInv = nx+1; 
    const R hyInv = ny+1;
    const R hzInv = nz+1;
    const R hxInvSquared = hxInv*hxInv;
    const R hyInvSquared = hyInv*hyInv;
    const R hzInvSquared = hzInv*hzInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared+hzInvSquared) - shift;

    const Int localHeight = H.LocalHeight();
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);
        const Int x = i % nx;
        const Int y = (i/nx) % ny;
        const Int z = i/(nx*ny);

        H.Set( i, i, mainTerm );
        if( x != 0 )
            H.Set( i, i-1, -hxInvSquared );
        if( x != nx-1 )
            H.Set( i, i+1, -hxInvSquared );
        if( y != 0 )
            H.Set( i, i-nx, -hyInvSquared );
        if( y != ny-1 )
            H.Set( i, i+nx, -hyInvSquared );
        if( z != 0 )
            H.Set( i, i-nx*ny, -hzInvSquared );
        if( z != nz-1 )
            H.Set( i, i+nx*ny, -hzInvSquared );
    }
}

} // namespace elem

#endif // ifndef ELEM_HELMHOLTZ_HPP
