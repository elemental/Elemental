/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_HELMHOLTZ_HPP
#define ELEM_MATRICES_HELMHOLTZ_HPP

#include "elemental/matrices/Zeros.hpp"

namespace elem {

// 1D Helmholtz
template<typename F> 
inline void
Helmholtz( Matrix<F>& H, Int n, F shift )
{
#ifndef RELEASE
    CallStackEntry cse("Helmholtz");
#endif
    typedef BASE(F) R;
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

template<typename F> 
inline Matrix<F>
Helmholtz( Int n, F shift )
{
    Matrix<F> H;
    Helmholtz( H, n, shift );
    return H;
}

// 2D Helmholtz
template<typename F> 
inline void
Helmholtz( Matrix<F>& H, Int nx, Int ny, F shift )
{
#ifndef RELEASE
    CallStackEntry cse("Helmholtz");
#endif
    typedef BASE(F) R;
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

template<typename F> 
inline Matrix<F>
Helmholtz( Int nx, Int ny, F shift )
{
    Matrix<F> H;
    Helmholtz( H, nx, ny, shift );
    return H;
}

// 3D Helmholtz
template<typename F> 
inline void
Helmholtz( Matrix<F>& H, Int nx, Int ny, Int nz, F shift )
{
#ifndef RELEASE
    CallStackEntry cse("Helmholtz");
#endif
    typedef BASE(F) R;
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

template<typename F> 
inline Matrix<F>
Helmholtz( Int nx, Int ny, Int nz, F shift )
{
    Matrix<F> H;
    Helmholtz( H, nx, ny, nz, shift );
    return H;
}

// 1D Helmholtz
template<typename F,Distribution U,Distribution V>
inline void
Helmholtz( DistMatrix<F,U,V>& H, Int n, F shift )
{
#ifndef RELEASE
    CallStackEntry cse("Helmholtz");
#endif
    typedef BASE(F) R;
    Zeros( H, n, n );

    const R hInv = n+1; 
    const R hInvSquared = hInv*hInv;
    const F mainTerm = 2*hInvSquared - shift;

    const Int colShift = H.ColShift();
    const Int colStride = H.ColStride();
    const Int localWidth = H.LocalWidth();
    for( Int iLoc=0; iLoc<localWidth; ++iLoc )
    {
        const Int i = colShift + iLoc*colStride;

        H.Set( i, i, mainTerm );
        if( i != 0 )
            H.Set( i, i-1, -hInvSquared );
        if( i != n-1 )
            H.Set( i, i+1, -hInvSquared );
    }
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Helmholtz( const Grid& g, Int n, F shift )
{
    DistMatrix<F,U,V> H(g);
    Helmholtz( H, n, shift );
    return H;
}

// 2D Helmholtz
template<typename F,Distribution U,Distribution V>
inline void
Helmholtz( DistMatrix<F,U,V>& H, Int nx, Int ny, F shift )
{
#ifndef RELEASE
    CallStackEntry cse("Helmholtz");
#endif
    typedef BASE(F) R;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const R hxInv = nx+1; 
    const R hyInv = ny+1;
    const R hxInvSquared = hxInv*hxInv;
    const R hyInvSquared = hyInv*hyInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared) - shift;

    const Int colShift = H.ColShift();
    const Int colStride = H.ColStride();
    const Int localWidth = H.LocalWidth();
    for( Int iLoc=0; iLoc<localWidth; ++iLoc )
    {
        const Int i = colShift + iLoc*colStride;
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

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Helmholtz( const Grid& g, Int nx, Int ny, F shift )
{
    DistMatrix<F,U,V> H(g);
    Helmholtz( H, nx, ny, shift );
    return H;
}

// 3D Helmholtz
template<typename F,Distribution U,Distribution V>
inline void
Helmholtz( DistMatrix<F,U,V>& H, Int nx, Int ny, Int nz, F shift )
{
#ifndef RELEASE
    CallStackEntry cse("Helmholtz");
#endif
    typedef BASE(F) R;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const R hxInv = nx+1; 
    const R hyInv = ny+1;
    const R hzInv = nz+1;
    const R hxInvSquared = hxInv*hxInv;
    const R hyInvSquared = hyInv*hyInv;
    const R hzInvSquared = hzInv*hzInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared+hzInvSquared) - shift;

    const Int colShift = H.ColShift();
    const Int colStride = H.ColStride();
    const Int localWidth = H.LocalWidth();
    for( Int iLoc=0; iLoc<localWidth; ++iLoc )
    {
        const Int i = colShift + iLoc*colStride;
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

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Helmholtz( const Grid& g, Int nx, Int ny, Int nz, F shift )
{
    DistMatrix<F,U,V> H(g);
    Helmholtz( H, nx, ny, nz, shift );
    return H;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_HELMHOLTZ_HPP
