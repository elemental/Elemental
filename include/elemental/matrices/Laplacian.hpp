/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_LAPLACIAN_HPP
#define ELEM_MATRICES_LAPLACIAN_HPP

#include "elemental/matrices/Helmholtz.hpp"

namespace elem {

// 1D Laplacian
template<typename F> 
inline void
Laplacian( Matrix<F>& L, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Laplacian");
#endif
    Helmholtz( L, n, F(0) );
}

template<typename F> 
inline Matrix<F>
Laplacian( Int n )
{ return Helmholtz( n, F(0) ); }

// 2D Laplacian
template<typename F> 
inline void
Laplacian( Matrix<F>& L, Int nx, Int ny )
{
#ifndef RELEASE
    CallStackEntry cse("Laplacian");
#endif
    Helmholtz( L, nx, ny, F(0) );
}

template<typename F> 
inline Matrix<F>
Laplacian( Int nx, Int ny )
{ return Helmholtz( nx, ny, F(0) ); }

// 3D Laplacian
template<typename F> 
inline void
Laplacian( Matrix<F>& L, Int nx, Int ny, Int nz )
{
#ifndef RELEASE
    CallStackEntry cse("Laplacian");
#endif
    Helmholtz( L, nx, ny, nz, F(0) );
}

template<typename F> 
inline Matrix<F>
Laplacian( Int nx, Int ny, Int nz )
{ return Helmholtz( nx, ny, nz, F(0) ); }

// 1D Laplacian
template<typename F,Distribution U,Distribution V>
inline void
Laplacian( DistMatrix<F,U,V>& L, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Laplacian");
#endif
    Helmholtz( L, n, F(0) );
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Laplacian( const Grid& g, Int n )
{ return Helmholtz( g, n, F(0) ); }

// 2D Laplacian
template<typename F,Distribution U,Distribution V>
inline void
Laplacian( DistMatrix<F,U,V>& L, Int nx, Int ny )
{
#ifndef RELEASE
    CallStackEntry cse("Laplacian");
#endif
    Helmholtz( L, nx, ny, F(0) );
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Laplacian( const Grid& g, Int nx, Int ny )
{ return Helmholtz( g, nx, ny, F(0) ); }

// 3D Laplacian
template<typename F,Distribution U,Distribution V>
inline void
Laplacian( DistMatrix<F,U,V>& L, Int nx, Int ny, Int nz )
{
#ifndef RELEASE
    CallStackEntry cse("Laplacian");
#endif
    Helmholtz( L, nx, ny, nz, F(0) );
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Laplacian( const Grid& g, Int nx, Int ny, Int nz )
{ return Helmholtz( g, nx, ny, nz, F(0) ); }

} // namespace elem

#endif // ifndef ELEM_MATRICES_LAPLACIAN_HPP
