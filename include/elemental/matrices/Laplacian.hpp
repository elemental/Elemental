/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPLACIAN_HPP
#define ELEM_LAPLACIAN_HPP

#include ELEM_HELMHOLTZ_INC

namespace elem {

// 1D Laplacian
template<typename F> 
inline void
Laplacian( Matrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
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
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
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
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, nz, F(0) );
}

template<typename F> 
inline Matrix<F>
Laplacian( Int nx, Int ny, Int nz )
{ return Helmholtz( nx, ny, nz, F(0) ); }

// 1D Laplacian
template<typename F,Dist U,Dist V>
inline void
Laplacian( DistMatrix<F,U,V>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, n, F(0) );
}

template<typename F,Dist U=MC,Dist V=MR>
inline DistMatrix<F,U,V>
Laplacian( const Grid& g, Int n )
{ return Helmholtz( g, n, F(0) ); }

// 2D Laplacian
template<typename F,Dist U,Dist V>
inline void
Laplacian( DistMatrix<F,U,V>& L, Int nx, Int ny )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, F(0) );
}

template<typename F,Dist U=MC,Dist V=MR>
inline DistMatrix<F,U,V>
Laplacian( const Grid& g, Int nx, Int ny )
{ return Helmholtz( g, nx, ny, F(0) ); }

// 3D Laplacian
template<typename F,Dist U,Dist V>
inline void
Laplacian( DistMatrix<F,U,V>& L, Int nx, Int ny, Int nz )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, nz, F(0) );
}

template<typename F,Dist U=MC,Dist V=MR>
inline DistMatrix<F,U,V>
Laplacian( const Grid& g, Int nx, Int ny, Int nz )
{ return Helmholtz( g, nx, ny, nz, F(0) ); }

} // namespace elem

#endif // ifndef ELEM_LAPLACIAN_HPP
