/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_LAPLACIAN_HPP
#define MATRICES_LAPLACIAN_HPP

#include "elemental/matrices/Helmholtz.hpp"

namespace elem {

// 1D Laplacian
template<typename F> 
inline void
Laplacian( Matrix<F>& L, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Laplacian");
#endif
    Helmholtz( L, n, F(0) );
}

// 2D Laplacian
template<typename F> 
inline void
Laplacian( Matrix<F>& L, int nx, int ny )
{
#ifndef RELEASE
    CallStackEntry entry("Laplacian");
#endif
    Helmholtz( L, nx, ny, F(0) );
}

// 3D Laplacian
template<typename F> 
inline void
Laplacian( Matrix<F>& L, int nx, int ny, int nz )
{
#ifndef RELEASE
    CallStackEntry entry("Laplacian");
#endif
    Helmholtz( L, nx, ny, nz, F(0) );
}

// 1D Laplacian
template<typename F,Distribution U,Distribution V>
inline void
Laplacian( DistMatrix<F,U,V>& L, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Laplacian");
#endif
    Helmholtz( L, n, F(0) );
}

// 2D Laplacian
template<typename F,Distribution U,Distribution V>
inline void
Laplacian( DistMatrix<F,U,V>& L, int nx, int ny )
{
#ifndef RELEASE
    CallStackEntry entry("Laplacian");
#endif
    Helmholtz( L, nx, ny, F(0) );
}

// 3D Laplacian
template<typename F,Distribution U,Distribution V>
inline void
Laplacian( DistMatrix<F,U,V>& L, int nx, int ny, int nz )
{
#ifndef RELEASE
    CallStackEntry entry("Laplacian");
#endif
    Helmholtz( L, nx, ny, nz, F(0) );
}

} // namespace elem

#endif // ifndef MATRICES_LAPLACIAN_HPP
