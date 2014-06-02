/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAPLACIAN_HPP
#define EL_LAPLACIAN_HPP

#include EL_HELMHOLTZ_INC

namespace El {

// 1D Laplacian
// ============

template<typename F> 
inline void
Laplacian( Matrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, n, F(0) );
}

template<typename F>
inline void
Laplacian( AbstractDistMatrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, n, F(0) );
}

template<typename F>
inline void
Laplacian( AbstractBlockDistMatrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, n, F(0) );
}

// 2D Laplacian
// ============

template<typename F> 
inline void
Laplacian( Matrix<F>& L, Int nx, Int ny )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, F(0) );
}

template<typename F>
inline void
Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, F(0) );
}

template<typename F>
inline void
Laplacian( AbstractBlockDistMatrix<F>& L, Int nx, Int ny )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, F(0) );
}

// 3D Laplacian
// ============

template<typename F> 
inline void
Laplacian( Matrix<F>& L, Int nx, Int ny, Int nz )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, nz, F(0) );
}

template<typename F>
inline void
Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny, Int nz )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, nz, F(0) );
}

template<typename F>
inline void
Laplacian( AbstractBlockDistMatrix<F>& L, Int nx, Int ny, Int nz )
{
    DEBUG_ONLY(CallStackEntry cse("Laplacian"))
    Helmholtz( L, nx, ny, nz, F(0) );
}

} // namespace El

#endif // ifndef EL_LAPLACIAN_HPP
