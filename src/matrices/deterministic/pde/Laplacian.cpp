/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/matrices.hpp>

namespace El {

// 1D Laplacian
// ============

template<typename F> 
void Laplacian( Matrix<F>& L, Int n )
{
    EL_DEBUG_CSE
    Helmholtz( L, n, F(0) );
    L *= -1;
}

template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int n )
{
    EL_DEBUG_CSE
    Helmholtz( L, n, F(0) );
    L *= -1;
}

template<typename F> 
void Laplacian( SparseMatrix<F>& L, Int n )
{
    EL_DEBUG_CSE
    Helmholtz( L, n, F(0) );
    L *= -1;
}

template<typename F> 
void Laplacian( DistSparseMatrix<F>& L, Int n )
{
    EL_DEBUG_CSE
    Helmholtz( L, n, F(0) );
    L *= -1;
}

// 2D Laplacian
// ============

template<typename F> 
void Laplacian( Matrix<F>& L, Int nx, Int ny )
{
    EL_DEBUG_CSE
    Helmholtz( L, nx, ny, F(0) );
    L *= -1;
}

template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny )
{
    EL_DEBUG_CSE
    Helmholtz( L, nx, ny, F(0) );
    L *= -1;
}

template<typename F> 
void Laplacian( SparseMatrix<F>& L, Int nx, Int ny )
{
    EL_DEBUG_CSE
    Helmholtz( L, nx, ny, F(0) );
    L *= -1;
}

template<typename F> 
void Laplacian( DistSparseMatrix<F>& L, Int nx, Int ny )
{
    EL_DEBUG_CSE
    Helmholtz( L, nx, ny, F(0) );
    L *= -1;
}

// 3D Laplacian
// ============

template<typename F> 
void Laplacian( Matrix<F>& L, Int nx, Int ny, Int nz )
{
    EL_DEBUG_CSE
    Helmholtz( L, nx, ny, nz, F(0) );
    L *= -1;
}

template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny, Int nz )
{
    EL_DEBUG_CSE
    Helmholtz( L, nx, ny, nz, F(0) );
    L *= -1;
}

template<typename F> 
void Laplacian( SparseMatrix<F>& L, Int nx, Int ny, Int nz )
{
    EL_DEBUG_CSE
    Helmholtz( L, nx, ny, nz, F(0) );
    L *= -1;
}

template<typename F> 
void Laplacian( DistSparseMatrix<F>& L, Int nx, Int ny, Int nz )
{
    EL_DEBUG_CSE
    Helmholtz( L, nx, ny, nz, F(0) );
    L *= -1;
}

#define PROTO(F) \
  template void Laplacian( Matrix<F>& L, Int nx ); \
  template void Laplacian( AbstractDistMatrix<F>& L, Int nx ); \
  template void Laplacian( SparseMatrix<F>& L, Int nx ); \
  template void Laplacian( DistSparseMatrix<F>& L, Int nx ); \
  template void Laplacian( Matrix<F>& L, Int nx, Int ny ); \
  template void Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny ); \
  template void Laplacian( SparseMatrix<F>& L, Int nx, Int ny ); \
  template void Laplacian( DistSparseMatrix<F>& L, Int nx, Int ny ); \
  template void Laplacian \
  ( Matrix<F>& L, Int nx, Int ny, Int nz ); \
  template void Laplacian \
  ( AbstractDistMatrix<F>& L, Int nx, Int ny, Int nz ); \
  template void Laplacian \
  ( SparseMatrix<F>& L, Int nx, Int ny, Int nz ); \
  template void Laplacian \
  ( DistSparseMatrix<F>& L, Int nx, Int ny, Int nz );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
