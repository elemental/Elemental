/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// 1D Helmholtz
// ============

template<typename F> 
void Helmholtz( Matrix<F>& H, Int n, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
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

template<typename F>
void Helmholtz( AbstractDistMatrix<F>& H, Int n, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    Zeros( H, n, n );

    const Real hInv = n+1; 
    const Real hInvSquared = hInv*hInv;
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

template<typename F>
void Helmholtz( SparseMatrix<F>& H, Int n, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    Zeros( H, n, n );

    const Real hInv = n+1; 
    const Real hInvSquared = hInv*hInv;
    const F mainTerm = 2*hInvSquared - shift;

    H.Reserve( 3*n );
    for( Int i=0; i<n; ++i )
    {
        H.QueueUpdate( i, i, mainTerm );
        if( i != 0 )
            H.QueueUpdate( i, i-1, -hInvSquared );
        if( i != n-1 )
            H.QueueUpdate( i, i+1, -hInvSquared );
    }
    H.ProcessQueues();
}

template<typename F>
void Helmholtz( DistSparseMatrix<F>& H, Int n, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    Zeros( H, n, n );

    const Real hInv = n+1; 
    const Real hInvSquared = hInv*hInv;
    const F mainTerm = 2*hInvSquared - shift;

    const Int localHeight = H.LocalHeight(); 
    H.Reserve( 3*localHeight );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);
        H.QueueUpdate( i, i, mainTerm );
        if( i != 0 )
            H.QueueUpdate( i, i-1, -hInvSquared );
        if( i != n-1 )
            H.QueueUpdate( i, i+1, -hInvSquared );
    }
    H.ProcessQueues();
}

// 2D Helmholtz
// ============

template<typename F> 
void Helmholtz( Matrix<F>& H, Int nx, Int ny, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const Real hxInv = nx+1; 
    const Real hyInv = ny+1;
    const Real hxInvSquared = hxInv*hxInv;
    const Real hyInvSquared = hyInv*hyInv;
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
void Helmholtz( AbstractDistMatrix<F>& H, Int nx, Int ny, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const Real hxInv = nx+1; 
    const Real hyInv = ny+1;
    const Real hxInvSquared = hxInv*hxInv;
    const Real hyInvSquared = hyInv*hyInv;
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

template<typename F>
void Helmholtz( SparseMatrix<F>& H, Int nx, Int ny, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const Real hxInv = nx+1;
    const Real hyInv = ny+1;
    const Real hxInvSquared = hxInv*hxInv;
    const Real hyInvSquared = hyInv*hyInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared) - shift;

    H.Reserve( 5*n );
    for( Int i=0; i<n; ++i )
    {
        const Int x = i % nx;
        const Int y = i/nx;

        H.QueueUpdate( i, i, mainTerm );
        if( x != 0 )
            H.QueueUpdate( i, i-1, -hxInvSquared );
        if( x != nx-1 )
            H.QueueUpdate( i, i+1, -hxInvSquared );
        if( y != 0 )
            H.QueueUpdate( i, i-nx, -hyInvSquared );
        if( y != ny-1 )
            H.QueueUpdate( i, i+nx, -hyInvSquared );
    }
    H.ProcessQueues();
}

template<typename F>
void Helmholtz( DistSparseMatrix<F>& H, Int nx, Int ny, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    const Int n = nx*ny;
    Zeros( H, n, n );

    const Real hxInv = nx+1;
    const Real hyInv = ny+1;
    const Real hxInvSquared = hxInv*hxInv;
    const Real hyInvSquared = hyInv*hyInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared) - shift;

    const Int localHeight = H.LocalHeight();
    H.Reserve( 5*localHeight );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);
        const Int x = i % nx;
        const Int y = i/nx;

        H.QueueUpdate( i, i, mainTerm );
        if( x != 0 )
            H.QueueUpdate( i, i-1, -hxInvSquared );
        if( x != nx-1 )
            H.QueueUpdate( i, i+1, -hxInvSquared );
        if( y != 0 )
            H.QueueUpdate( i, i-nx, -hyInvSquared );
        if( y != ny-1 )
            H.QueueUpdate( i, i+nx, -hyInvSquared );
    }
    H.ProcessQueues();
}

// 3D Helmholtz
// ============

template<typename F> 
void Helmholtz( Matrix<F>& H, Int nx, Int ny, Int nz, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const Real hxInv = nx+1; 
    const Real hyInv = ny+1;
    const Real hzInv = nz+1;
    const Real hxInvSquared = hxInv*hxInv;
    const Real hyInvSquared = hyInv*hyInv;
    const Real hzInvSquared = hzInv*hzInv;
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
void Helmholtz( AbstractDistMatrix<F>& H, Int nx, Int ny, Int nz, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const Real hxInv = nx+1; 
    const Real hyInv = ny+1;
    const Real hzInv = nz+1;
    const Real hxInvSquared = hxInv*hxInv;
    const Real hyInvSquared = hyInv*hyInv;
    const Real hzInvSquared = hzInv*hzInv;
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

template<typename F> 
void Helmholtz( SparseMatrix<F>& H, Int nx, Int ny, Int nz, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const Real hxInv = nx+1; 
    const Real hyInv = ny+1;
    const Real hzInv = nz+1;
    const Real hxInvSquared = hxInv*hxInv;
    const Real hyInvSquared = hyInv*hyInv;
    const Real hzInvSquared = hzInv*hzInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared+hzInvSquared) - shift;

    H.Reserve( 7*n );
    for( Int i=0; i<n; ++i )
    {
        const Int x = i % nx;
        const Int y = (i/nx) % ny;
        const Int z = i/(nx*ny);

        H.QueueUpdate( i, i, mainTerm );
        if( x != 0 )
            H.QueueUpdate( i, i-1, -hxInvSquared );
        if( x != nx-1 )
            H.QueueUpdate( i, i+1, -hxInvSquared );
        if( y != 0 )
            H.QueueUpdate( i, i-nx, -hyInvSquared );
        if( y != ny-1 )
            H.QueueUpdate( i, i+nx, -hyInvSquared );
        if( z != 0 )
            H.QueueUpdate( i, i-nx*ny, -hzInvSquared );
        if( z != nz-1 )
            H.QueueUpdate( i, i+nx*ny, -hzInvSquared );
    }
    H.ProcessQueues();
}

template<typename F> 
void Helmholtz( DistSparseMatrix<F>& H, Int nx, Int ny, Int nz, F shift )
{
    DEBUG_ONLY(CSE cse("Helmholtz"))
    typedef Base<F> Real;
    const Int n = nx*ny*nz;
    Zeros( H, n, n );

    const Real hxInv = nx+1; 
    const Real hyInv = ny+1;
    const Real hzInv = nz+1;
    const Real hxInvSquared = hxInv*hxInv;
    const Real hyInvSquared = hyInv*hyInv;
    const Real hzInvSquared = hzInv*hzInv;
    const F mainTerm = 2*(hxInvSquared+hyInvSquared+hzInvSquared) - shift;

    const Int localHeight = H.LocalHeight();
    H.Reserve( 7*localHeight );
    for( Int iLoc=0; iLoc<localHeight; ++iLoc )
    {
        const Int i = H.GlobalRow(iLoc);
        const Int x = i % nx;
        const Int y = (i/nx) % ny;
        const Int z = i/(nx*ny);

        H.QueueUpdate( i, i, mainTerm );
        if( x != 0 )
            H.QueueUpdate( i, i-1, -hxInvSquared );
        if( x != nx-1 )
            H.QueueUpdate( i, i+1, -hxInvSquared );
        if( y != 0 )
            H.QueueUpdate( i, i-nx, -hyInvSquared );
        if( y != ny-1 )
            H.QueueUpdate( i, i+nx, -hyInvSquared );
        if( z != 0 )
            H.QueueUpdate( i, i-nx*ny, -hzInvSquared );
        if( z != nz-1 )
            H.QueueUpdate( i, i+nx*ny, -hzInvSquared );
    }
    H.ProcessQueues();
}

#define PROTO(F) \
  template void Helmholtz \
  ( Matrix<F>& H, Int nx, F shift ); \
  template void Helmholtz \
  ( AbstractDistMatrix<F>& H, Int nx, F shift ); \
  template void Helmholtz \
  ( SparseMatrix<F>& H, Int nx, F shift ); \
  template void Helmholtz \
  ( DistSparseMatrix<F>& H, Int nx, F shift ); \
  template void Helmholtz \
  ( Matrix<F>& H, Int nx, Int ny, F shift ); \
  template void Helmholtz \
  ( AbstractDistMatrix<F>& H, Int nx, Int ny, F shift ); \
  template void Helmholtz \
  ( SparseMatrix<F>& H, Int nx, Int ny, F shift ); \
  template void Helmholtz \
  ( DistSparseMatrix<F>& H, Int nx, Int ny, F shift ); \
  template void Helmholtz \
  ( Matrix<F>& H, Int nx, Int ny, Int nz, F shift ); \
  template void Helmholtz \
  ( AbstractDistMatrix<F>& H, Int nx, Int ny, Int nz, F shift ); \
  template void Helmholtz \
  ( SparseMatrix<F>& H, Int nx, Int ny, Int nz, F shift ); \
  template void Helmholtz \
  ( DistSparseMatrix<F>& H, Int nx, Int ny, Int nz, F shift );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
