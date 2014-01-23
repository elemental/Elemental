/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_QR_PANEL_HPP
#define ELEM_LAPACK_QR_PANEL_HPP

#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace qr {

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t )
{
    DEBUG_ONLY(CallStackEntry cse("qr::PanelHouseholder"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );

    Matrix<F> z;

    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11   = ViewRange( A, k,   k,   k+1, k+1 );
        auto a21       = ViewRange( A, k+1, k,   m,   k+1 );
        auto aLeftCol  = ViewRange( A, k,   k,   m,   k+1 );
        auto ARightPan = ViewRange( A, k,   k+1, m,   n   );

        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a21 );
        t.Set( k, 0, tau );

        // Apply the Householder reflector
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Zeros( z, ARightPan.Width(), 1 );
        Gemv( ADJOINT, F(1), ARightPan, aLeftCol, F(0), z );
        Ger( -Conj(tau), aLeftCol, z, ARightPan );
        alpha11.Set(0,0,alpha);
    }
}

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("qr::PanelHouseholder"))
    Matrix<F> t;
    PanelHouseholder( A, t );
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("qr::PanelHouseholder");
        if( A.Grid() != t.Grid() )
            LogicError("{A,t} must be distributed over the same grid");
        if( !t.AlignedWithDiagonal( A, 0 ) )
            LogicError("t must be aligned with A's main diagonal");
    )
    const Grid& g = A.Grid();
    DistMatrix<F,MC,STAR> aLeftCol_MC_STAR(g);
    DistMatrix<F,MR,STAR> z_MR_STAR(g);

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.Resize( minDim, 1 );

    for( Int k=0; k<minDim; ++k )
    {
        auto alpha11   = ViewRange( A, k,   k,   k+1, k+1 );
        auto a21       = ViewRange( A, k+1, k,   m,   k+1 );
        auto aLeftCol  = ViewRange( A, k,   k,   m,   k+1 );
        auto ARightPan = ViewRange( A, k,   k+1, m,   n   );

        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a21 );
        t.Set( k, 0, tau );

        // Apply the Householder reflector
        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlign() && 
                                       g.Col() == alpha11.RowAlign() );
        F alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }
        aLeftCol_MC_STAR.AlignWith( ARightPan );
        aLeftCol_MC_STAR = aLeftCol;
        z_MR_STAR.AlignWith( ARightPan );
        Zeros( z_MR_STAR, ARightPan.Width(), 1 );
        LocalGemv
        ( ADJOINT, F(1), ARightPan, aLeftCol_MC_STAR, F(0), z_MR_STAR );
        z_MR_STAR.SumOverCol(); 
        Ger
        ( -Conj(tau), 
          aLeftCol_MC_STAR.LockedMatrix(), 
          z_MR_STAR.LockedMatrix(),
          ARightPan.Matrix() );
        if( myDiagonalEntry )
            alpha11.SetLocal(0,0,alpha);
    }
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("qr::PanelHouseholder"))
    DistMatrix<F,MD,STAR> t(A.Grid());
    PanelHouseholder( A, t );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_LAPACK_QR_PANEL_HPP
