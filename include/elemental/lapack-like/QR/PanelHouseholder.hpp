/*
   Copyright (c) 2009-2013, Jack Poulson
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
#ifndef RELEASE
    CallStackEntry entry("qr::PanelHouseholder");
#endif
    t.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  aLeftCol, ARightPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;
    Matrix<F> z;

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        View2x1( aLeftCol, alpha11,
                           a21 );

        View2x1( ARightPan, a12,
                            A22 );

        //--------------------------------------------------------------------//
        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a21 );
        t.Set( A00.Width(), 0, tau );

        // Apply the Householder reflector
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Zeros( z, ARightPan.Width(), 1 );
        Gemv( ADJOINT, F(1), ARightPan, aLeftCol, F(0), z );
        Ger( -Conj(tau), aLeftCol, z, ARightPan );
        alpha11.Set(0,0,alpha);
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
}

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("qr::PanelHouseholder");
#endif
    Matrix<F> t;
    PanelHouseholder( A, t );
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("qr::PanelHouseholder");
    if( A.Grid() != t.Grid() )
        LogicError("{A,t} must be distributed over the same grid");
    if( !t.AlignedWithDiagonal( A, 0 ) )
        LogicError("t must be aligned with A's main diagonal");
#endif
    const Grid& g = A.Grid();
    t.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aLeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);

    // Temporary distributions
    DistMatrix<F,MC,STAR> aLeftCol_MC_STAR(g);
    DistMatrix<F,MR,STAR> z_MR_STAR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        View2x1( aLeftCol, alpha11,
                           a21 );

        View2x1( ARightPan, a12,
                            A22 );

        aLeftCol_MC_STAR.AlignWith( ARightPan );
        z_MR_STAR.AlignWith( ARightPan );
        //--------------------------------------------------------------------//
        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a21 );
        t.Set( A00.Width(), 0, tau );

        // Apply the Householder reflector
        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlignment() && 
                                       g.Col() == alpha11.RowAlignment() );
        F alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }
        aLeftCol_MC_STAR = aLeftCol;
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
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("qr::PanelHouseholder");
#endif
    DistMatrix<F,MD,STAR> t(A.Grid());
    PanelHouseholder( A, t );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_LAPACK_QR_PANEL_HPP
