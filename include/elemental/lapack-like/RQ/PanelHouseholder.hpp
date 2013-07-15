/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_RQ_PANEL_HPP
#define LAPACK_RQ_PANEL_HPP

#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace rq {

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t )
{
#ifndef RELEASE
    CallStackEntry entry("rq::PanelHouseholder");
    if( t.Viewing() && 
        (t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1) )
        throw std::logic_error
        ("t must be a vector of height equal to the minimum dimension of A");
#endif
    if( !t.Viewing() )
        t.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  ATopPan, aBottomRow,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;
    Matrix<F>
        tT,  t0,
        tB,  t1,
             t2;
    Matrix<F> z, aBottomRowConj;

    PartitionUpOffsetDiagonal
    ( A.Width()-A.Height(),
      A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    while( ABR.Height() < A.Height() && ABR.Width() < A.Width() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22, 1 );

        RepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, 1 );

        View1x2( ATopPan, A00, a01 );
        View1x2( aBottomRow, a10, alpha11 );
        //--------------------------------------------------------------------//
        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a10 );
        t1.Set( 0, 0, tau );

        // Apply the Householder reflector
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Conjugate( aBottomRow, aBottomRowConj );
        Zeros( z, ATopPan.Height(), 1 );
        Gemv( NORMAL, F(1), ATopPan, aBottomRowConj, F(0), z );
        Ger( -Conj(tau), z, aBottomRowConj, ATopPan );
        alpha11.Set(0,0,alpha);
        //--------------------------------------------------------------------//

        SlidePartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );
    }
}

template<typename F> 
inline void
PanelHouseholder( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("rq::PanelHouseholder");
#endif
    Matrix<F> t;
    PanelHouseholder( A, t );
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("rq::PanelHouseholder");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
    if( t.Viewing() && 
        (t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1) )
        throw std::logic_error
        ("t must be a vector of height equal to the minimum dimension of A");
    if( !t.AlignedWithDiagonal( A, A.Width()-A.Height() ) )
        throw std::logic_error("t must be aligned with A's main diagonal");
#endif
    const Grid& g = A.Grid();
    if( !t.Viewing() )
        t.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ATopPan(g), 
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  aBottomRow(g),
                         A20(g), a21(g),     A22(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    // Temporary distributions
    DistMatrix<F> aBottomRowConj(g);
    DistMatrix<F,STAR,MR  > aBottomRowConj_STAR_MR(g);
    DistMatrix<F,MC,  STAR> z_MC_STAR(g);

    PartitionUpOffsetDiagonal
    ( A.Width()-A.Height(),
      A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    while( ABR.Height() < A.Height() && ABR.Width() < A.Width() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22, 1 );

        RepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, 1 );

        View1x2( ATopPan, A00, a01 );
        View1x2( aBottomRow, a10, alpha11 );

        aBottomRowConj_STAR_MR.AlignWith( ATopPan );
        z_MC_STAR.AlignWith( ATopPan );
        //--------------------------------------------------------------------//
        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a10 );
        t1.Set( 0, 0, tau );

        // Apply the Householder reflector
        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlignment() && 
                                       g.Col() == alpha11.RowAlignment() );
        F alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }
        Conjugate( aBottomRow, aBottomRowConj );
        aBottomRowConj_STAR_MR = aBottomRowConj;
        Zeros( z_MC_STAR, ATopPan.Height(), 1 );
        LocalGemv
        ( NORMAL, F(1), ATopPan, aBottomRowConj_STAR_MR, F(0), z_MC_STAR );
        z_MC_STAR.SumOverRow(); 
        Ger
        ( -Conj(tau), 
          z_MC_STAR.LockedMatrix(),
          aBottomRowConj_STAR_MR.LockedMatrix(),
          ATopPan.Matrix() ); 
        if( myDiagonalEntry )
            alpha11.SetLocal(0,0,alpha);
        //--------------------------------------------------------------------//

        SlidePartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );
    }
}

template<typename F> 
inline void
PanelHouseholder( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("rq::PanelHouseholder");
#endif
    DistMatrix<F,MD,STAR> t(A.Grid());
    PanelHouseholder( A, t );
}

} // namespace rq
} // namespace elem

#endif // ifndef LAPACK_RQ_PANEL_HPP
