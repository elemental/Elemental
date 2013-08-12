/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LQ_PANELHOUSEHOLDER_HPP
#define ELEM_LAPACK_LQ_PANELHOUSEHOLDER_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace lq {

template<typename F>
inline void
PanelHouseholder( Matrix<F>& A, Matrix<F>& t )
{
#ifndef RELEASE
    CallStackEntry entry("lq::PanelHouseholder");
    if( t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1 )
        LogicError
        ("t must be a vector of height equal to the minimum dimension of A");
#endif
    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  aTopRow, ABottomPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;
    Matrix<F>
        tT,  t0,
        tB,  tau1,
             t2;

    Matrix<F> z, aTopRowConj;

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        RepartitionDown
        ( tT,  t0,
         /**/ /****/
               tau1,
          tB,  t2, 1 );

        View1x2( aTopRow, alpha11, a12 );
        View1x2( ABottomPan, a21, A22 );
        //--------------------------------------------------------------------//
        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a12 );
        tau1.Set( 0, 0, tau );

        // Apply the Householder reflector
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Conjugate( aTopRow, aTopRowConj );
        Zeros( z, ABottomPan.Height(), 1 );
        Gemv( NORMAL, F(1), ABottomPan, aTopRowConj, F(0), z );
        Ger( -Conj(tau), z, aTopRowConj, ABottomPan );
        alpha11.Set(0,0,alpha);
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );

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
    CallStackEntry entry("lq::PanelHouseholder");
#endif
    Matrix<F> t;
    PanelHouseholder( A, t );
}

template<typename F>
inline void
PanelHouseholder( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("lq::PanelHouseholder");
    if( A.Grid() != t.Grid() )
        LogicError("{A,t} must be distributed over the same grid");
    if( t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1 )
        LogicError
        ("t must be a vector of height equal to the minimum dimension of A");
    if( !t.AlignedWithDiagonal( A, 0 ) )
        LogicError("t must be aligned with A's main diagonal");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aTopRow(g), ABottomPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    DistMatrix<F> aTopRowConj(g);
    DistMatrix<F,STAR,MR  > aTopRowConj_STAR_MR(g);
    DistMatrix<F,MC,  STAR> z_MC_STAR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        RepartitionDown
        ( tT,  t0,
         /**/ /****/
               tau1,
          tB,  t2, 1 );

        View1x2( aTopRow, alpha11, a12 );
        View1x2( ABottomPan, a21, A22 );

        aTopRowConj_STAR_MR.AlignWith( ABottomPan );
        z_MC_STAR.AlignWith( ABottomPan );
        //--------------------------------------------------------------------//
        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a12 );
        tau1.Set( 0, 0, tau );

        // Apply the Householder reflector
        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlignment() &&
                                       g.Col() == alpha11.RowAlignment() );
        F alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }
        Conjugate( aTopRow, aTopRowConj );
        aTopRowConj_STAR_MR = aTopRowConj;
        Zeros( z_MC_STAR, ABottomPan.Height(), 1 );
        LocalGemv
        ( NORMAL, F(1), ABottomPan, aTopRowConj_STAR_MR, F(0), z_MC_STAR );
        z_MC_STAR.SumOverRow();
        Ger
        ( -Conj(tau),
          z_MC_STAR.LockedMatrix(),
          aTopRowConj_STAR_MR.LockedMatrix(),
          ABottomPan.Matrix() );
        if( myDiagonalEntry )
            alpha11.SetLocal(0,0,alpha);
        //--------------------------------------------------------------------//

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );

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
    CallStackEntry entry("lq::PanelHouseholder");
#endif
    DistMatrix<F,MD,STAR> t(A.Grid());
    PanelHouseholder( A, t );
}

} // namespace lq
} // namespace elem

#endif // ifndef ELEM_LAPACK_LQ_PANELHOUSEHOLDER_HPP
