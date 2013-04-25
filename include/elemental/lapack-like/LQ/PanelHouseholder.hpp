/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_LQ_PANELHOUSEHOLDER_HPP
#define LAPACK_LQ_PANELHOUSEHOLDER_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace lq {

template<typename Real> 
inline void
PanelHouseholder( Matrix<Real>& A )
{
#ifndef RELEASE
    PushCallStack("lq::PanelHouseholder");
#endif
    Matrix<Real>
        ATL, ATR,  A00, a01,     A02,  aTopRow, ABottomPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;

    Matrix<Real> z;

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        View1x2( aTopRow, alpha11, a12 );
        View1x2( ABottomPan, a21, A22 );

        //--------------------------------------------------------------------//
        const Real tau = Reflector( alpha11, a12 );
        const Real alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);

        Zeros( z, ABottomPan.Height(), 1 );
        Gemv( NORMAL, Real(1), ABottomPan, aTopRow, Real(0), z );
        Ger( -tau, z, aTopRow, ABottomPan );
        alpha11.Set(0,0,alpha);
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real> 
inline void
PanelHouseholder( DistMatrix<Real>& A )
{
#ifndef RELEASE
    PushCallStack("lq::PanelHouseholder");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<Real>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aTopRow(g), ABottomPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);

    // Temporary distributions
    DistMatrix<Real,STAR,MR> aTopRow_STAR_MR(g);
    DistMatrix<Real,MC,STAR> z_MC_STAR(g);

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        View1x2( aTopRow, alpha11, a12 );
        View1x2( ABottomPan, a21, A22 );

        aTopRow_STAR_MR.AlignWith( ABottomPan );
        z_MC_STAR.AlignWith( ABottomPan );
        //--------------------------------------------------------------------//
        const Real tau = Reflector( alpha11, a12 );

        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlignment() &&
                                       g.Col() == alpha11.RowAlignment() );
        Real alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }

        aTopRow_STAR_MR = aTopRow;

        Zeros( z_MC_STAR, ABottomPan.Height(), 1 );
        Gemv
        ( NORMAL,
          Real(1), ABottomPan.LockedMatrix(), aTopRow_STAR_MR.LockedMatrix(),
          Real(0), z_MC_STAR.Matrix() );
        z_MC_STAR.SumOverRow();

        Ger
        ( -tau,
          z_MC_STAR.LockedMatrix(), aTopRow_STAR_MR.LockedMatrix(),
          ABottomPan.Matrix() );

        if( myDiagonalEntry )
            alpha11.SetLocal(0,0,alpha);
        //--------------------------------------------------------------------//
        aTopRow_STAR_MR.FreeAlignments();
        z_MC_STAR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real>
inline void
PanelHouseholder
( Matrix<Complex<Real> >& A,
  Matrix<Complex<Real> >& t )
{
#ifndef RELEASE
    PushCallStack("lq::PanelHouseholder");
    if( t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1 )
        throw std::logic_error
        ("t must be a vector of height equal to the minimum dimension of A");
#endif
    typedef Complex<Real> C;

    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  aTopRow, ABottomPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;
    Matrix<C>
        tT,  t0,
        tB,  tau1,
             t2;

    Matrix<C> z, aTopRowConj;

    PartitionDownLeftDiagonal
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
        const C tau = Reflector( alpha11, a12 );
        tau1.Set( 0, 0, tau );
        const C alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);

        Conjugate( aTopRow, aTopRowConj );
        Zeros( z, ABottomPan.Height(), 1 );
        Gemv( NORMAL, C(1), ABottomPan, aTopRowConj, C(0), z );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename Real>
inline void
PanelHouseholder
( DistMatrix<Complex<Real> >& A,
  DistMatrix<Complex<Real>,MD,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("lq::PanelHouseholder");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
    if( t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1 )
        throw std::logic_error
        ("t must be a vector of height equal to the minimum dimension of A");
    if( !t.AlignedWithDiagonal( A, 0 ) )
        throw std::logic_error("t must be aligned with A's main diagonal");
#endif
    typedef Complex<Real> C;
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<C>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aTopRow(g), ABottomPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    DistMatrix<C> aTopRowConj(g);
    DistMatrix<C,STAR,MR  > aTopRowConj_STAR_MR(g);
    DistMatrix<C,MC,  STAR> z_MC_STAR(g);

    PartitionDownLeftDiagonal
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
        const C tau = Reflector( alpha11, a12 );
        tau1.Set( 0, 0, tau );

        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlignment() &&
                                       g.Col() == alpha11.RowAlignment() );
        C alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }

        Conjugate( aTopRow, aTopRowConj );
        aTopRowConj_STAR_MR = aTopRowConj;

        Zeros( z_MC_STAR, ABottomPan.Height(), 1 );
        Gemv
        ( NORMAL,
          C(1), ABottomPan.LockedMatrix(),
                aTopRowConj_STAR_MR.LockedMatrix(),
          C(0), z_MC_STAR.Matrix() );
        z_MC_STAR.SumOverRow();

        Ger
        ( -Conj(tau),
          z_MC_STAR.LockedMatrix(),
          aTopRowConj_STAR_MR.LockedMatrix(),
          ABottomPan.Matrix() );

        if( myDiagonalEntry )
            alpha11.SetLocal(0,0,alpha);
        //--------------------------------------------------------------------//
        aTopRowConj_STAR_MR.FreeAlignments();
        z_MC_STAR.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace lq
} // namespace elem

#endif // ifndef LAPACK_LQ_PANELHOUSEHOLDER_HPP
