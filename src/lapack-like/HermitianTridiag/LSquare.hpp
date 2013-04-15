/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef LAPACK_HERMITIANTRIDIAG_LSQUARE_HPP
#define LAPACK_HERMITIANTRIDIAG_LSQUARE_HPP

#include "./PanelLSquare.hpp"

namespace elem {
namespace hermitian_tridiag {

template<typename R> 
void LSquare( DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("hermitian_tridiag::LSquare");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Grid().Height() != A.Grid().Width() )
        throw std::logic_error("The process grid must be square");
#endif
    const Grid& g = A.Grid();

    // Matrix views 
    DistMatrix<R> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<R> WPan(g);
    DistMatrix<R,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<R,MC,  STAR> APan_MC_STAR(g),  A11_MC_STAR(g),
                                              A21_MC_STAR(g);
    DistMatrix<R,MR,  STAR> APan_MR_STAR(g),  A11_MR_STAR(g),
                                              A21_MR_STAR(g);
    DistMatrix<R,MC,  STAR> WPan_MC_STAR(g),  W11_MC_STAR(g),
                                              W21_MC_STAR(g);
    DistMatrix<R,MR,  STAR> WPan_MR_STAR(g),  W11_MR_STAR(g),
                                              W21_MR_STAR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        if( A22.Height() > 0 )
        {
            WPan.AlignWith( A11 );
            APan_MC_STAR.AlignWith( A11 );
            WPan_MC_STAR.AlignWith( A11 );
            APan_MR_STAR.AlignWith( A11 );
            WPan_MR_STAR.AlignWith( A11 );
            //----------------------------------------------------------------//
            WPan.ResizeTo( ABR.Height(), A11.Width() );
            APan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            WPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            APan_MR_STAR.ResizeTo( ABR.Height(), A11.Width() );
            WPan_MR_STAR.ResizeTo( ABR.Height(), A11.Width() );

            hermitian_tridiag::PanelLSquare
            ( ABR, WPan, 
              APan_MC_STAR, APan_MR_STAR, WPan_MC_STAR, WPan_MR_STAR );

            PartitionDown
            ( APan_MC_STAR, A11_MC_STAR,
                            A21_MC_STAR, A11.Height() );
            PartitionDown
            ( APan_MR_STAR, A11_MR_STAR,
                            A21_MR_STAR, A11.Height() );
            PartitionDown
            ( WPan_MC_STAR, W11_MC_STAR,
                            W21_MC_STAR, A11.Height() );
            PartitionDown
            ( WPan_MR_STAR, W11_MR_STAR,
                            W21_MR_STAR, A11.Height() );

            LocalTrr2k
            ( LOWER, TRANSPOSE, TRANSPOSE,
              R(-1), A21_MC_STAR, W21_MR_STAR,
                     W21_MC_STAR, A21_MR_STAR,
              R(1), A22 );
            //----------------------------------------------------------------//
            WPan_MR_STAR.FreeAlignments();
            APan_MR_STAR.FreeAlignments();
            WPan_MC_STAR.FreeAlignments();
            APan_MC_STAR.FreeAlignments();
            WPan.FreeAlignments();
        }
        else
        {
            A11_STAR_STAR = A11;
            HermitianTridiag( LOWER, A11_STAR_STAR.Matrix() );
            A11 = A11_STAR_STAR;
        }

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
void LSquare
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("hermitian_tridiag::LSquare");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
#endif
    const Grid& g = A.Grid();
#ifndef RELEASE
    if( g.Height() != g.Width() )
        throw std::logic_error("The process grid must be square");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() )
        throw std::logic_error("t must not be a view");
#endif
    typedef Complex<R> C;

    DistMatrix<C,MD,STAR> tDiag(g);
    tDiag.AlignWithDiagonal( A, -1 );
    tDiag.ResizeTo( A.Height()-1, 1 );

    // Matrix views 
    DistMatrix<C> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<C,MD,STAR> tT(g),  t0(g), 
                          tB(g),  t1(g),
                                  t2(g);

    // Temporary distributions
    DistMatrix<C> WPan(g);
    DistMatrix<C,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<C,MC,  STAR> APan_MC_STAR(g),  A11_MC_STAR(g),
                                              A21_MC_STAR(g);
    DistMatrix<C,MR,  STAR> APan_MR_STAR(g),  A11_MR_STAR(g),
                                              A21_MR_STAR(g);
    DistMatrix<C,MC,  STAR> WPan_MC_STAR(g),  W11_MC_STAR(g),
                                              W21_MC_STAR(g);
    DistMatrix<C,MR,  STAR> WPan_MR_STAR(g),  W11_MR_STAR(g),
                                              W21_MR_STAR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( tDiag, tT,
             tB, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );
            
        if( A22.Height() > 0 )
        {
            WPan.AlignWith( A11 );
            APan_MC_STAR.AlignWith( A11 );
            WPan_MC_STAR.AlignWith( A11 );
            APan_MR_STAR.AlignWith( A11 );
            WPan_MR_STAR.AlignWith( A11 );
            //----------------------------------------------------------------//
            WPan.ResizeTo( ABR.Height(), A11.Width() );
            APan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            WPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            APan_MR_STAR.ResizeTo( ABR.Height(), A11.Width() );
            WPan_MR_STAR.ResizeTo( ABR.Height(), A11.Width() );

            hermitian_tridiag::PanelLSquare
            ( ABR, WPan, t1,
              APan_MC_STAR, APan_MR_STAR, WPan_MC_STAR, WPan_MR_STAR );

            PartitionDown
            ( APan_MC_STAR, A11_MC_STAR,
                            A21_MC_STAR, A11.Height() );
            PartitionDown
            ( APan_MR_STAR, A11_MR_STAR,
                            A21_MR_STAR, A11.Height() );
            PartitionDown
            ( WPan_MC_STAR, W11_MC_STAR,
                            W21_MC_STAR, A11.Height() );
            PartitionDown
            ( WPan_MR_STAR, W11_MR_STAR,
                            W21_MR_STAR, A11.Height() );

            LocalTrr2k
            ( LOWER, ADJOINT, ADJOINT,
              C(-1), A21_MC_STAR, W21_MR_STAR,
                     W21_MC_STAR, A21_MR_STAR,
              C(1), A22 );
            //----------------------------------------------------------------//
            WPan_MR_STAR.FreeAlignments();
            APan_MR_STAR.FreeAlignments();
            WPan_MC_STAR.FreeAlignments();
            APan_MC_STAR.FreeAlignments();
            WPan.FreeAlignments();
        }
        else
        {
            A11_STAR_STAR = A11;
            t1_STAR_STAR.ResizeTo( t1.Height(), 1 );

            HermitianTridiag
            ( LOWER, A11_STAR_STAR.Matrix(), t1_STAR_STAR.Matrix() );

            A11 = A11_STAR_STAR;
            t1 = t1_STAR_STAR;
        }

        SlidePartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }

    // Redistribute from matrix-diagonal form to fully replicated
    t = tDiag;
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace hermitian_tridiag
} // namespace elem

#endif // ifndef LAPACK_HERMITIANTRIDIAG_LSQUARE_HPP
