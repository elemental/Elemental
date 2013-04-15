/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef LAPACK_HERMITIANTRIDIAG_L_HPP
#define LAPACK_HERMITIANTRIDIAG_L_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Dot.hpp"
#include "elemental/blas-like/level2/Hemv.hpp"
#include "elemental/blas-like/level2/Her2.hpp"
#include "elemental/blas-like/level2/Symv.hpp"
#include "elemental/blas-like/level2/Syr2.hpp"
#include "elemental/lapack-like/Reflector.hpp"

#include "./PanelL.hpp"

namespace elem {
namespace hermitian_tridiag {

template<typename R>
void L( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("hermitian_tridiag::L");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height()+1 < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        const R tau = Reflector( alpha21T, a21B );
        const R epsilon1 = alpha21T.Get(0,0);
        alpha21T.Set(0,0,R(1));

        Symv( LOWER, tau, A22, a21, R(0), w21 );
        const R alpha = -tau*Dot( w21, a21 )/R(2);
        Axpy( alpha, a21, w21 );
        Syr2( LOWER, R(-1), a21, w21, A22 );
        alpha21T.Set(0,0,epsilon1);
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void L( DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("hermitian_tridiag::L");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square.");
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

            hermitian_tridiag::PanelL
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
              R(1),  A22 );
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
void L( Matrix<Complex<R> >& A, Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("hermitian_tridiag::L");
#endif
    const int tHeight = std::max(A.Height()-1,0);
#ifndef RELEASE
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() && (t.Height() != tHeight || t.Width() != 1) )
        throw std::logic_error("t is of the wrong size");
#endif
    typedef Complex<R> C;

    if( !t.Viewing() )
        t.ResizeTo( tHeight, 1 );

    // Matrix views 
    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<C> w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height()+1 < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//
        const C tau = Reflector( alpha21T, a21B );
        const R epsilon1 = alpha21T.GetRealPart(0,0);
        t.Set(A00.Height(),0,tau);
        alpha21T.Set(0,0,C(1));

        Hemv( LOWER, tau, A22, a21, C(0), w21 );
        const C alpha = -tau*Dot( w21, a21 )/C(2);
        Axpy( alpha, a21, w21 );
        Her2( LOWER, C(-1), a21, w21, A22 );
        alpha21T.Set(0,0,epsilon1);
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
void L
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("hermitian_tridiag::L");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() )
        throw std::logic_error("t must not be a view");
#endif
    typedef Complex<R> C;

    const Grid& g = A.Grid();
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

            hermitian_tridiag::PanelL
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
              C(1),  A22 );
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

#endif // ifndef LAPACK_HERMITIANTRIDIAG_L_HPP
