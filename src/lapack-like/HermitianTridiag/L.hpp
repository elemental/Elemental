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

template<typename F>
void L( Matrix<F>& A, Matrix<F>& t )
{
#ifndef RELEASE
    CallStackEntry entry("hermitian_tridiag::L");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
#endif
    typedef BASE(F) R;
    const Int tHeight = Max(A.Height()-1,0);
    t.ResizeTo( tHeight, 1 );

    // Matrix views 
    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  alpha21T,
        ABL, ABR,  a10, alpha11, a12,  a21B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<F> w21;

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height()+1 < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        //--------------------------------------------------------------------//
        const F tau = Reflector( alpha21T, a21B );
        const R epsilon1 = alpha21T.GetRealPart(0,0);
        t.Set(A00.Height(),0,tau);
        alpha21T.Set(0,0,F(1));

        Zeros( w21, a21.Height(), 1 );
        Hemv( LOWER, tau, A22, a21, F(0), w21 );
        const F alpha = -tau*Dot( w21, a21 )/F(2);
        Axpy( alpha, a21, w21 );
        Her2( LOWER, F(-1), a21, w21, A22 );
        alpha21T.Set(0,0,epsilon1);
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
}

template<typename F> 
void L( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("hermitian_tridiag::L");
    if( A.Grid() != t.Grid() )
        LogicError("{A,t} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        LogicError("A must be square");
    if( t.Viewing() )
        LogicError("t must not be a view");
#endif
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> tDiag(g);
    tDiag.AlignWithDiagonal( A, -1 );
    tDiag.ResizeTo( A.Height()-1, 1 );

    // Matrix views 
    DistMatrix<F> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F,MD,STAR> tT(g),  t0(g), 
                          tB(g),  t1(g),
                                  t2(g);

    // Temporary distributions
    DistMatrix<F> WPan(g);
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> APan_MC_STAR(g),  A11_MC_STAR(g),
                                              A21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> APan_MR_STAR(g),  A11_MR_STAR(g),
                                              A21_MR_STAR(g);
    DistMatrix<F,MC,  STAR> WPan_MC_STAR(g),  W11_MC_STAR(g),
                                              W21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> WPan_MR_STAR(g),  W11_MR_STAR(g),
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
              F(-1), A21_MC_STAR, W21_MR_STAR,
                     W21_MC_STAR, A21_MR_STAR,
              F(1),  A22 );
            //----------------------------------------------------------------//
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
}

} // namespace hermitian_tridiag
} // namespace elem

#endif // ifndef LAPACK_HERMITIANTRIDIAG_L_HPP
