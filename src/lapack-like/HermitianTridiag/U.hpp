/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef LAPACK_HERMITIANTRIDIAG_U_HPP
#define LAPACK_HERMITIANTRIDIAG_U_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level1/Dot.hpp"
#include "elemental/blas-like/level2/Hemv.hpp"
#include "elemental/blas-like/level2/Her2.hpp"
#include "elemental/blas-like/level2/Symv.hpp"
#include "elemental/blas-like/level2/Syr2.hpp"
#include "elemental/lapack-like/Reflector.hpp"

#include "./PanelU.hpp"

namespace elem {
namespace hermitian_tridiag {

template<typename F>
void U( Matrix<F>& A, Matrix<F>& t )
{
#ifndef RELEASE
    CallStackEntry entry("hermitian_tridiag::U");
#endif
    const int tHeight = std::max(A.Height()-1,0);
#ifndef RELEASE
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() && (t.Height() != tHeight || t.Width() != 1) )
        throw std::logic_error("t is of the wrong size");
#endif
    typedef BASE(F) R;
    if( !t.Viewing() )
        t.ResizeTo( tHeight, 1 );

    // Matrix views 
    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  a01T,
        ABL, ABR,  a10, alpha11, a12,  alpha01B,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<F> w01;

    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ABR.Height()+1 < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22, 1 );

        PartitionUp
        ( a01, a01T,
               alpha01B, 1 );

        //--------------------------------------------------------------------//
        const F tau = Reflector( alpha01B, a01T );
        const R epsilon1 = alpha01B.GetRealPart(0,0);
        t.Set(t.Height()-1-A22.Height(),0,tau);
        alpha01B.Set(0,0,F(1));

        Zeros( w01, a01.Height(), 1 );
        Hemv( UPPER, tau, A00, a01, F(0), w01 );
        const F alpha = -tau*Dot( w01, a01 )/F(2);
        Axpy( alpha, a01, w01 );
        Her2( UPPER, F(-1), a01, w01, A00 );
        alpha01B.Set(0,0,epsilon1);
        //--------------------------------------------------------------------//

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );
    }
}

template<typename F>
void U( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
#ifndef RELEASE
    CallStackEntry entry("hermitian_tridiag::U");
    if( A.Grid() != t.Grid() )
        throw std::logic_error("{A,t} must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( t.Viewing() )
        throw std::logic_error("t must not be a view");
#endif
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> tDiag(g);
    tDiag.AlignWithDiagonal( A, 1 );
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
    DistMatrix<F,MC,  STAR> APan_MC_STAR(g),  A01_MC_STAR(g),
                                              A11_MC_STAR(g);
    DistMatrix<F,MR,  STAR> APan_MR_STAR(g),  A01_MR_STAR(g),
                                              A11_MR_STAR(g);
    DistMatrix<F,MC,  STAR> WPan_MC_STAR(g),  W01_MC_STAR(g),
                                              W11_MC_STAR(g);
    DistMatrix<F,MR,  STAR> WPan_MR_STAR(g),  W01_MR_STAR(g),
                                              W11_MR_STAR(g);

    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( tDiag, tT,
             tB, 0 );
    while( ABR.Height() < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        RepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );
        
        if( A00.Height() > 0 )
        {
            WPan.AlignWith( A01 );
            APan_MC_STAR.AlignWith( A00 );
            WPan_MC_STAR.AlignWith( A00 );
            APan_MR_STAR.AlignWith( A00 );
            WPan_MR_STAR.AlignWith( A00 );
           //-----------------------------------------------------------------//
            WPan.ResizeTo( ATL.Height(), A11.Width() );
            APan_MC_STAR.ResizeTo( ATL.Height(), A11.Width() );
            WPan_MC_STAR.ResizeTo( ATL.Height(), A11.Width() );
            APan_MR_STAR.ResizeTo( ATL.Height(), A11.Width() );
            WPan_MR_STAR.ResizeTo( ATL.Height(), A11.Width() );

            hermitian_tridiag::PanelU
            ( ATL, WPan, t1,
              APan_MC_STAR, APan_MR_STAR, WPan_MC_STAR, WPan_MR_STAR );

            PartitionUp
            ( APan_MC_STAR, A01_MC_STAR,
                            A11_MC_STAR, A11.Height() );
            PartitionUp
            ( APan_MR_STAR, A01_MR_STAR,
                            A11_MR_STAR, A11.Height() );
            PartitionUp
            ( WPan_MC_STAR, W01_MC_STAR,
                            W11_MC_STAR, A11.Height() );
            PartitionUp
            ( WPan_MR_STAR, W01_MR_STAR,
                            W11_MR_STAR, A11.Height() );

            LocalTrr2k
            ( UPPER, ADJOINT, ADJOINT,
              F(-1), A01_MC_STAR, W01_MR_STAR,
                     W01_MC_STAR, A01_MR_STAR,
              F(1),  A00 );
            //----------------------------------------------------------------//
        }
        else
        {
            A11_STAR_STAR = A11;
            t1_STAR_STAR.ResizeTo( t1.Height(), 1 );

            HermitianTridiag
            ( UPPER, A11_STAR_STAR.Matrix(), t1_STAR_STAR.Matrix() );

            A11 = A11_STAR_STAR;
            t1 = t1_STAR_STAR;
        }

        SlidePartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );
    }

    // Redistribute from matrix-diagonal form to fully replicated
    t = tDiag;
}

} // namespace hermitian_tridiag
} // namespace elem

#endif // ifndef LAPACK_HERMITIANTRIDIAG_U_HPP
