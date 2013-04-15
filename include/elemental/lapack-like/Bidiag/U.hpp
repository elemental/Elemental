/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_BIDIAG_U_HPP
#define LAPACK_BIDIAG_U_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/Bidiag/PanelU.hpp"
#include "elemental/lapack-like/Bidiag/UUnb.hpp"
#include "elemental/lapack-like/Reflector.hpp"

namespace elem {
namespace bidiag {

template<typename R>
inline void U( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("bidiag::U");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  alpha12L, a12R,
        ABL, ABR,  a10, alpha11, a12,  aB1, AB2,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> x12Trans, w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        View2x1
        ( aB1, alpha11,
               a21 );
        View2x1
        ( AB2, a12,
               A22 );

        x12Trans.ResizeTo( a12.Width(), 1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//

        // Find tauQ, u, and epsilonQ such that
        //     I - tauQ | 1 | | 1, u^T | | alpha11 | = | epsilonQ |
        //              | u |            |   a21   | = |    0     |
        const R tauQ = Reflector( alpha11, a21 );
        const R epsilonQ = alpha11.Get(0,0);

        // Set aB1 = | 1 | and form x12^T := (aB1^T AB2)^T = AB2^T aB1
        //           | u |
        alpha11.Set(0,0,R(1));
        Gemv( TRANSPOSE, R(1), AB2, aB1, R(0), x12Trans );

        // Update AB2 := AB2 - tauQ aB1 x12
        //             = AB2 - tauQ aB1 aB1^T AB2
        //             = (I - tauQ aB1 aB1^T) AB2
        Ger( -tauQ, aB1, x12Trans, AB2 );

        // Put epsilonQ back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonQ);

        if( A22.Width() != 0 )
        {
            // Expose the subvector we seek to zero, a12R
            PartitionRight( a12, alpha12L, a12R );

            // Find tauP, v, and epsilonP such that
            //     I - tauP | 1 | | 1, v^T | | alpha12L | = | epsilonP |
            //              | v |            |  a12R^T  | = |    0     |
            const R tauP = Reflector( alpha12L, a12R );
            const R epsilonP = alpha12L.Get(0,0);

            // Set a12^T = | 1 | and form w21 := A22 a12^T = A22 | 1 |
            //             | v |                                 | v |
            alpha12L.Set(0,0,R(1));
            Gemv( NORMAL, R(1), A22, a12, R(0), w21 );
            // A22 := A22 - tauP w21 a12
            //      = A22 - tauP A22 a12^T a12
            //      = A22 (I - tauP a12^T a12)
            Ger( -tauP, w21, a12, A22 );

            // Put epsilonP back instead of the temporary value, 1
            alpha12L.Set(0,0,epsilonP);
        }
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
inline void 
U( DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("bidiag::U");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
#endif
    const Grid& g = A.Grid();

    // Matrix views 
    DistMatrix<R> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<R> X(g), X11(g),
                        X21(g);
    DistMatrix<R> Y(g), Y11(g),
                        Y21(g);
    DistMatrix<R,MC,  STAR> X21_MC_STAR(g);
    DistMatrix<R,MR,  STAR> Y21_MR_STAR(g);
    DistMatrix<R,MC,  STAR> AColPan_MC_STAR(g), A11_MC_STAR(g),
                                                A21_MC_STAR(g);
    DistMatrix<R,STAR,MR  > ARowPan_STAR_MR(g), A11_STAR_MR(g), A12_STAR_MR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        if( A22.Width() > 0 )
        {
            X.AlignWith( A11 );
            Y.AlignWith( A11 );
            X21_MC_STAR.AlignWith( A21 );
            Y21_MR_STAR.AlignWith( A12 );
            AColPan_MC_STAR.AlignWith( A11 );
            ARowPan_STAR_MR.AlignWith( A11 );
            //----------------------------------------------------------------//
            X.ResizeTo( ABR.Height(), A11.Width() );
            Y.ResizeTo( ABR.Width(),  A11.Height() );
            AColPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            ARowPan_STAR_MR.ResizeTo( A11.Height(), ABR.Width() );

            bidiag::PanelU( ABR, X, Y, AColPan_MC_STAR, ARowPan_STAR_MR );

            PartitionDown
            ( AColPan_MC_STAR, A11_MC_STAR,
                               A21_MC_STAR, A11.Height() );
            PartitionRight
            ( ARowPan_STAR_MR, A11_STAR_MR, A12_STAR_MR, A11.Width() );

            PartitionDown
            ( X, X11,
                 X21, A11.Height() );
            PartitionDown
            ( Y, Y11,
                 Y21, A11.Width() );
            X21_MC_STAR = X21;
            Y21_MR_STAR = Y21;

            LocalGemm
            ( NORMAL, TRANSPOSE, R(-1), A21_MC_STAR, Y21_MR_STAR, R(1), A22 );
            LocalGemm
            ( NORMAL, NORMAL, R(-1), X21_MC_STAR, A12_STAR_MR, R(1), A22 );
            //----------------------------------------------------------------//
            ARowPan_STAR_MR.FreeAlignments();
            AColPan_MC_STAR.FreeAlignments();
            Y21_MR_STAR.FreeAlignments();
            X21_MC_STAR.FreeAlignments();
            Y.FreeAlignments();
            X.FreeAlignments();
        }
        else
        {
            bidiag::UUnb( ABR );
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
inline void U
( Matrix<Complex<R> >& A,
  Matrix<Complex<R> >& tP,
  Matrix<Complex<R> >& tQ )
{
#ifndef RELEASE
    PushCallStack("bidiag::U");
#endif
    const int tPHeight = std::max(A.Width()-1,0);
    const int tQHeight = A.Width();
#ifndef RELEASE
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( tP.Viewing() && (tP.Height() != tPHeight || tP.Width() != 1) )
        throw std::logic_error("tP is the wrong height");
    if( tQ.Viewing() && (tQ.Height() != tQHeight || tQ.Width() != 1) )
        throw std::logic_error("tQ is the wrong height");
#endif
    typedef Complex<R> C;

    if( !tP.Viewing() )
        tP.ResizeTo( tPHeight, 1 );
    if( !tQ.Viewing() )
        tQ.ResizeTo( tQHeight, 1 );

    // Matrix views 
    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  alpha12L, a12R,
        ABL, ABR,  a10, alpha11, a12,  aB1, AB2,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<C> x12Adj, w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        View2x1
        ( aB1, alpha11,
               a21 );
        View2x1
        ( AB2, a12,
               A22 );

        x12Adj.ResizeTo( a12.Width(),  1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//

        // Find tauQ, u, and epsilonQ such that
        //     I - conj(tauQ) | 1 | | 1, u^H | | alpha11 | = | epsilonQ |
        //                    | u |            |    a21  |   |    0     |
        const C tauQ = Reflector( alpha11, a21 );
        const C epsilonQ = alpha11.Get(0,0);
        tQ.Set(A00.Height(),0,tauQ );

        // Set aB1 = | 1 | and form x12^H := (aB1^H AB2)^H = AB2^H aB1
        //           | u |
        alpha11.Set(0,0,C(1));
        Gemv( ADJOINT, C(1), AB2, aB1, C(0), x12Adj );

        // Update AB2 := AB2 - conj(tauQ) aB1 x12
        //             = AB2 - conj(tauQ) aB1 aB1^H AB2 
        //             = (I - conj(tauQ) aB1 aB1^H) AB2
        Ger( -Conj(tauQ), aB1, x12Adj, AB2 );
        // Put epsilonQ back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonQ);

        if( A22.Width() != 0 )
        {
            // Due to the deficiencies in the BLAS ?gemv routines, this section
            // is easier if we temporarily conjugate a12
            Conjugate( a12 );

            // Expose the subvector we seek to zero, a12R
            PartitionRight( a12, alpha12L, a12R );

            // Find tauP, v, and epsilonP such that
            //     I - conj(tauP) | 1 | | 1, v^H | | alpha12L | = | epsilonP |
            //                    | v |            |  a12R^T  |   |    0     |
            const C tauP = Reflector( alpha12L, a12R );
            const C epsilonP = alpha12L.Get(0,0);
            tP.Set(A00.Height(),0,tauP);

            // Set a12^T = | 1 | and form w21 := A22 a12^T = A22 | 1 |
            //             | v |                                 | v |
            alpha12L.Set(0,0,C(1));
            Gemv( NORMAL, C(1), A22, a12, C(0), w21 );

            // A22 := A22 - tauP w21 conj(a12)
            //      = A22 - tauP A22 a12^T conj(a12)
            //      = A22 (I - tauP a12^T conj(a12))
            //      = A22 conj(I - conj(tauP) a12^H a12)
            // which compensates for the fact that the reflector was generated
            // on the conjugated a12.
            Ger( -tauP, w21, a12, A22 );

            // Put epsilonP back instead of the temporary value, 1
            alpha12L.Set(0,0,epsilonP);

            // Undue the temporary conjugation
            Conjugate( a12 );
        }
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
inline void
U
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& tP,
  DistMatrix<Complex<R>,STAR,STAR>& tQ )
{
#ifndef RELEASE
    PushCallStack("bidiag::U");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
        throw std::logic_error
        ("{A,tP,tQ} must be distributed over the same grid");
    if( A.Height() < A.Width() )
        throw std::logic_error("A must be at least as tall as it is wide");
    if( tP.Viewing() || tQ.Viewing() )
        throw std::logic_error("tP and tQ must not be views");
#endif
    typedef Complex<R> C;

    const Grid& g = A.Grid();
    const int tPHeight = std::max(A.Width()-1,0);
    const int tQHeight = A.Width();
    DistMatrix<C,MD,STAR> tPDiag(g), tQDiag(g);
    tPDiag.AlignWithDiagonal( A, 1 );
    tQDiag.AlignWithDiagonal( A, 0 );
    tPDiag.ResizeTo( tPHeight, 1 );
    tQDiag.ResizeTo( tQHeight, 1 );

    // Matrix views 
    DistMatrix<C> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<C,MD,STAR> tPT(g),  tP0(g), 
                          tPB(g),  tP1(g),
                                   tP2(g);
    DistMatrix<C,MD,STAR> tQT(g),  tQ0(g), 
                          tQB(g),  tQ1(g),
                                   tQ2(g);

    // Temporary distributions
    DistMatrix<C> X(g), X11(g),
                        X21(g);
    DistMatrix<C> Y(g), Y11(g),
                        Y21(g);
    DistMatrix<C,MC,  STAR> X21_MC_STAR(g);
    DistMatrix<C,MR,  STAR> Y21_MR_STAR(g);
    DistMatrix<C,MC,  STAR> AColPan_MC_STAR(g), A11_MC_STAR(g),
                                                A21_MC_STAR(g);
    DistMatrix<C,STAR,MR  > ARowPan_STAR_MR(g), A11_STAR_MR(g), A12_STAR_MR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( tPDiag, tPT,
              tPB, 0 );
    PartitionDown
    ( tQDiag, tQT,
              tQB, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( tPT,  tP0,
         /***/ /***/
                tP1,
          tPB,  tP2 );
        
        RepartitionDown
        ( tQT,  tQ0,
         /***/ /***/
                tQ1,
          tQB,  tQ2 );
        
        if( A22.Width() > 0 )
        {
            X.AlignWith( A11 );
            Y.AlignWith( A11 );
            X21_MC_STAR.AlignWith( A21 );
            Y21_MR_STAR.AlignWith( A12 );
            AColPan_MC_STAR.AlignWith( A11 );
            ARowPan_STAR_MR.AlignWith( A11 );
            //----------------------------------------------------------------//
            X.ResizeTo( ABR.Height(), A11.Width() );
            Y.ResizeTo( ABR.Width(), A11.Height() );
            AColPan_MC_STAR.ResizeTo( ABR.Height(), A11.Width() );
            ARowPan_STAR_MR.ResizeTo( A11.Height(), ABR.Width() );

            bidiag::PanelU
            ( ABR, tP1, tQ1, X, Y, AColPan_MC_STAR, ARowPan_STAR_MR );

            PartitionDown
            ( AColPan_MC_STAR, A11_MC_STAR,
                               A21_MC_STAR, A11.Height() );
            PartitionRight
            ( ARowPan_STAR_MR, A11_STAR_MR, A12_STAR_MR, A11.Width() );

            PartitionDown
            ( X, X11,
                 X21, A11.Height() );
            PartitionDown
            ( Y, Y11,
                 Y21, A11.Width() );
            X21_MC_STAR = X21;
            Y21_MR_STAR = Y21;

            LocalGemm
            ( NORMAL, ADJOINT, C(-1), A21_MC_STAR, Y21_MR_STAR, C(1), A22 );
            LocalGemm
            ( NORMAL, NORMAL, C(-1), X21_MC_STAR, A12_STAR_MR, C(1), A22 );
            //----------------------------------------------------------------//
            ARowPan_STAR_MR.FreeAlignments();
            AColPan_MC_STAR.FreeAlignments();
            Y21_MR_STAR.FreeAlignments();
            X21_MC_STAR.FreeAlignments();
            Y.FreeAlignments();
            X.FreeAlignments();
        }
        else
        {
            bidiag::UUnb( ABR, tP1, tQ1 );
        }

        SlidePartitionDown
        ( tQT,  tQ0,
                tQ1,
         /***/ /***/
          tQB,  tQ2 );

        SlidePartitionDown
        ( tPT,  tP0,
                tP1,
         /***/ /***/
          tPB,  tP2 );
        
        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }

    // Redistribute from matrix-diagonal form to fully replicated
    tP = tPDiag;
    tQ = tQDiag;
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace bidiag
} // namespace elem

#endif // ifndef LAPACK_BIDIAG_U_HPP
