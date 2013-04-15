/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_BIDIAG_L_HPP
#define LAPACK_BIDIAG_L_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/Bidiag/LUnb.hpp"
#include "elemental/lapack-like/Bidiag/PanelL.hpp"
#include "elemental/lapack-like/Reflector.hpp"

namespace elem {
namespace bidiag {

template<typename R>
inline void L( Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("bidiag::L");
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
#endif
    // Matrix views 
    Matrix<R>
        ATL, ATR,  A00, a01,     A02,  alpha21T,  a1R,
        ABL, ABR,  a10, alpha11, a12,  a21B,      A2R,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<R> x12Trans, w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        View1x2( a1R, alpha11, a12 );
        View1x2( A2R, a21, A22 );

        x12Trans.ResizeTo( a12.Width(), 1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//

        // Find tauP, v, and epsilonP such that
        //     I - tauP | 1 | | 1, v^T | | alpha11 | = | epsilonP |
        //              | v |            |  a12^T  | = |    0     |
        const R tauP = Reflector( alpha11, a12 );
        const R epsilonP = alpha11.Get(0,0);

        // Set a1R^T = | 1 | and form w21 := A2R a1R^T = A2R | 1 |
        //             | v |                                 | v |
        alpha11.Set(0,0,R(1));
        Gemv( NORMAL, R(1), A2R, a1R, R(0), w21 );

        // A2R := A2R - tauP w21 a1R
        //      = A2R - tauP A2R a1R^T a1R
        //      = A2R (I - tauP a1R^T a1R)
        Ger( -tauP, w21, a1R, A2R );

        // Put epsilonP back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonP);

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            PartitionDown
            ( a21, alpha21T,
                   a21B );

            // Find tauQ, u, and epsilonQ such that
            //     I - tauQ | 1 | | 1, u^T | | alpha21T | = | epsilonQ |
            //              | u |            |   a21B   | = |    0     |
            const R tauQ = Reflector( alpha21T, a21B );
            const R epsilonQ = alpha21T.Get(0,0);

            // Set a21 = | 1 | and form x12^T = (a21^T A22)^T = A22^T a21
            //           | u |
            alpha21T.Set(0,0,R(1));
            Gemv( TRANSPOSE, R(1), A22, a21, R(0), x12Trans );

            // A22 := A22 - tauQ a21 x12
            //      = A22 - tauQ a21 a21^T A22
            //      = (I - tauQ a21 a21^T) A22
            Ger( -tauQ, a21, x12Trans, A22 );

            // Put epsilonQ back instead of the temporary value, 1
            alpha21T.Set(0,0,epsilonQ);
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
L( DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("bidiag::L");
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
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
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        if( A22.Height() > 0 )
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

            bidiag::PanelL( ABR, X, Y, AColPan_MC_STAR, ARowPan_STAR_MR );

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
            bidiag::LUnb( ABR );
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
inline void L
( Matrix<Complex<R> >& A,
  Matrix<Complex<R> >& tP,
  Matrix<Complex<R> >& tQ )
{
#ifndef RELEASE
    PushCallStack("bidiag::L");
#endif
    const int tPHeight = A.Height();
    const int tQHeight = std::max(A.Height()-1,0);
#ifndef RELEASE
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
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
        ATL, ATR,  A00, a01,     A02,  alpha21T,  a1R,
        ABL, ABR,  a10, alpha11, a12,  a21B,      A2R,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<C> x12Adj, w21;

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        View1x2( a1R, alpha11, a12 );
        View1x2( A2R, a21, A22 );

        x12Adj.ResizeTo( a12.Width(), 1 );
        w21.ResizeTo( a21.Height(), 1 );
        //--------------------------------------------------------------------//

        // Due to deficiencies in the BLAS ?gemv routines, this section is 
        // easier if we temporarily conjugate a1R = | alpha11, a12 |
        Conjugate( a1R );

        // Find tauP, v, and epsilonP such that
        //     I - conj(tauP) | 1 | | 1, v^H | | alpha11 | = | epsilonP | 
        //                    | v |            |   a12^T |   |    0     |
        const C tauP = Reflector( alpha11, a12 );
        const C epsilonP = alpha11.Get(0,0);
        tP.Set(A00.Height(),0,tauP);

        // Set a1R^T = | 1 | and form w21 := A2R a1R^T = A2R | 1 |
        //             | v |                                 | v |
        alpha11.Set(0,0,C(1));
        Gemv( NORMAL, C(1), A2R, a1R, C(0), w21 );

        // A2R := A2R - tauP w21 conj(a1R)
        //      = A2R - tauP A2R a1R^T conj(a1R)
        //      = A22 (I - tauP a1R^T conj(a1R))
        //      = A22 conj(I - conj(tauP) a1R^H a1R)
        // which compensates for the fact that the reflector was generated
        // on the conjugated a1R.
        Ger( -tauP, w21, a1R, A2R );

        // Put epsilonP back instead of the temporary value, 1
        alpha11.Set(0,0,epsilonP);

        // Undo the temporary conjugation
        Conjugate( a1R );

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            PartitionDown
            ( a21, alpha21T,
                   a21B );

            // Find tauQ, u, and epsilonQ such that
            //     I - conj(tauQ) | 1 | | 1, u^H | | alpha21T | = | epsilonQ |
            //                    | u |            |   a21B   | = |    0     |
            const C tauQ = Reflector( alpha21T, a21B );
            const C epsilonQ = alpha21T.Get(0,0);
            tQ.Set(A00.Height(),0,tauQ);

            // Set a21 = | 1 | and form x12^H = (a21^H A22)^H = A22^H a21
            //           | u |
            alpha21T.Set(0,0,C(1));
            Gemv( ADJOINT, C(1), A22, a21, C(0), x12Adj );

            // A22 := A22 - conj(tauQ) a21 x12 
            //      = A22 - conj(tauQ) a21 a21^H A22
            //      = (I - conj(tauQ) a21 a21^H) A22
            Ger( -Conj(tauQ), a21, x12Adj, A22 );

            // Put epsilonQ back instead of the temporary value, 1
            alpha21T.Set(0,0,epsilonQ);
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
L
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R>,STAR,STAR>& tP,
  DistMatrix<Complex<R>,STAR,STAR>& tQ )
{
#ifndef RELEASE
    PushCallStack("bidiag::L");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
        throw std::logic_error
        ("{A,tP,tQ} must be distributed over the same grid");
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
    if( tP.Viewing() || tQ.Viewing() )
        throw std::logic_error("tP and tQ must not be views");
#endif
    typedef Complex<R> C;

    const Grid& g = A.Grid();
    const int tPHeight = std::max(A.Height()-1,0);
    const int tQHeight = A.Height();
    DistMatrix<C,MD,STAR> tPDiag(g), tQDiag(g);
    tPDiag.AlignWithDiagonal( A, -1 );
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
    DistMatrix<C,STAR,STAR> ABR_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> tP1_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> tQ1_STAR_STAR(g);
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
    while( ATL.Height() < A.Height() )
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
        
        if( A22.Height() > 0 )
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

            bidiag::PanelL
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
            bidiag::LUnb( ABR, tP1, tQ1 );
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

#endif // ifndef LAPACK_BIDIAG_L_HPP
