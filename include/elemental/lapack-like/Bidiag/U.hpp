/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_BIDIAG_U_HPP
#define ELEM_LAPACK_BIDIAG_U_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/Bidiag/PanelU.hpp"
#include "elemental/lapack-like/Bidiag/UUnb.hpp"
#include "elemental/lapack-like/Reflector.hpp"

namespace elem {
namespace bidiag {

template<typename F>
inline void U( Matrix<F>& A, Matrix<F>& tP, Matrix<F>& tQ )
{
#ifndef RELEASE
    CallStackEntry entry("bidiag::U");
#endif
    const Int tPHeight = Max(A.Width()-1,0);
    const Int tQHeight = A.Width();
#ifndef RELEASE
    if( A.Height() < A.Width() )
        LogicError("A must be at least as tall as it is wide");
#endif
    tP.ResizeTo( tPHeight, 1 );
    tQ.ResizeTo( tQHeight, 1 );

    // Matrix views 
    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  alpha12L, a12R,
        ABL, ABR,  a10, alpha11, a12,  aB1, AB2,
                   A20, a21,     A22;

    // Temporary matrices
    Matrix<F> x12Adj, w21;

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

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
        const F tauQ = Reflector( alpha11, a21 );
        const F epsilonQ = alpha11.Get(0,0);
        tQ.Set(A00.Height(),0,tauQ );

        // Set aB1 = | 1 | and form x12^H := (aB1^H AB2)^H = AB2^H aB1
        //           | u |
        alpha11.Set(0,0,F(1));
        Gemv( ADJOINT, F(1), AB2, aB1, F(0), x12Adj );

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
            PartitionRight( a12, alpha12L, a12R, 1 );

            // Find tauP, v, and epsilonP such that
            //     I - conj(tauP) | 1 | | 1, v^H | | alpha12L | = | epsilonP |
            //                    | v |            |  a12R^T  |   |    0     |
            const F tauP = Reflector( alpha12L, a12R );
            const F epsilonP = alpha12L.Get(0,0);
            tP.Set(A00.Height(),0,tauP);

            // Set a12^T = | 1 | and form w21 := A22 a12^T = A22 | 1 |
            //             | v |                                 | v |
            alpha12L.Set(0,0,F(1));
            Gemv( NORMAL, F(1), A22, a12, F(0), w21 );

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
}

template<typename F> 
inline void
U( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& tP, DistMatrix<F,STAR,STAR>& tQ )
{
#ifndef RELEASE
    CallStackEntry entry("bidiag::U");
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
        LogicError
        ("{A,tP,tQ} must be distributed over the same grid");
    if( A.Height() < A.Width() )
        LogicError("A must be at least as tall as it is wide");
    // Are these requirements necessary?!?
    if( tP.Viewing() || tQ.Viewing() )
        LogicError("tP and tQ must not be views");
#endif
    const Grid& g = A.Grid();
    const Int tPHeight = Max(A.Width()-1,0);
    const Int tQHeight = A.Width();
    DistMatrix<F,MD,STAR> tPDiag(g), tQDiag(g);
    tPDiag.AlignWithDiagonal( A, 1 );
    tQDiag.AlignWithDiagonal( A, 0 );
    tPDiag.ResizeTo( tPHeight, 1 );
    tQDiag.ResizeTo( tQHeight, 1 );

    // Matrix views 
    DistMatrix<F> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F,MD,STAR> tPT(g),  tP0(g), 
                          tPB(g),  tP1(g),
                                   tP2(g);
    DistMatrix<F,MD,STAR> tQT(g),  tQ0(g), 
                          tQB(g),  tQ1(g),
                                   tQ2(g);

    // Temporary distributions
    DistMatrix<F> X(g), X11(g),
                        X21(g);
    DistMatrix<F> Y(g), Y11(g),
                        Y21(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Y21_MR_STAR(g);
    DistMatrix<F,MC,  STAR> AColPan_MC_STAR(g), A11_MC_STAR(g),
                                                A21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > ARowPan_STAR_MR(g), A11_STAR_MR(g), A12_STAR_MR(g);

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
            ( NORMAL, ADJOINT, F(-1), A21_MC_STAR, Y21_MR_STAR, F(1), A22 );
            LocalGemm
            ( NORMAL, NORMAL, F(-1), X21_MC_STAR, A12_STAR_MR, F(1), A22 );
            //----------------------------------------------------------------//
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
}

} // namespace bidiag
} // namespace elem

#endif // ifndef ELEM_LAPACK_BIDIAG_U_HPP
