/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_BIDIAG_LUNB_HPP
#define LAPACK_BIDIAG_LUNB_HPP

#include "elemental/blas-like/level1/Conjugate.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace bidiag {

template<typename F> 
inline void LUnb
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& tP, DistMatrix<F,MD,STAR>& tQ )
{
#ifndef RELEASE
    CallStackEntry entry("bidiag::LUnb");
#endif
    const int tPHeight = A.Height();
    const int tQHeight = std::max(A.Height()-1,0);
#ifndef RELEASE
    if( A.Grid() != tP.Grid() || tP.Grid() != tQ.Grid() )
        throw std::logic_error("Process grids do not match");
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
    if( tP.Viewing() && (tP.Height() != tPHeight || tP.Width() != 1) )
        throw std::logic_error("tP is the wrong height");
    if( tQ.Viewing() && (tQ.Height() != tQHeight || tQ.Width() != 1) )
        throw std::logic_error("tQ is the wrong height");
#endif
    const Grid& g = A.Grid();

    if( !tP.Viewing() )
        tP.ResizeTo( tPHeight, 1 );
    if( !tQ.Viewing() )
        tQ.ResizeTo( tQHeight, 1 );

    // Matrix views 
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  alpha21T(g), a1R(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  a21B(g),     A2R(g),
                         A20(g), a21(g),     A22(g);

    // Temporary matrices
    DistMatrix<F,MC,  STAR> a21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > a1R_STAR_MR(g);
    DistMatrix<F,MR,  STAR> x12Adj_MR_STAR(g);
    DistMatrix<F,MC,  STAR> w21_MC_STAR(g);

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

        View1x2( a1R, alpha11, a12 );
        View1x2( A2R, a21, A22 );

        a21_MC_STAR.AlignWith( A22 );
        a1R_STAR_MR.AlignWith( A2R );
        x12Adj_MR_STAR.AlignWith( A22 );
        w21_MC_STAR.AlignWith( A2R );

        const bool thisIsMyRow = ( g.Row() == alpha11.ColAlignment() );
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        const bool nextIsMyRow = ( g.Row() == a21.ColAlignment() );
        //--------------------------------------------------------------------//

        // Due to deficiencies in the BLAS ?gemv routines, this section is
        // easier if we temporary conjugate a1R = | alpha11, a12 |
        Conjugate( a1R );

        // Find tauP, u, and epsilonP such that
        //     I - conj(tauP) | 1 | | 1, v^H | | alpha11 | = | epsilonP |
        //                    | v |            |   a12^T | = |    0     |
        const F tauP = Reflector( alpha11, a12 );
        tP.Set(A00.Height(),0,tauP);
        F epsilonP=0;
        if( thisIsMyCol && thisIsMyRow )
            epsilonP = alpha11.GetLocal(0,0);

        // Set a1R^T = | 1 | and form w21 := A2R a1R^T = A2R | 1 |
        //             | v |                                 | v |
        alpha11.Set(0,0,F(1));
        a1R_STAR_MR = a1R;
        Zeros( w21_MC_STAR, a21.Height(), 1 );
        LocalGemv( NORMAL, F(1), A2R, a1R_STAR_MR, F(0), w21_MC_STAR );
        w21_MC_STAR.SumOverRow();

        // A2R := A2R - tauP w21 conj(a1R)
        //      = A2R - tauP A2R a1R^T conj(a1R)
        //      = A2R conj(I - conj(tauP) a1R^H a1R)
        // which compensates for the fact that the reflector was generated
        // on the conjugated a1R
        LocalGer( -tauP, w21_MC_STAR, a1R_STAR_MR, A2R );

        // Put epsilonP back instead of the temporary value, 1
        if( thisIsMyCol && thisIsMyRow )
            alpha11.SetLocal(0,0,epsilonP);

        // Undo the temporary conjugation
        Conjugate( a1R );

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            PartitionDown
            ( a21, alpha21T,
                   a21B, 1 );

            // Find tauQ, u, and epsilonQ such that
            //     I - conj(tauQ) | 1 | | 1, u^H | | alpha21T | = | epsilonQ |
            //                    | u |            | a21B     | = |    0     |
            const F tauQ = Reflector( alpha21T, a21B );
            tQ.Set(A00.Height(),0,tauQ);
            F epsilonQ=0;
            if( nextIsMyRow && thisIsMyCol )
                epsilonQ = alpha21T.GetLocal(0,0);

            // Set a21 = | 1 | and form x12^H = (a21^H A22)^H = A22^H a21
            //           | u |  
            alpha21T.Set(0,0,F(1));
            a21_MC_STAR = a21;
            Zeros( x12Adj_MR_STAR, a12.Width(), 1 );
            LocalGemv( ADJOINT, F(1), A22, a21_MC_STAR, F(0), x12Adj_MR_STAR );
            x12Adj_MR_STAR.SumOverCol();

            // A22 := A22 - conj(tauQ) a21 x12
            //      = A22 - conj(tauQ) a21 a21^H A22
            //      = (I - conj(tauQ) a21 a21^H) A22
            LocalGer( -Conj(tauQ), a21_MC_STAR, x12Adj_MR_STAR, A22 );

            // Put epsilonQ back instead of the temporary value, 1
            if( nextIsMyRow && thisIsMyCol )
                alpha21T.SetLocal(0,0,epsilonQ);
        }
        //--------------------------------------------------------------------//
        a21_MC_STAR.FreeAlignments();
        a1R_STAR_MR.FreeAlignments();
        x12Adj_MR_STAR.FreeAlignments();
        w21_MC_STAR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
}

} // namespace bidiag
} // namespace elem

#endif // ifndef LAPACK_BIDIAG_LUNB_HPP
