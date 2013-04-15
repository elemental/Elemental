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

template<typename R>
inline void LUnb( DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("bidiag::LUnb");
    if( A.Height() > A.Width() )
        throw std::logic_error("A must be at least as wide as it is tall");
#endif
    const Grid& g = A.Grid();

    // Matrix views 
    DistMatrix<R>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  alpha21T(g), a1R(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  a21B(g),     A2R(g),
                         A20(g), a21(g),     A22(g);

    // Temporary matrices
    DistMatrix<R,MC,  STAR> a21_MC_STAR(g);
    DistMatrix<R,STAR,MR  > a1R_STAR_MR(g);
    DistMatrix<R,MR,  STAR> x12Trans_MR_STAR(g);
    DistMatrix<R,MC,  STAR> w21_MC_STAR(g);

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

        View1x2( a1R, alpha11, a12 );
        View1x2( A2R, a21, A22 );

        a21_MC_STAR.AlignWith( A22 );
        a1R_STAR_MR.AlignWith( A2R );
        x12Trans_MR_STAR.AlignWith( A22 );
        w21_MC_STAR.AlignWith( A2R );
        Zeros( a12.Width(), 1, x12Trans_MR_STAR );
        Zeros( a21.Height(), 1, w21_MC_STAR );
        const bool thisIsMyRow = ( g.Row() == alpha11.ColAlignment() );
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        const bool nextIsMyRow = ( g.Row() == a21.ColAlignment() );
        //--------------------------------------------------------------------//

        // Find tauP, u, and epsilonP such that
        //     I - tauP | 1 | | 1, v^T | | alpha11 | = | epsilonP |
        //              | v |            |   a12^T | = |    0     |
        const R tauP = Reflector( alpha11, a12 );
        R epsilonP=0;
        if( thisIsMyCol && thisIsMyRow )
            epsilonP = alpha11.GetLocal(0,0);

        // Set a1R^T = | 1 | and form w21 := A2R a1R^T = A2R | 1 |
        //             | v |                                 | v |
        alpha11.Set(0,0,R(1));
        a1R_STAR_MR = a1R;
        LocalGemv( NORMAL, R(1), A2R, a1R_STAR_MR, R(0), w21_MC_STAR );
        w21_MC_STAR.SumOverRow();

        // A2R := A2R - tauP w21 a1R
        //      = A2R - tauP A2R a1R^T a1R
        //      = A2R (I - tauP a1R^T a1R)
        LocalGer( -tauP, w21_MC_STAR, a1R_STAR_MR, A2R );

        // Put epsilonP back instead of the temporary value, 1
        if( thisIsMyCol && thisIsMyRow )
            alpha11.SetLocal(0,0,epsilonP);

        if( A22.Height() != 0 )
        {
            // Expose the subvector we seek to zero, a21B
            PartitionDown
            ( a21, alpha21T,
                   a21B );

            // Find tauQ, u, and epsilonQ such that
            //     I - tauQ | 1 | | 1, u^T | | alpha21T | = | epsilonQ |
            //              | u |            | a21B     | = |    0     |
            const R tauQ = Reflector( alpha21T, a21B );
            R epsilonQ=0;
            if( nextIsMyRow && thisIsMyCol )
                epsilonQ = alpha21T.GetLocal(0,0);

            // Set a21 = | 1 | and form x12^T = (a21^T A22)^T = A22^T a21
            //           | u |  
            alpha21T.Set(0,0,R(1));
            a21_MC_STAR = a21;
            LocalGemv
            ( TRANSPOSE, R(1), A22, a21_MC_STAR, R(0), x12Trans_MR_STAR );
            x12Trans_MR_STAR.SumOverCol();

            // A22 := A22 - tauQ a21 x12
            //      = A22 - tauQ a21 a21^T A22
            //      = (I - tauQ a21 a21^T) A22
            LocalGer( -tauQ, a21_MC_STAR, x12Trans_MR_STAR, A22 );

            // Put epsilonQ back instead of the temporary value, 1
            if( nextIsMyRow && thisIsMyCol )
                alpha21T.SetLocal(0,0,epsilonQ);
        }
        //--------------------------------------------------------------------//
        a1R_STAR_MR.FreeAlignments();
        a21_MC_STAR.FreeAlignments();
        x12Trans_MR_STAR.FreeAlignments();
        w21_MC_STAR.FreeAlignments();

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
inline void LUnb
( DistMatrix<Complex<R> >& A, 
  DistMatrix<Complex<R>,MD,STAR>& tP,
  DistMatrix<Complex<R>,MD,STAR>& tQ )
{
#ifndef RELEASE
    PushCallStack("bidiag::LUnb");
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
    typedef Complex<R> C;
    const Grid& g = A.Grid();

    if( !tP.Viewing() )
        tP.ResizeTo( tPHeight, 1 );
    if( !tQ.Viewing() )
        tQ.ResizeTo( tQHeight, 1 );

    // Matrix views 
    DistMatrix<C>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  alpha21T(g), a1R(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  a21B(g),     A2R(g),
                         A20(g), a21(g),     A22(g);

    // Temporary matrices
    DistMatrix<C,MC,  STAR> a21_MC_STAR(g);
    DistMatrix<C,STAR,MR  > a1R_STAR_MR(g);
    DistMatrix<C,MR,  STAR> x12Adj_MR_STAR(g);
    DistMatrix<C,MC,  STAR> w21_MC_STAR(g);

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

        View1x2( a1R, alpha11, a12 );
        View1x2( A2R, a21, A22 );

        a21_MC_STAR.AlignWith( A22 );
        a1R_STAR_MR.AlignWith( A2R );
        x12Adj_MR_STAR.AlignWith( A22 );
        w21_MC_STAR.AlignWith( A2R );
        Zeros( a12.Width(), 1, x12Adj_MR_STAR );
        Zeros( a21.Height(), 1, w21_MC_STAR );
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
        const C tauP = Reflector( alpha11, a12 );
        tP.Set(A00.Height(),0,tauP);
        C epsilonP=0;
        if( thisIsMyCol && thisIsMyRow )
            epsilonP = alpha11.GetLocal(0,0);

        // Set a1R^T = | 1 | and form w21 := A2R a1R^T = A2R | 1 |
        //             | v |                                 | v |
        alpha11.Set(0,0,C(1));
        a1R_STAR_MR = a1R;
        LocalGemv( NORMAL, C(1), A2R, a1R_STAR_MR, C(0), w21_MC_STAR );
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
                   a21B );

            // Find tauQ, u, and epsilonQ such that
            //     I - conj(tauQ) | 1 | | 1, u^H | | alpha21T | = | epsilonQ |
            //                    | u |            | a21B     | = |    0     |
            const C tauQ = Reflector( alpha21T, a21B );
            tQ.Set(A00.Height(),0,tauQ);
            C epsilonQ=0;
            if( nextIsMyRow && thisIsMyCol )
                epsilonQ = alpha21T.GetLocal(0,0);

            // Set a21 = | 1 | and form x12^H = (a21^H A22)^H = A22^H a21
            //           | u |  
            alpha21T.Set(0,0,C(1));
            a21_MC_STAR = a21;
            LocalGemv( ADJOINT, C(1), A22, a21_MC_STAR, C(0), x12Adj_MR_STAR );
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
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace bidiag
} // namespace elem

#endif // ifndef LAPACK_BIDIAG_LUNB_HPP
