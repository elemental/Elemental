/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_CHOLESKY_LVAR2_HPP
#define LAPACK_CHOLESKY_LVAR2_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"

#include "elemental/lapack-like/Cholesky/LVar3.hpp"

// TODO: Reverse variants

namespace elem {
namespace cholesky {

template<typename F> 
inline void
LVar2( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::LVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
#endif
    // Matrix views
    Matrix<F> 
        ATL, ATR,   A00, A01, A02,
        ABL, ABR,   A10, A11, A12,
                    A20, A21, A22;

    // Start the algorithm
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

        //--------------------------------------------------------------------//
        Herk( LOWER, NORMAL, F(-1), A10, F(1), A11 );
        cholesky::LVar3Unb( A11 );
        Gemm( NORMAL, ADJOINT, F(-1), A20, A10, F(1), A21 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11, A21 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
}

template<typename F> 
inline void
LVar2( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("cholesky::LVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error
        ("Can only compute Cholesky factor of square matrices");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F> 
        ATL(g), ATR(g),   A00(g), A01(g), A02(g),
        ABL(g), ABR(g),   A10(g), A11(g), A12(g),
                          A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,MR,  STAR> A10Adj_MR_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,MC,  STAR> X11_MC_STAR(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);

    // Start the algorithm
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

        A10Adj_MR_STAR.AlignWith( A10 );
        X11_MC_STAR.AlignWith( A10 );
        X21_MC_STAR.AlignWith( A20 );
        //--------------------------------------------------------------------//
        A10Adj_MR_STAR.AdjointFrom( A10 );
        LocalGemm( NORMAL, NORMAL, F(1), A10, A10Adj_MR_STAR, X11_MC_STAR );
        A11.SumScatterUpdate( F(-1), X11_MC_STAR );

        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );
        A11 = A11_STAR_STAR;

        LocalGemm( NORMAL, NORMAL, F(1), A20, A10Adj_MR_STAR, X21_MC_STAR );
        A21.SumScatterUpdate( F(-1), X21_MC_STAR );

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A10Adj_MR_STAR.FreeAlignments();
        X11_MC_STAR.FreeAlignments();
        X21_MC_STAR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
}

} // namespace cholesky
} // namespace elem

#endif // ifndef LAPACK_CHOLESKY_LVAR2_HPP
