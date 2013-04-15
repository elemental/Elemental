/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HPDINVERSE_CHOLESKYLVAR2_HPP
#define LAPACK_HPDINVERSE_CHOLESKYLVAR2_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/blas-like/level3/Trtrmm.hpp"
#include "elemental/lapack-like/Cholesky.hpp"
#include "elemental/lapack-like/TriangularInverse.hpp"

namespace elem {
namespace hpd_inverse {

//
// This approach is based upon a (conjugate)-transposition of the reordered 
// Variant 2 algorithm from Fig. 9 in Bientinesi et al.'s "Families of 
// Algorithms Related to the Inversion of a Symmetric Positive Definite Matrix".
//

template<typename F> 
inline void
CholeskyLVar2( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("hpd_inverse::CholeskyLVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    // Matrix views
    Matrix<F> 
        ATL, ATR,  A00, A01, A02,
        ABL, ABR,  A10, A11, A12,
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
        Cholesky( LOWER, A11 );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A11, A10 );
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11, A21 );
        Herk( LOWER, ADJOINT, F(1), A10, F(1), A00 );
        Gemm( NORMAL, NORMAL, F(-1), A21, A10, F(1), A20 );
        Herk( LOWER, NORMAL, F(-1), A21, F(1), A22 );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A11, A10 );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(-1), A11, A21 );
        TriangularInverse( LOWER, NON_UNIT, A11 );
        Trtrmm( ADJOINT, LOWER, A11 );
        //--------------------------------------------------------------------//

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

template<typename F> 
inline void
CholeskyLVar2( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("hpd_inverse::CholeskyLVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error("Nonsquare matrices cannot be triangular");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,STAR,MC  > A10_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);

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

        A10_STAR_VR.AlignWith( A00 );
        A21_VC_STAR.AlignWith( A20 );
        A10_STAR_MC.AlignWith( A00 );
        A10_STAR_MR.AlignWith( A00 );
        A21Trans_STAR_MC.AlignWith( A20 );
        A21_VR_STAR.AlignWith( A22 );
        A21Adj_STAR_MR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        LocalCholesky( LOWER, A11_STAR_STAR );

        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A10_STAR_VR );

        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A21_VC_STAR );

        A10_STAR_MC = A10_STAR_VR;
        A10_STAR_MR = A10_STAR_VR;
        LocalTrrk
        ( LOWER, ADJOINT,
          F(1), A10_STAR_MC, A10_STAR_MR, F(1), A00 );

        A21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        LocalGemm
        ( TRANSPOSE, NORMAL, F(-1), A21Trans_STAR_MC, A10_STAR_MR, F(1), A20 );

        A21_VR_STAR = A21_VC_STAR;
        A21Adj_STAR_MR.AdjointFrom( A21_VR_STAR );
        LocalTrrk
        ( LOWER, TRANSPOSE,
          F(-1), A21Trans_STAR_MC, A21Adj_STAR_MR, F(1), A22 );

        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A10_STAR_VR );

        LocalTrsm
        ( RIGHT, LOWER, NORMAL, NON_UNIT, F(-1), A11_STAR_STAR, A21_VC_STAR );

        LocalTriangularInverse( LOWER, NON_UNIT, A11_STAR_STAR );

        LocalTrtrmm( ADJOINT, LOWER, A11_STAR_STAR );

        A11 = A11_STAR_STAR;
        A10 = A10_STAR_VR;
        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A10_STAR_VR.FreeAlignments();
        A21_VC_STAR.FreeAlignments();
        A10_STAR_MC.FreeAlignments();
        A10_STAR_MR.FreeAlignments();
        A21Trans_STAR_MC.FreeAlignments();
        A21_VR_STAR.FreeAlignments();
        A21Adj_STAR_MR.FreeAlignments();

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

} // namespace hpd_inverse
} // namespace elem

#endif // ifndef LAPACK_HPDINVERSE_CHOLESKYLVAR2_HPP
