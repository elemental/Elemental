/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_HPDINVERSE_CHOLESKYUVAR2_HPP
#define LAPACK_HPDINVERSE_CHOLESKYUVAR2_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/blas-like/level3/Trtrmm.hpp"
#include "elemental/lapack-like/Cholesky.hpp"
#include "elemental/lapack-like/TriangularInverse.hpp"

namespace elem {
namespace hpd_inverse {

//
// This approach is based upon the reordered Variant 2 algorithm from Fig. 9 in 
// Bientinesi et al.'s "Families of Algorithms Related to the Inversion of 
// a Symmetric Positive Definite Matrix".
//

template<typename F> 
inline void
CholeskyUVar2( Matrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("hpd_inverse::CholeskyUVar2");
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
        Cholesky( UPPER, A11 );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11, A01 );
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11, A12 );
        Herk( UPPER, NORMAL, F(1), A01, F(1), A00 );
        Gemm( NORMAL, NORMAL, F(-1), A01, A12, F(1), A02 );
        Herk( UPPER, ADJOINT, F(-1), A12, F(1), A22 );
        Trsm( RIGHT, UPPER, ADJOINT, NON_UNIT, F(1), A11, A01 );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(-1), A11, A12 );
        TriangularInverse( UPPER, NON_UNIT, A11 );
        Trtrmm( ADJOINT, UPPER, A11 );
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
CholeskyUVar2( DistMatrix<F>& A )
{
#ifndef RELEASE
    PushCallStack("hpd_inverse::CholeskyUVar2");
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
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A01_VR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,MC  > A01Trans_STAR_MC(g);
    DistMatrix<F,MR,  STAR> A01_MR_STAR(g);
    DistMatrix<F,STAR,MR  > A01Adj_STAR_MR(g);
    DistMatrix<F,STAR,MR  > A12_STAR_MR(g);
    DistMatrix<F,STAR,MC  > A12_STAR_MC(g);

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

        A01_VC_STAR.AlignWith( A00 );
        A12_STAR_VR.AlignWith( A02 );
        A01Trans_STAR_MC.AlignWith( A00 );
        A01_VR_STAR.AlignWith( A00 );
        A01Adj_STAR_MR.AlignWith( A00 );
        A12_STAR_MR.AlignWith( A02 );
        A12_STAR_MC.AlignWith( A22 );
        //--------------------------------------------------------------------//
        A11_STAR_STAR = A11;
        LocalCholesky( UPPER, A11_STAR_STAR );

        A01_VC_STAR = A01;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), A11_STAR_STAR, A01_VC_STAR );

        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A12_STAR_VR );

        A01Trans_STAR_MC.TransposeFrom( A01_VC_STAR );
        A01_VR_STAR = A01_VC_STAR;
        A01Adj_STAR_MR.AdjointFrom( A01_VR_STAR );
        LocalTrrk
        ( UPPER, TRANSPOSE,
          F(1), A01Trans_STAR_MC, A01Adj_STAR_MR, F(1), A00 );

        A12_STAR_MR = A12_STAR_VR;
        LocalGemm
        ( TRANSPOSE, NORMAL, F(-1), A01Trans_STAR_MC, A12_STAR_MR, F(1), A02 );

        A12_STAR_MC = A12_STAR_VR;
        LocalTrrk
        ( UPPER, ADJOINT,
          F(-1), A12_STAR_MC, A12_STAR_MR, F(1), A22 );

        LocalTrsm
        ( RIGHT, UPPER, ADJOINT, NON_UNIT, F(1), A11_STAR_STAR, A01_VC_STAR );

        LocalTrsm
        ( LEFT, UPPER, NORMAL, NON_UNIT, F(-1), A11_STAR_STAR, A12_STAR_VR );

        LocalTriangularInverse( UPPER, NON_UNIT, A11_STAR_STAR );

        LocalTrtrmm( ADJOINT, UPPER, A11_STAR_STAR );

        A11 = A11_STAR_STAR;
        A01 = A01_VC_STAR;
        A12 = A12_STAR_VR;
        //--------------------------------------------------------------------//
        A01_VC_STAR.FreeAlignments();
        A12_STAR_VR.FreeAlignments();
        A01Trans_STAR_MC.FreeAlignments();
        A01_VR_STAR.FreeAlignments();
        A01Adj_STAR_MR.FreeAlignments();
        A12_STAR_MR.FreeAlignments();
        A12_STAR_MC.FreeAlignments();

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

#endif // ifndef LAPACK_HPDINVERSE_CHOLESKYUVAR2_HPP
