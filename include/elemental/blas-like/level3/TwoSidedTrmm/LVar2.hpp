/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TWOSIDEDTRMM_LVAR2_HPP
#define ELEM_TWOSIDEDTRMM_LVAR2_HPP

#include ELEM_AXPY_INC
#include ELEM_GEMM_INC
#include ELEM_HEMM_INC
#include ELEM_HER2K_INC
#include ELEM_TRMM_INC
#include ELEM_ZEROS_INC

namespace elem {
namespace internal {

// The only reason a field is required is for the existence of 1/2, which is 
// an artifact of the algorithm...
template<typename F> 
inline void
TwoSidedTrmmLVar2( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TwoSidedTrmmLVar2");
        if( A.Height() != A.Width() )
            LogicError( "A must be square." );
        if( L.Height() != L.Width() )
            LogicError( "Triangular matrices must be square." );
        if( A.Height() != L.Height() )
            LogicError( "A and L must be the same size." );
    )
    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;
    Matrix<F>
        LTL, LTR,  L00, L01, L02,
        LBL, LBR,  L10, L11, L12,
                   L20, L21, L22;

    // Temporary products
    Matrix<F> Y21;

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        //--------------------------------------------------------------------//
        // A10 := L11' A10
        Trmm( LEFT, LOWER, ADJOINT, diag, F(1), L11, A10 );

        // A10 := A10 + L21' A20
        Gemm( ADJOINT, NORMAL, F(1), L21, A20, F(1), A10 );

        // Y21 := A22 L21
        Zeros( Y21, A21.Height(), A21.Width() );
        Hemm( LEFT, LOWER, F(1), A22, L21, F(0), Y21 );

        // A21 := A21 L11
        Trmm( RIGHT, LOWER, NORMAL, diag, F(1), L11, A21 );

        // A21 := A21 + 1/2 Y21
        Axpy( F(1)/F(2), Y21, A21 );

        // A11 := L11' A11 L11
        TwoSidedTrmmLUnb( diag, A11, L11 );

        // A11 := A11 + (A21' L21 + L21' A21)
        Her2k( LOWER, ADJOINT, F(1), A21, L21, F(1), A11 );

        // A21 := A21 + 1/2 Y21
        Axpy( F(1)/F(2), Y21, A21 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
}

template<typename F> 
inline void
TwoSidedTrmmLVar2
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& L )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TwoSidedTrmmLVar2");
        if( A.Height() != A.Width() )
            LogicError( "A must be square." );
        if( L.Height() != L.Width() )
            LogicError( "Triangular matrices must be square." );
        if( A.Height() != L.Height() )
            LogicError( "A and L must be the same size." );
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > L21Adj_STAR_MR(g);
    DistMatrix<F,VC,  STAR> L21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > X10_STAR_MR(g);
    DistMatrix<F,STAR,STAR> X11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> Z21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Z21_MR_STAR(g);
    DistMatrix<F,MR,  MC  > Z21_MR_MC(g);
    DistMatrix<F> Y21(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        A21_VC_STAR.AlignWith( A22 );
        L21_MC_STAR.AlignWith( A20 );
        L21_VC_STAR.AlignWith( A22 );
        L21_VR_STAR.AlignWith( A22 );
        L21Adj_STAR_MR.AlignWith( A22 );
        X10_STAR_MR.AlignWith( A10 );
        Y21.AlignWith( A21 );
        Z21_MC_STAR.AlignWith( A22 );
        Z21_MR_STAR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        // A10 := L11' A10
        L11_STAR_STAR = L11;
        A10_STAR_VR = A10;
        LocalTrmm
        ( LEFT, LOWER, ADJOINT, diag, F(1), L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;

        // A10 := A10 + L21' A20
        L21_MC_STAR = L21;
        LocalGemm( ADJOINT, NORMAL, F(1), L21_MC_STAR, A20, X10_STAR_MR );
        A10.ColSumScatterUpdate( F(1), X10_STAR_MR );

        // Y21 := A22 L21
        L21_VC_STAR = L21_MC_STAR;
        L21_VR_STAR = L21_VC_STAR;
        L21_VR_STAR.AdjointPartialColAllGather( L21Adj_STAR_MR );
        Zeros( Z21_MC_STAR, A21.Height(), A21.Width() );
        Zeros( Z21_MR_STAR, A21.Height(), A21.Width() );
        LocalSymmetricAccumulateLL
        ( ADJOINT, 
          F(1), A22, L21_MC_STAR, L21Adj_STAR_MR, Z21_MC_STAR, Z21_MR_STAR );
        Z21_MR_MC.RowSumScatterFrom( Z21_MR_STAR );
        Y21 = Z21_MR_MC;
        Y21.RowSumScatterUpdate( F(1), Z21_MC_STAR ); 

        // A21 := A21 L11
        A21_VC_STAR = A21;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;

        // A21 := A21 + 1/2 Y21
        Axpy( F(1)/F(2), Y21, A21 );

        // A11 := L11' A11 L11
        A11_STAR_STAR = A11;
        LocalTwoSidedTrmm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A11 := A11 + (A21' L21 + L21' A21)
        A21_VC_STAR = A21;
        Zeros( X11_STAR_STAR, A11.Height(), A11.Width() );
        Her2k
        ( LOWER, ADJOINT,
          F(1), A21_VC_STAR.Matrix(), L21_VC_STAR.Matrix(),
          F(0), X11_STAR_STAR.Matrix() );
        A11.SumScatterUpdate( F(1), X11_STAR_STAR );

        // A21 := A21 + 1/2 Y21
        Axpy( F(1)/F(2), Y21, A21 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TWOSIDEDTRMM_LVAR2_HPP
