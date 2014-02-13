/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TWOSIDEDTRMM_UVAR2_HPP
#define ELEM_TWOSIDEDTRMM_UVAR2_HPP

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
TwoSidedTrmmUVar2( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TwoSidedTrmmUVar2");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( U.Height() != U.Width() )
            LogicError("Triangular matrices must be square");
        if( A.Height() != U.Height() )
            LogicError("A and U must be the same size");
    )
    // Matrix views
    Matrix<F>
        ATL, ATR,  A00, A01, A02,
        ABL, ABR,  A10, A11, A12,
                   A20, A21, A22;
    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;

    // Temporary products
    Matrix<F> Y12;

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        //--------------------------------------------------------------------//
        // A01 := A01 U11'
        Trmm( RIGHT, UPPER, ADJOINT, diag, F(1), U11, A01 );

        // A01 := A01 + A02 U12'
        Gemm( NORMAL, ADJOINT, F(1), A02, U12, F(1), A01 );

        // Y12 := U12 A22
        Zeros( Y12, A12.Height(), A12.Width() );
        Hemm( RIGHT, UPPER, F(1), A22, U12, F(0), Y12 );

        // A12 := U11 A12
        Trmm( LEFT, UPPER, NORMAL, diag, F(1), U11, A12 );

        // A12 := A12 + 1/2 Y12
        Axpy( F(1)/F(2), Y12, A12 );

        // A11 := U11 A11 U11'
        TwoSidedTrmmUUnb( diag, A11, U11 );

        // A11 := A11 + (A12 U12' + U12 A12')
        Her2k( UPPER, NORMAL, F(1), A12, U12, F(1), A11 );

        // A12 := A12 + 1/2 Y12
        Axpy( F(1)/F(2), Y12, A12 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
}

template<typename F> 
inline void
TwoSidedTrmmUVar2
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& U )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TwoSidedTrmmUVar2");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( U.Height() != U.Width() )
            LogicError("Triangular matrices must be square");
        if( A.Height() != U.Height() )
            LogicError("A and U must be the same size");
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,STAR,VR  > U12_STAR_VR(g);
    DistMatrix<F,MR,  STAR> U12Adj_MR_STAR(g);
    DistMatrix<F,VC,  STAR> U12Adj_VC_STAR(g);
    DistMatrix<F,MC,  STAR> X01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> X11_STAR_STAR(g);
    DistMatrix<F,MR,  MC  > Z12Adj_MR_MC(g);
    DistMatrix<F,MC,  STAR> Z12Adj_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Z12Adj_MR_STAR(g);
    DistMatrix<F> Y12(g);
    DistMatrix<F> Z12Adj(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        A12_STAR_VR.AlignWith( A12 );
        U12_STAR_MC.AlignWith( A22 );
        U12_STAR_VR.AlignWith( A12 );
        U12Adj_MR_STAR.AlignWith( A22 );
        U12Adj_VC_STAR.AlignWith( A22 );
        X01_MC_STAR.AlignWith( A01 );
        Y12.AlignWith( A12 );
        Z12Adj.AlignWith( A12 );
        Z12Adj_MR_MC.AlignWith( A12 );
        Z12Adj_MC_STAR.AlignWith( A22 );
        Z12Adj_MR_STAR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        // A01 := A01 U11'
        U11_STAR_STAR = U11;
        A01_VC_STAR = A01;
        LocalTrmm
        ( RIGHT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // A01 := A01 + A02 U12'
        U12.AdjointColAllGather( U12Adj_MR_STAR );
        LocalGemm( NORMAL, NORMAL, F(1), A02, U12Adj_MR_STAR, X01_MC_STAR );
        A01.RowSumScatterUpdate( F(1), X01_MC_STAR );

        // Y12 := U12 A22
        U12Adj_VC_STAR = U12Adj_MR_STAR;
        U12Adj_VC_STAR.AdjointPartialColAllGather( U12_STAR_MC );
        Zeros( Z12Adj_MC_STAR, A12.Width(), A12.Height() );
        Zeros( Z12Adj_MR_STAR, A12.Width(), A12.Height() );
        LocalSymmetricAccumulateRU
        ( ADJOINT, 
          F(1), A22, U12_STAR_MC, U12Adj_MR_STAR, 
          Z12Adj_MC_STAR, Z12Adj_MR_STAR );
        Z12Adj.RowSumScatterFrom( Z12Adj_MC_STAR );
        Z12Adj_MR_MC = Z12Adj;
        Z12Adj_MR_MC.RowSumScatterUpdate( F(1), Z12Adj_MR_STAR );
        Y12.Resize( A12.Height(), A12.Width() );
        Adjoint( Z12Adj_MR_MC.LockedMatrix(), Y12.Matrix() );

        // A12 := U11 A12
        A12_STAR_VR = A12;
        U11_STAR_STAR = U11;
        LocalTrmm
        ( LEFT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;

        // A12 := A12 + 1/2 Y12
        Axpy( F(1)/F(2), Y12, A12 );

        // A11 := U11 A11 U11'
        A11_STAR_STAR = A11;
        LocalTwoSidedTrmm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A11 := A11 + (A12 U12' + U12 A12')
        A12_STAR_VR = A12;
        U12_STAR_VR = U12;
        Zeros( X11_STAR_STAR, A11.Height(), A11.Width() );
        Her2k
        ( UPPER, NORMAL,
          F(1), A12_STAR_VR.Matrix(), U12_STAR_VR.Matrix(),
          F(0), X11_STAR_STAR.Matrix() );
        A11.SumScatterUpdate( F(1), X11_STAR_STAR );

        // A12 := A12 + 1/2 Y12
        Axpy( F(1)/F(2), Y12, A12 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TWOSIDEDTRMM_UVAR2_HPP
