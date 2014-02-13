/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TWOSIDEDTRSM_UVAR3_HPP
#define ELEM_TWOSIDEDTRSM_UVAR3_HPP

#include ELEM_AXPY_INC
#include ELEM_MAKEHERMITIAN_INC
#include ELEM_MAKETRIANGULAR_INC
#include ELEM_GEMM_INC
#include ELEM_HEMM_INC
#include ELEM_HER2K_INC
#include ELEM_TRSM_INC
#include ELEM_ZEROS_INC

namespace elem {
namespace internal {

template<typename F> 
inline void
TwoSidedTrsmUVar3( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TwoSidedTrsmUVar4");
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
        YTL, YTR,  Y00, Y01, Y02,
        YBL, YBR,  Y10, Y11, Y12,
                   Y20, Y21, Y22;
    Matrix<F>
        UTL, UTR,  U00, U01, U02,
        UBL, UBR,  U10, U11, U12,
                   U20, U21, U22;

    // We will use an entire extra matrix as temporary storage. If this is not
    // acceptable, use TwoSidedTrsmUVar4 instead.
    Matrix<F> Y;
    Zeros( Y, A.Height(), A.Width() );

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDownDiagonal
    ( Y, YTL, YTR,
         YBL, YBR, 0 );
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

        RepartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, /**/ Y01, Y02,
         /*************/ /******************/
               /**/       Y10, /**/ Y11, Y12,
          YBL, /**/ YBR,  Y20, /**/ Y21, Y22 );

        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        //--------------------------------------------------------------------//
        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A11 := A11 - (A01' U01 + U01' A01)
        Her2k( UPPER, ADJOINT, F(-1), A01, U01, F(1), A11 );

        // A11 := inv(U11)' A11 inv(U11)
        TwoSidedTrsmUUnb( diag, A11, U11 );

        // A12 := A12 - U01' A02
        Gemm( ADJOINT, NORMAL, F(-1), U01, A02, F(1), A12 );

        // A12 := inv(U11)' A12
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U11, A12 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        Trsm( RIGHT, UPPER, NORMAL, diag, F(1), U11, A01 );

        // Y02 := Y02 + A01 U12
        Gemm( NORMAL, NORMAL, F(1), A01, U12, F(1), Y02 );

        // Y12 := Y12 + A11 U12
        Hemm( LEFT, UPPER, F(1), A11, U12, F(1), Y12 );

        // Y12 := Y12 + A01' U02
        Gemm( ADJOINT, NORMAL, F(1), A01, U02, F(1), Y12 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, Y01, /**/ Y02,
               /**/       Y10, Y11, /**/ Y12,
         /*************/ /******************/
          YBL, /**/ YBR,  Y20, Y21, /**/ Y22 );

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /**********************************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
}

template<typename F> 
inline void
TwoSidedTrsmUVar3
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& U )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TwoSidedTrsmUVar4");
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
        YTL(g), YTR(g),  Y00(g), Y01(g), Y02(g),
        YBL(g), YBR(g),  Y10(g), Y11(g), Y12(g),
                         Y20(g), Y21(g), Y22(g);
    DistMatrix<F>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions
    DistMatrix<F,MC,  STAR> A11_MC_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,MR,  STAR> U12Adj_MR_STAR(g);
    DistMatrix<F,STAR,STAR> X11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X12_STAR_MR(g);
    DistMatrix<F,STAR,MR  > Z12_STAR_MR(g);

    // We will use an entire extra matrix as temporary storage. If this is not
    // acceptable, use TwoSidedTrsmUVar4 instead.
    DistMatrix<F> Y(g);
    Y.AlignWith( A );
    Zeros( Y, A.Height(), A.Width() );

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDownDiagonal
    ( Y, YTL, YTR,
         YBL, YBR, 0 );
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

        RepartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, /**/ Y01, Y02,
         /*************/ /******************/
               /**/       Y10, /**/ Y11, Y12,
          YBL, /**/ YBR,  Y20, /**/ Y21, Y22 );

        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        A11_MC_STAR.AlignWith( Y12 );
        A12_STAR_VR.AlignWith( A12 );
        A01_VC_STAR.AlignWith( A01 );
        A01_MC_STAR.AlignWith( A01 );
        U01_VC_STAR.AlignWith( A01 );
        U01_MC_STAR.AlignWith( A01 );
        U12Adj_MR_STAR.AlignWith( Y12 );
        X12_STAR_MR.AlignWith( A02 );
        Z12_STAR_MR.AlignWith( U02 );
        //--------------------------------------------------------------------//
        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A11 := A11 - (A01' U01 + U01' A01)
        A01_VC_STAR = A01;
        U01_VC_STAR = U01;
        Zeros( X11_STAR_STAR, A11.Height(), A11.Width() );
        Her2k
        ( UPPER, ADJOINT, 
          F(1), A01_VC_STAR.Matrix(), U01_VC_STAR.Matrix(),
          F(0), X11_STAR_STAR.Matrix() );
        MakeTriangular( UPPER, X11_STAR_STAR );
        A11.SumScatterUpdate( F(-1), X11_STAR_STAR );

        // A11 := inv(U11)' A11 inv(U11)
        A11_STAR_STAR = A11;
        U11_STAR_STAR = U11;
        LocalTwoSidedTrsm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A12 := A12 - U01' A02
        U01_MC_STAR = U01;
        LocalGemm( ADJOINT, NORMAL, F(1), U01_MC_STAR, A02, X12_STAR_MR );
        A12.ColSumScatterUpdate( F(-1), X12_STAR_MR );

        // A12 := inv(U11)' A12
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        A01_VC_STAR = A01;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // Y02 := Y02 + A01 U12
        A01_MC_STAR = A01;
        U12.AdjointColAllGather( U12Adj_MR_STAR );
        LocalGemm
        ( NORMAL, ADJOINT, F(1), A01_MC_STAR, U12Adj_MR_STAR, F(1), Y02 );

        // Y12 := Y12 + A11 U12
        MakeHermitian( UPPER, A11_STAR_STAR );
        A11_MC_STAR = A11_STAR_STAR;
        LocalGemm
        ( NORMAL, ADJOINT, F(1), A11_MC_STAR, U12Adj_MR_STAR, F(0), Y12 );

        // Y12 := Y12 + A01' U02
        LocalGemm( ADJOINT, NORMAL, F(1), A01_MC_STAR, U02, Z12_STAR_MR );
        Y12.ColSumScatterUpdate( F(1), Z12_STAR_MR );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, Y01, /**/ Y02,
               /**/       Y10, Y11, /**/ Y12,
         /*************/ /******************/
          YBL, /**/ YBR,  Y20, Y21, /**/ Y22 );

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /**********************************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TWOSIDEDTRSM_UVAR3_HPP
