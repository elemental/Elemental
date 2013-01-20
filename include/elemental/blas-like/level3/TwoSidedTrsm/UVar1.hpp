/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TWOSIDEDTRSM_UVAR1_HPP
#define BLAS_TWOSIDEDTRSM_UVAR1_HPP

namespace elem {
namespace internal {

template<typename F> 
inline void
TwoSidedTrsmUVar1( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmUVar1");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( U.Height() != U.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != U.Height() )
        throw std::logic_error("A and U must be the same size");
#endif
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
    Matrix<F> Y01;

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
        // Y01 := A00 U01
        Zeros( A01.Height(), A01.Width(), Y01 );
        Hemm( LEFT, UPPER, F(1), A00, U01, F(0), Y01 );

        // A01 := inv(U00)' A01
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U00, A01 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A11 := A11 - (U01' A01 + A01' U01)
        Her2k( UPPER, ADJOINT, F(-1), U01, A01, F(1), A11 );

        // A11 := inv(U11)' A11 inv(U11)
        TwoSidedTrsmUUnb( diag, A11, U11 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        Trsm( RIGHT, UPPER, NORMAL, diag, F(1), U11, A01 );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
TwoSidedTrsmUVar1
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmUVar1");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( U.Height() != U.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != U.Height() )
        throw std::logic_error("A and U must be the same size");
#endif
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
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,VC,  STAR> U01_VC_STAR(g);
    DistMatrix<F,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<F,STAR,MR  > U01Adj_STAR_MR(g);
    DistMatrix<F,STAR,STAR> X11_STAR_STAR(g);
    DistMatrix<F,MR,  MC  > Z01_MR_MC(g);
    DistMatrix<F,MC,  STAR> Z01_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Z01_MR_STAR(g);
    DistMatrix<F> Y01(g);

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

        A01_VC_STAR.AlignWith( A01 );
        U01_MC_STAR.AlignWith( A00 );
        U01_VR_STAR.AlignWith( A00 );
        U01_VC_STAR.AlignWith( A00 );
        U01Adj_STAR_MR.AlignWith( A00 );
        Y01.AlignWith( A01 );
        Z01_MR_MC.AlignWith( A01 );
        Z01_MC_STAR.AlignWith( A00 );
        Z01_MR_STAR.AlignWith( A00 );
        //--------------------------------------------------------------------//
        // Y01 := A00 U01
        U01_MC_STAR = U01;
        U01_VR_STAR = U01_MC_STAR;
        U01Adj_STAR_MR.AdjointFrom( U01_VR_STAR );
        Z01_MC_STAR.ResizeTo( A01.Height(), A01.Width() );
        Z01_MR_STAR.ResizeTo( A01.Height(), A01.Width() );
        Zero( Z01_MC_STAR );
        Zero( Z01_MR_STAR );
        LocalSymmetricAccumulateLU
        ( ADJOINT, 
          F(1), A00, U01_MC_STAR, U01Adj_STAR_MR, Z01_MC_STAR, Z01_MR_STAR );
        Z01_MR_MC.SumScatterFrom( Z01_MR_STAR );
        Y01 = Z01_MR_MC;
        Y01.SumScatterUpdate( F(1), Z01_MC_STAR );

        // A01 := inv(U00)' A01
        //
        // This is the bottleneck because A01 only has blocksize columns
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U00, A01 );

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A11 := A11 - (U01' A01 + A01' U01)
        A01_VC_STAR = A01;
        U01_VC_STAR = U01_MC_STAR;
        X11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        Her2k
        ( UPPER, ADJOINT,
          F(-1), A01_VC_STAR.LocalMatrix(), U01_VC_STAR.LocalMatrix(),
          F(0), X11_STAR_STAR.LocalMatrix() );
        A11.SumScatterUpdate( F(1), X11_STAR_STAR );

        // A11 := inv(U11)' A11 inv(U11)
        A11_STAR_STAR = A11;
        U11_STAR_STAR = U11;
        LocalTwoSidedTrsm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );

        // A01 := A01 inv(U11)
        A01_VC_STAR = A01;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;
        //--------------------------------------------------------------------//
        A01_VC_STAR.FreeAlignments();
        U01_MC_STAR.FreeAlignments();
        U01_VR_STAR.FreeAlignments();
        U01_VC_STAR.FreeAlignments();
        U01Adj_STAR_MR.FreeAlignments();
        Y01.FreeAlignments();
        Z01_MR_MC.FreeAlignments();
        Z01_MC_STAR.FreeAlignments();
        Z01_MR_STAR.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TWOSIDEDTRSM_UVAR1_HPP
