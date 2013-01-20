/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TWOSIDEDTRSM_UVAR2_HPP
#define BLAS_TWOSIDEDTRSM_UVAR2_HPP

namespace elem {
namespace internal {

template<typename F> 
inline void
TwoSidedTrsmUVar2( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmUVar2");
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

        // A01 := A01 - 1/2 Y01
        Axpy( F(-1)/F(2), Y01, A01 );
        
        // A11 := A11 - (U01' A01 + A01' U01)
        Her2k( UPPER, ADJOINT, F(-1), U01, A01, F(1), A11 );

        // A11 := inv(U11)' A11 inv(U11)
        TwoSidedTrsmUUnb( diag, A11, U11 );

        // A12 := A12 - A02' U01
        Gemm( ADJOINT, NORMAL, F(-1), A02, U01, F(1), A12 );

        // A12 := inv(U11)' A12
        Trsm( LEFT, UPPER, ADJOINT, diag, F(1), U11, A12 );
        
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

// This routine has only partially been optimized. The ReduceScatter operations
// need to be (conjugate-)transposed in order to play nice with cache.
template<typename F> 
inline void
TwoSidedTrsmUVar2
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& U )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmUVar2");
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
    DistMatrix<F,MC,  STAR> A01_MC_STAR(g);
    DistMatrix<F,VC,  STAR> A01_VC_STAR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A12_STAR_VR(g);
    DistMatrix<F,MC,  STAR> F01_MC_STAR(g);
    DistMatrix<F,MC,  STAR> U01_MC_STAR(g);
    DistMatrix<F,VR,  STAR> U01_VR_STAR(g);
    DistMatrix<F,STAR,MR  > U01Adj_STAR_MR(g);
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > X11_STAR_MR(g);
    DistMatrix<F,MR,  STAR> X12Adj_MR_STAR(g);
    DistMatrix<F,MR,  MC  > X12Adj_MR_MC(g);
    DistMatrix<F,MR,  MC  > Y01_MR_MC(g);
    DistMatrix<F,MR,  STAR> Y01_MR_STAR(g);
    DistMatrix<F> X11(g);
    DistMatrix<F> Y01(g);

    Matrix<F> X12Local;

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

        A01_MC_STAR.AlignWith( U01 );
        Y01.AlignWith( A01 );
        Y01_MR_STAR.AlignWith( A00 );
        U01_MC_STAR.AlignWith( A00 );
        U01_VR_STAR.AlignWith( A00 );
        U01Adj_STAR_MR.AlignWith( A00 );
        X11_STAR_MR.AlignWith( U01 );
        X11.AlignWith( A11 );
        X12Adj_MR_STAR.AlignWith( A02 );
        X12Adj_MR_MC.AlignWith( A12 );
        F01_MC_STAR.AlignWith( A00 );
        //--------------------------------------------------------------------//
        // Y01 := A00 U01
        U01_MC_STAR = U01;
        U01_VR_STAR = U01_MC_STAR;
        U01Adj_STAR_MR.AdjointFrom( U01_VR_STAR );
        Y01_MR_STAR.ResizeTo( A01.Height(), A01.Width() ); 
        F01_MC_STAR.ResizeTo( A01.Height(), A01.Width() );
        Zero( Y01_MR_STAR );
        Zero( F01_MC_STAR );
        LocalSymmetricAccumulateLU
        ( ADJOINT, 
          F(1), A00, U01_MC_STAR, U01Adj_STAR_MR, F01_MC_STAR, Y01_MR_STAR );
        Y01_MR_MC.SumScatterFrom( Y01_MR_STAR );
        Y01 = Y01_MR_MC;
        Y01.SumScatterUpdate( F(1), F01_MC_STAR );

        // X11 := U01' A01
        X11_STAR_MR.ResizeTo( A11.Height(), A11.Width() );
        LocalGemm( ADJOINT, NORMAL, F(1), U01_MC_STAR, A01, F(0), X11_STAR_MR );

        // A01 := A01 - Y01
        Axpy( F(-1), Y01, A01 );
        A01_MC_STAR = A01;
        
        // A11 := A11 - triu(X11 + A01' U01) = A11 - (U01 A01 + A01' U01)
        LocalGemm( ADJOINT, NORMAL, F(1), A01_MC_STAR, U01, F(1), X11_STAR_MR );
        X11.SumScatterFrom( X11_STAR_MR );
        MakeTrapezoidal( LEFT, UPPER, 0, X11 );
        Axpy( F(-1), X11, A11 );

        // A01 := A01 inv(U11)
        U11_STAR_STAR = U11;
        A01_VC_STAR = A01_MC_STAR;
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, F(1), U11_STAR_STAR, A01_VC_STAR );
        A01 = A01_VC_STAR;

        // A11 := inv(U11)' A11 inv(U11)
        A11_STAR_STAR = A11;
        LocalTwoSidedTrsm( UPPER, diag, A11_STAR_STAR, U11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A12 := A12 - A02' U01
        X12Adj_MR_STAR.ResizeTo( A12.Width(), A12.Height() );
        LocalGemm
        ( ADJOINT, NORMAL, F(1), A02, U01_MC_STAR, F(0), X12Adj_MR_STAR );
        X12Adj_MR_MC.SumScatterFrom( X12Adj_MR_STAR );
        Adjoint( X12Adj_MR_MC.LockedLocalMatrix(), X12Local );
        Axpy( F(-1), X12Local, A12.LocalMatrix() );

        // A12 := inv(U11)' A12
        A12_STAR_VR = A12;
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, diag, F(1), U11_STAR_STAR, A12_STAR_VR );
        A12 = A12_STAR_VR;
        //--------------------------------------------------------------------//
        A01_MC_STAR.FreeAlignments();
        F01_MC_STAR.FreeAlignments();
        U01_MC_STAR.FreeAlignments();
        U01_VR_STAR.FreeAlignments();
        U01Adj_STAR_MR.FreeAlignments();
        X11.FreeAlignments();
        X11_STAR_MR.FreeAlignments();
        X12Adj_MR_STAR.FreeAlignments();
        X12Adj_MR_MC.FreeAlignments();
        Y01.FreeAlignments();
        Y01_MR_STAR.FreeAlignments();

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

#endif // ifndef BLAS_TWOSIDEDTRSM_UVAR2_HPP
