/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TWOSIDEDTRSM_LVAR2_HPP
#define BLAS_TWOSIDEDTRSM_LVAR2_HPP

namespace elem {
namespace internal {

template<typename F>
inline void
TwoSidedTrsmLVar2( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmLVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( L.Height() != L.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != L.Height() )
        throw std::logic_error("A and L must be the same size");
#endif
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
    Matrix<F> X11;
    Matrix<F> Y10;

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
        // Y10 := L10 A00
        Zeros( L10.Height(), A00.Width(), Y10 );
        Hemm( RIGHT, LOWER, F(1), A00, L10, F(0), Y10 );

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A11 := A11 - (A10 L10' + L10 A10')
        Her2k( LOWER, NORMAL, F(-1), A10, L10, F(1), A11 );

        // A11 := inv(L11) A11 inv(L11)'
        TwoSidedTrsmLUnb( diag, A11, L11 );

        // A21 := A21 - A20 L10'
        Gemm( NORMAL, ADJOINT, F(-1), A20, L10, F(1), A21 );

        // A21 := A21 inv(L11)'
        Trsm( RIGHT, LOWER, ADJOINT, diag, F(1), L11, A21 );

        // A10 := A10 - 1/2 Y10
        Axpy( F(-1)/F(2), Y10, A10 );

        // A10 := inv(L11) A10
        Trsm( LEFT, LOWER, NORMAL, diag, F(1), L11, A10 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /**********************************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// This routine has only partially been optimized. The ReduceScatter operations
// need to be (conjugate-)transposed in order to play nice with cache.
template<typename F>
inline void
TwoSidedTrsmLVar2
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmLVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( L.Height() != L.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != L.Height() )
        throw std::logic_error("A and L must be the same size");
#endif
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
    DistMatrix<F,MR,  STAR> A10Adj_MR_STAR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,MR,  STAR> F10Adj_MR_STAR(g);
    DistMatrix<F,MR,  STAR> L10Adj_MR_STAR(g);
    DistMatrix<F,VC,  STAR> L10Adj_VC_STAR(g);
    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> X11_MC_STAR(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);
    DistMatrix<F,MC,  STAR> Y10Adj_MC_STAR(g);
    DistMatrix<F,MR,  MC  > Y10Adj_MR_MC(g);
    DistMatrix<F> X11(g);
    DistMatrix<F> Y10Adj(g);

    Matrix<F> Y10Local;

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

        A10Adj_MR_STAR.AlignWith( L10 );
        F10Adj_MR_STAR.AlignWith( A00 );
        L10Adj_MR_STAR.AlignWith( A00 );
        L10Adj_VC_STAR.AlignWith( A00 );
        L10_STAR_MC.AlignWith( A00 );
        X11.AlignWith( A11 );
        X11_MC_STAR.AlignWith( L10 );
        X21_MC_STAR.AlignWith( A20 );
        Y10Adj_MC_STAR.AlignWith( A00 );
        Y10Adj_MR_MC.AlignWith( A10 );
        //--------------------------------------------------------------------//
        // Y10 := L10 A00
        L10Adj_MR_STAR.AdjointFrom( L10 );
        L10Adj_VC_STAR = L10Adj_MR_STAR;
        L10_STAR_MC.AdjointFrom( L10Adj_VC_STAR );
        Y10Adj_MC_STAR.ResizeTo( A10.Width(), A10.Height() );
        F10Adj_MR_STAR.ResizeTo( A10.Width(), A10.Height() );
        Zero( Y10Adj_MC_STAR );
        Zero( F10Adj_MR_STAR );
        LocalSymmetricAccumulateRL
        ( ADJOINT,
          F(1), A00, L10_STAR_MC, L10Adj_MR_STAR, 
          Y10Adj_MC_STAR, F10Adj_MR_STAR );
        Y10Adj.SumScatterFrom( Y10Adj_MC_STAR );
        Y10Adj_MR_MC = Y10Adj;
        Y10Adj_MR_MC.SumScatterUpdate( F(1), F10Adj_MR_STAR );
        Adjoint( Y10Adj_MR_MC.LockedLocalMatrix(), Y10Local );

        // X11 := A10 L10'
        X11_MC_STAR.ResizeTo( A11.Height(), A11.Width() );
        LocalGemm
        ( NORMAL, NORMAL, F(1), A10, L10Adj_MR_STAR, F(0), X11_MC_STAR );

        // A10 := A10 - Y10
        Axpy( F(-1), Y10Local, A10.LocalMatrix() );
        A10Adj_MR_STAR.AdjointFrom( A10 );
        
        // A11 := A11 - (X11 + L10 A10') = A11 - (A10 L10' + L10 A10')
        LocalGemm
        ( NORMAL, NORMAL, F(1), L10, A10Adj_MR_STAR, F(1), X11_MC_STAR );
        X11.SumScatterFrom( X11_MC_STAR );
        MakeTrapezoidal( LEFT, LOWER, 0, X11 );
        Axpy( F(-1), X11, A11 );

        // A10 := inv(L11) A10
        L11_STAR_STAR = L11;
        A10_STAR_VR.AdjointFrom( A10Adj_MR_STAR );
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, F(1), L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;

        // A11 := inv(L11) A11 inv(L11)'
        A11_STAR_STAR = A11;
        LocalTwoSidedTrsm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A21 := A21 - A20 L10'
        X21_MC_STAR.ResizeTo( A21.Height(), A21.Width() );
        LocalGemm
        ( NORMAL, NORMAL, F(1), A20, L10Adj_MR_STAR, F(0), X21_MC_STAR );
        A21.SumScatterUpdate( F(-1), X21_MC_STAR );

        // A21 := A21 inv(L11)'
        A21_VC_STAR =  A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, diag, F(1), L11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A10Adj_MR_STAR.FreeAlignments();
        F10Adj_MR_STAR.FreeAlignments();
        L10Adj_MR_STAR.FreeAlignments();
        L10Adj_VC_STAR.FreeAlignments();
        L10_STAR_MC.FreeAlignments();
        X11.FreeAlignments();
        X11_MC_STAR.FreeAlignments();
        X21_MC_STAR.FreeAlignments();
        Y10Adj_MC_STAR.FreeAlignments();
        Y10Adj_MR_MC.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /**********************************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TWOSIDEDTRSM_LVAR2_HPP
