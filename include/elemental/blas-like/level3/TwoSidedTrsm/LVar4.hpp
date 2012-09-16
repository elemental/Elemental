/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {
namespace internal {

template<typename F> 
inline void
TwoSidedTrsmLVar4( UnitOrNonUnit diag, Matrix<F>& A, const Matrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmLVar4");
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
        // A10 := inv(L11) A10
        Trsm( LEFT, LOWER, NORMAL, diag, (F)1, L11, A10 );

        // A11 := inv(L11) A11 inv(L11)'
        TwoSidedTrsmLUnb( diag, A11, L11 );

        // A20 := A20 - L21 A10
        Gemm( NORMAL, NORMAL, (F)-1, L21, A10, (F)1, A20 );

        // Y21 := L21 A11
        Zeros( A21.Height(), A21.Width(), Y21 );
        Hemm( RIGHT, LOWER, (F)1, A11, L21, (F)0, Y21 );

        // A21 := A21 inv(L11)'
        Trsm( RIGHT, LOWER, ADJOINT, diag, (F)1, L11, A21 );

        // A21 := A21 - 1/2 Y21
        Axpy( (F)-0.5, Y21, A21 );

        // A22 := A22 - (L21 A21' + A21 L21')
        Her2k( LOWER, NORMAL, (F)-1, L21, A21, (F)1, A22 );

        // A21 := A21 - 1/2 Y21
        Axpy( (F)-0.5, Y21, A21 );
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

template<typename F> 
inline void
TwoSidedTrsmLVar4
( UnitOrNonUnit diag, DistMatrix<F>& A, const DistMatrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::TwoSidedTrsmLVar4");
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
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> A21_VR_STAR(g);
    DistMatrix<F,STAR,MC  > A21Trans_STAR_MC(g);
    DistMatrix<F,STAR,MR  > A21Adj_STAR_MR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> L21_VC_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > L21Adj_STAR_MR(g);
    DistMatrix<F,VC,  STAR> Y21_VC_STAR(g);

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

        A10_STAR_MR.AlignWith( A20 );
        A10_STAR_VR.AlignWith( A20 );
        A21_VC_STAR.AlignWith( A22 );
        A21_VR_STAR.AlignWith( A22 );
        A21Trans_STAR_MC.AlignWith( A22 );
        A21Adj_STAR_MR.AlignWith( A22 );
        L21_VC_STAR.AlignWith( A22 );
        L21_VR_STAR.AlignWith( A22 );
        L21_MC_STAR.AlignWith( A22 );
        L21Adj_STAR_MR.AlignWith( A22 );
        Y21_VC_STAR.AlignWith( A22 );
        //--------------------------------------------------------------------//
        // A10 := inv(L11) A10
        L11_STAR_STAR = L11;
        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, diag, (F)1, L11_STAR_STAR, A10_STAR_VR );

        // A11 := inv(L11) A11 inv(L11)'
        A11_STAR_STAR = A11; 
        LocalTwoSidedTrsm( LOWER, diag, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A20 := A20 - L21 A10
        L21_MC_STAR = L21;
        A10_STAR_MR = A10_STAR_VR;
        LocalGemm( NORMAL, NORMAL, (F)-1, L21_MC_STAR, A10_STAR_MR, (F)1, A20 );
        A10 = A10_STAR_MR; // delayed write from  A10 := inv(L11) A10

        // Y21 := L21 A11
        L21_VC_STAR = L21_MC_STAR;
        Y21_VC_STAR.ResizeTo( A21.Height(), A21.Width() );
        Hemm
        ( RIGHT, LOWER, 
          (F)1, A11_STAR_STAR.LocalMatrix(), L21_VC_STAR.LocalMatrix(), 
          (F)0, Y21_VC_STAR.LocalMatrix() );

        // A21 := A21 inv(L11)'
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, diag, (F)1, L11_STAR_STAR, A21_VC_STAR );

        // A21 := A21 - 1/2 Y21
        Axpy( (F)-0.5, Y21_VC_STAR, A21_VC_STAR );

        // A22 := A22 - (L21 A21' + A21 L21')
        A21Trans_STAR_MC.TransposeFrom( A21_VC_STAR );
        A21_VR_STAR = A21_VC_STAR;
        L21_VR_STAR = L21_VC_STAR;
        A21Adj_STAR_MR.AdjointFrom( A21_VR_STAR );
        L21Adj_STAR_MR.AdjointFrom( L21_VR_STAR );
        LocalTrr2k
        ( LOWER, TRANSPOSE,
          (F)-1, L21_MC_STAR,      A21Adj_STAR_MR, 
                 A21Trans_STAR_MC, L21Adj_STAR_MR,
          (F)1,  A22 );

        // A21 := A21 - 1/2 Y21
        Axpy( (F)-0.5, Y21_VC_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A10_STAR_MR.FreeAlignments();
        A10_STAR_VR.FreeAlignments();
        A21_VC_STAR.FreeAlignments();
        A21_VR_STAR.FreeAlignments();
        A21Trans_STAR_MC.FreeAlignments();
        A21Adj_STAR_MR.FreeAlignments();
        L21_VC_STAR.FreeAlignments();
        L21_VR_STAR.FreeAlignments();
        L21_MC_STAR.FreeAlignments();
        L21Adj_STAR_MR.FreeAlignments();
        Y21_VC_STAR.FreeAlignments();

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
