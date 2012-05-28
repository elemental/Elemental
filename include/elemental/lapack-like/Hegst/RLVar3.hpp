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
HegstRLVar3( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("internal::HegstRLVar4");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( L.Height() != L.Width() )
        throw std::logic_error("Triangular matrices must be square");
    if( A.Height() != L.Height() )
        throw std::logic_error("A and L must be the same size");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<F,MC,MR>
        YTL(g), YTR(g),  Y00(g), Y01(g), Y02(g),
        YBL(g), YBR(g),  Y10(g), Y11(g), Y12(g),
                         Y20(g), Y21(g), Y22(g);
    DistMatrix<F,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<F,STAR,MR  > A11_STAR_MR(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,STAR,MR  > L10_STAR_MR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,STAR> X11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> X21_MC_STAR(g);
    DistMatrix<F,MC,  STAR> Z21_MC_STAR(g);

    // We will use an entire extra matrix as temporary storage. If this is not
    // acceptable, use HegstRLVar4 instead.
    DistMatrix<F,MC,MR> Y(g);
    Y.AlignWith( A );
    Y.ResizeTo( A.Height(), A.Width() );
    Zero( Y );

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDownDiagonal
    ( Y, YTL, YTR,
         YBL, YBR, 0 );
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

        RepartitionDownDiagonal
        ( YTL, /**/ YTR,  Y00, /**/ Y01, Y02,
         /*************/ /******************/
               /**/       Y10, /**/ Y11, Y12,
          YBL, /**/ YBR,  Y20, /**/ Y21, Y22 );

        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        A11_STAR_MR.AlignWith( Y21 );
        A21_VC_STAR.AlignWith( A21 );
        A10_STAR_VR.AlignWith( A10 );
        A10_STAR_MR.AlignWith( A10 );
        L10_STAR_VR.AlignWith( A10 );
        L10_STAR_MR.AlignWith( A10 );
        L21_MC_STAR.AlignWith( Y21 );
        X21_MC_STAR.AlignWith( A20 );
        Z21_MC_STAR.AlignWith( L20 );
        //--------------------------------------------------------------------//
        // A10 := A10 - 1/2 Y10
        Axpy( (F)-0.5, Y10, A10 );

        // A11 := A11 - (A10 L10' + L10 A10')
        A10_STAR_VR = A10;
        L10_STAR_VR = L10;
        X11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        Her2k
        ( LOWER, NORMAL, 
          (F)1, A10_STAR_VR.LocalMatrix(), L10_STAR_VR.LocalMatrix(),
          (F)0, X11_STAR_STAR.LocalMatrix() );
        MakeTrapezoidal( LEFT, LOWER, 0, X11_STAR_STAR );
        A11.SumScatterUpdate( (F)-1, X11_STAR_STAR );

        // A11 := inv(L11) A11 inv(L11)'
        A11_STAR_STAR = A11;
        L11_STAR_STAR = L11;
        LocalHegst( RIGHT, LOWER, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A21 := A21 - A20 L10'
        L10_STAR_MR = L10_STAR_VR;
        X21_MC_STAR.ResizeTo( A21.Height(), A21.Width() );
        LocalGemm
        ( NORMAL, ADJOINT, (F)1, A20, L10_STAR_MR, (F)0, X21_MC_STAR );
        A21.SumScatterUpdate( (F)-1, X21_MC_STAR );

        // A21 := A21 inv(L11)'
        A21_VC_STAR = A21;
        LocalTrsm
        ( RIGHT, LOWER, ADJOINT, NON_UNIT, (F)1, L11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;

        // A10 := A10 - 1/2 Y10
        Axpy( (F)-0.5, Y10, A10 );

        // A10 := inv(L11) A10
        A10_STAR_VR = A10;
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT,
          (F)1, L11_STAR_STAR, A10_STAR_VR );

        // Y20 := Y20 + L21 A10
        A10_STAR_MR = A10_STAR_VR;
        A10 = A10_STAR_MR;
        L21_MC_STAR = L21;
        LocalGemm
        ( NORMAL, NORMAL, (F)1, L21_MC_STAR, A10_STAR_MR, (F)1, Y20 );

        // Y21 := L21 A11
        //
        // Symmetrize A11[* ,* ] by copying the lower triangle into the upper
        // so that we can call a local gemm instead of worrying about
        // reproducing a hemm with nonsymmetric local matrices.
        {
            const int height = A11_STAR_STAR.LocalHeight();
            const int ldim = A11_STAR_STAR.LocalLDim();
            F* A11Buffer = A11_STAR_STAR.LocalBuffer();
            for( int i=1; i<height; ++i )
                for( int j=0; j<i; ++j )
                    A11Buffer[j+i*ldim] = Conj(A11Buffer[i+j*ldim]);
        }
        A11_STAR_MR = A11_STAR_STAR;
        LocalGemm
        ( NORMAL, NORMAL, (F)1, L21_MC_STAR, A11_STAR_MR, (F)0, Y21 );

        // Y21 := Y21 + L20 A10'
        Z21_MC_STAR.ResizeTo( A21.Height(), A21.Width() );
        LocalGemm
        ( NORMAL, ADJOINT, (F)1, L20, A10_STAR_MR, (F)0, Z21_MC_STAR );
        Y21.SumScatterUpdate( (F)1, Z21_MC_STAR );
        //--------------------------------------------------------------------//
        A11_STAR_MR.FreeAlignments();
        A21_VC_STAR.FreeAlignments();
        A10_STAR_VR.FreeAlignments();
        A10_STAR_MR.FreeAlignments();
        L10_STAR_VR.FreeAlignments();
        L10_STAR_MR.FreeAlignments();
        L21_MC_STAR.FreeAlignments();
        X21_MC_STAR.FreeAlignments();
        Z21_MC_STAR.FreeAlignments();

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
