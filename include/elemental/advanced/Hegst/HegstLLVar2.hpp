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

namespace elemental {

template<typename F> 
inline void
internal::HegstLLVar2( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("internal::HegstLLVar2");
    if( A.Height() != A.Width() )
        throw std::logic_error( "A must be square." );
    if( L.Height() != L.Width() )
        throw std::logic_error( "Triangular matrices must be square." );
    if( A.Height() != L.Height() )
        throw std::logic_error( "A and L must be the same size." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<F,MC,MR>
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
    DistMatrix<F,MC,  MR  > Y21(g);
    DistMatrix<F,MC,  STAR> Z21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Z21_MR_STAR(g);
    DistMatrix<F,MR,  MC  > Z21_MR_MC(g);

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
        internal::LocalTrmm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, (F)1, L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;

        // A10 := A10 + L21' A20
        L21_MC_STAR = L21;
        X10_STAR_MR.ResizeTo( A10.Height(), A10.Width() );
        internal::LocalGemm
        ( ADJOINT, NORMAL, (F)1, L21_MC_STAR, A20, (F)0, X10_STAR_MR );
        A10.SumScatterUpdate( (F)1, X10_STAR_MR );

        // Y21 := A22 L21
        L21_VC_STAR = L21_MC_STAR;
        L21_VR_STAR = L21_VC_STAR;
        L21Adj_STAR_MR.AdjointFrom( L21_VR_STAR );
        Z21_MC_STAR.ResizeTo( A21.Height(), A21.Width() );
        Z21_MR_STAR.ResizeTo( A21.Height(), A21.Width() );
        Z21_MC_STAR.SetToZero();
        Z21_MR_STAR.SetToZero();
        internal::LocalSymmetricAccumulateLL
        ( ADJOINT, 
          (F)1, A22, L21_MC_STAR, L21Adj_STAR_MR, Z21_MC_STAR, Z21_MR_STAR );
        Z21_MR_MC.SumScatterFrom( Z21_MR_STAR );
        Y21 = Z21_MR_MC;
        Y21.SumScatterUpdate( (F)1, Z21_MC_STAR ); 

        // A21 := A21 L11
        A21_VC_STAR = A21;
        internal::LocalTrmm
        ( RIGHT, LOWER, NORMAL, NON_UNIT,
          (F)1, L11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;

        // A21 := A21 + 1/2 Y21
        Axpy( (F)0.5, Y21, A21 );

        // A11 := L11' A11 L11
        A11_STAR_STAR = A11;
        internal::LocalHegst( LEFT, LOWER, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A11 := A11 + (A21' L21 + L21' A21)
        A21_VC_STAR = A21;
        X11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        Her2k
        ( LOWER, ADJOINT,
          (F)1, A21_VC_STAR.LocalMatrix(), L21_VC_STAR.LocalMatrix(),
          (F)0, X11_STAR_STAR.LocalMatrix() );
        A11.SumScatterUpdate( (F)1, X11_STAR_STAR );

        // A21 := A21 + 1/2 Y21
        Axpy( (F)0.5, Y21, A21 );
        //--------------------------------------------------------------------//
        A21_VC_STAR.FreeAlignments();
        L21_MC_STAR.FreeAlignments();
        L21_VC_STAR.FreeAlignments();
        L21_VR_STAR.FreeAlignments();
        L21Adj_STAR_MR.FreeAlignments();
        X10_STAR_MR.FreeAlignments();
        Y21.FreeAlignments();
        Z21_MC_STAR.FreeAlignments();
        Z21_MR_STAR.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
