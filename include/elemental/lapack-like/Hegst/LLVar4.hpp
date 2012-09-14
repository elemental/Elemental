/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   Copyright (c) 2012, The University of Texas at Austin
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
HegstLLVar4( DistMatrix<F>& A, const DistMatrix<F>& L )
{
#ifndef RELEASE
    PushCallStack("internal::HegstLLVar4");
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
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,MR  > A10_STAR_MR(g);
    DistMatrix<F,STAR,MC  > A10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,VC,  STAR> A21_VC_STAR(g);
    DistMatrix<F,MC,  STAR> A21_MC_STAR(g);
    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,MR,  STAR> L10Adj_MR_STAR(g);
    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > Y10_STAR_VR(g);

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

        A10_STAR_VR.AlignWith( A00 );
        A10_STAR_MR.AlignWith( A00 );
        A10_STAR_MC.AlignWith( A00 );
        A21_MC_STAR.AlignWith( A20 );
        L10_STAR_VR.AlignWith( A00 );
        L10Adj_MR_STAR.AlignWith( A00 );
        L10_STAR_MC.AlignWith( A00 );
        Y10_STAR_VR.AlignWith( A10 );
        //--------------------------------------------------------------------//
        // Y10 := A11 L10
        A11_STAR_STAR = A11;
        L10Adj_MR_STAR.AdjointFrom( L10 );
        L10_STAR_VR.AdjointFrom( L10Adj_MR_STAR );
        Y10_STAR_VR.ResizeTo( A10.Height(), A10.Width() );
        Zero( Y10_STAR_VR );
        Hemm
        ( LEFT, LOWER,
          (F)0.5, A11_STAR_STAR.LockedLocalMatrix(),
                  L10_STAR_VR.LockedLocalMatrix(),
          (F)0,   Y10_STAR_VR.LocalMatrix() );

        // A10 := A10 + 1/2 Y10
        A10_STAR_VR = A10;
        Axpy( (F)1, Y10_STAR_VR, A10_STAR_VR );

        // A00 := A00 + (A10' L10 + L10' A10)
        A10_STAR_MR = A10_STAR_VR;
        A10_STAR_MC = A10_STAR_VR;
        L10_STAR_MC = L10_STAR_VR;
        LocalTrr2k
        ( LOWER, ADJOINT, ADJOINT, ADJOINT,
          (F)1, A10_STAR_MC, L10Adj_MR_STAR, 
                L10_STAR_MC, A10_STAR_MR, 
          (F)1, A00 );

        // A10 := A10 + 1/2 Y10
        Axpy( (F)1, Y10_STAR_VR, A10_STAR_VR );

        // A10 := L11' A10
        L11_STAR_STAR = L11;
        LocalTrmm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, (F)1, L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;

        // A20 := A20 + A21 L10
        A21_MC_STAR = A21;
        LocalGemm
        ( NORMAL, ADJOINT, (F)1, A21_MC_STAR, L10Adj_MR_STAR, (F)1, A20 );

        // A11 := L11' A11 L11
        LocalHegst( LEFT, LOWER, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A21 := A21 L11
        A21_VC_STAR = A21_MC_STAR;
        LocalTrmm
        ( RIGHT, LOWER, NORMAL, NON_UNIT, (F)1, L11_STAR_STAR, A21_VC_STAR );
        A21 = A21_VC_STAR;
        //--------------------------------------------------------------------//
        A10_STAR_VR.FreeAlignments();
        A10_STAR_MR.FreeAlignments();
        A10_STAR_MC.FreeAlignments();
        A21_MC_STAR.FreeAlignments();
        L10_STAR_VR.FreeAlignments();
        L10Adj_MR_STAR.FreeAlignments();
        L10_STAR_MC.FreeAlignments();
        Y10_STAR_VR.FreeAlignments();

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

} // namespace internal
} // namespace elem
