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

template<typename F> 
inline void
internal::HegstRLVar1( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("internal::HegstRLVar1");
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
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g);
    DistMatrix<F,STAR,VR  > A10_STAR_VR(g);
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > L10_STAR_MC(g);
    DistMatrix<F,STAR,VR  > L10_STAR_VR(g);
    DistMatrix<F,MR,  STAR> L10Adj_MR_STAR(g);
    DistMatrix<F,VC,  STAR> L10Adj_VC_STAR(g);
    DistMatrix<F,STAR,STAR> X11_STAR_STAR(g);
    DistMatrix<F,MC,  MR  > Y10(g);
    DistMatrix<F,MC,  MR  > Z10Adj(g);
    DistMatrix<F,MR,  MC  > Z10Adj_MR_MC(g);
    DistMatrix<F,MC,  STAR> Z10Adj_MC_STAR(g);
    DistMatrix<F,MR,  STAR> Z10Adj_MR_STAR(g);

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

        A10_STAR_VR.AlignWith( A10 );
        L10_STAR_MC.AlignWith( A00 );
        L10_STAR_VR.AlignWith( A00 );
        L10Adj_MR_STAR.AlignWith( A00 );
        L10Adj_VC_STAR.AlignWith( A00 );
        Y10.AlignWith( A10 );
        Z10Adj.AlignWith( A10 );
        Z10Adj_MR_MC.AlignWith( A10 );
        Z10Adj_MC_STAR.AlignWith( A00 );
        Z10Adj_MR_STAR.AlignWith( A00 );
        //--------------------------------------------------------------------//
        // Y10 := L10 A00
        L10Adj_MR_STAR.AdjointFrom( L10 );
        L10Adj_VC_STAR = L10Adj_MR_STAR;
        L10_STAR_MC.AdjointFrom( L10Adj_VC_STAR );
        Z10Adj_MC_STAR.ResizeTo( A10.Width(), A10.Height() );
        Z10Adj_MR_STAR.ResizeTo( A10.Width(), A10.Height() );
        Z10Adj_MC_STAR.SetToZero();
        Z10Adj_MR_STAR.SetToZero();
        internal::LocalSymmetricAccumulateRL
        ( ADJOINT,
          (F)1, A00, L10_STAR_MC, L10Adj_MR_STAR, 
          Z10Adj_MC_STAR, Z10Adj_MR_STAR );
        Z10Adj.SumScatterFrom( Z10Adj_MC_STAR );
        Z10Adj_MR_MC = Z10Adj;
        Z10Adj_MR_MC.SumScatterUpdate( (F)1, Z10Adj_MR_STAR );
        Y10.ResizeTo( A10.Height(), A10.Width() );
        Adjoint( Z10Adj_MR_MC.LocalMatrix(), Y10.LocalMatrix() );

        // A10 := A10 inv(L00)'
        // This is the bottleneck because A10 only has blocksize rows
        Trsm( RIGHT, LOWER, ADJOINT, NON_UNIT, (F)1, L00, A10 );

        // A10 := A10 - 1/2 Y10
        Axpy( (F)-0.5, Y10, A10 );

        // A11 := A11 - (A10 L10' + L10 A10')
        A10_STAR_VR = A10;
        L10_STAR_VR = L10;
        X11_STAR_STAR.ResizeTo( A11.Height(), A11.Width() );
        Her2k
        ( LOWER, NORMAL,
          (F)-1, A10_STAR_VR.LocalMatrix(), L10_STAR_VR.LocalMatrix(), 
          (F)0, X11_STAR_STAR.LocalMatrix() );
        A11.SumScatterUpdate( (F)1, X11_STAR_STAR );

        // A11 := inv(L11) A11 inv(L11)'
        A11_STAR_STAR = A11;
        L11_STAR_STAR = L11;
        internal::LocalHegst( RIGHT, LOWER, A11_STAR_STAR, L11_STAR_STAR );
        A11 = A11_STAR_STAR;

        // A10 := A10 - 1/2 Y10
        Axpy( (F)-0.5, Y10, A10 );

        // A10 := inv(L11) A10
        A10_STAR_VR = A10;
        internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, (F)1, L11_STAR_STAR, A10_STAR_VR );
        A10 = A10_STAR_VR;
        //--------------------------------------------------------------------//
        A10_STAR_VR.FreeAlignments();
        L10_STAR_MC.FreeAlignments();
        L10_STAR_VR.FreeAlignments();
        L10Adj_MR_STAR.FreeAlignments();
        L10Adj_VC_STAR.FreeAlignments();
        Y10.FreeAlignments();
        Z10Adj.FreeAlignments();
        Z10Adj_MR_MC.FreeAlignments();
        Z10Adj_MC_STAR.FreeAlignments();
        Z10Adj_MR_STAR.FreeAlignments();

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

} // namespace elem
