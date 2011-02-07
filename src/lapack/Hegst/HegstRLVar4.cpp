/*
   Copyright (c) 2009-2011, Jack Poulson
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
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

template<typename F> // F represents a real or complex field
void
elemental::lapack::internal::HegstRLVar4
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstRLVar4");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( L.Height() != L.Width() )
        throw logic_error( "Triangular matrices must be square." );
    if( A.Height() != L.Height() )
        throw logic_error( "A and L must be the same size." );
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
    DistMatrix<F,Star,MR  > A10_Star_MR(g);
    DistMatrix<F,Star,VR  > A10_Star_VR(g);
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,VC,  Star> A21_VC_Star(g);
    DistMatrix<F,VR,  Star> A21_VR_Star(g);
    DistMatrix<F,Star,MC  > A21Trans_Star_MC(g);
    DistMatrix<F,Star,MR  > A21Herm_Star_MR(g);
    DistMatrix<F,Star,Star> L11_Star_Star(g);
    DistMatrix<F,VC,  Star> L21_VC_Star(g);
    DistMatrix<F,VR,  Star> L21_VR_Star(g);
    DistMatrix<F,MC,  Star> L21_MC_Star(g);
    DistMatrix<F,Star,MR  > L21Herm_Star_MR(g);
    DistMatrix<F,VC,  Star> Y21_VC_Star(g);

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

        A10_Star_MR.AlignWith( A20 );
        A10_Star_VR.AlignWith( A20 );
        A21_VC_Star.AlignWith( A22 );
        A21_VR_Star.AlignWith( A22 );
        A21Trans_Star_MC.AlignWith( A22 );
        A21Herm_Star_MR.AlignWith( A22 );
        L21_VC_Star.AlignWith( A22 );
        L21_VR_Star.AlignWith( A22 );
        L21_MC_Star.AlignWith( A22 );
        L21Herm_Star_MR.AlignWith( A22 );
        Y21_VC_Star.AlignWith( A22 );
        //--------------------------------------------------------------------//
        // A10 := inv(L11) A10
        L11_Star_Star = L11;
        A10_Star_VR = A10;
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, NonUnit,
          (F)1, L11_Star_Star, A10_Star_VR );
        A10 = A10_Star_VR;

        // A11 := inv(L11) A11 inv(L11)'
        A11_Star_Star = A11; 
        lapack::internal::LocalHegst
        ( Right, Lower, A11_Star_Star, L11_Star_Star );
        A11 = A11_Star_Star;

        // A20 := A20 - L21 A10
        L21_MC_Star = L21;
        A10_Star_MR = A10_Star_VR;
        blas::internal::LocalGemm
        ( Normal, Normal, (F)-1, L21_MC_Star, A10_Star_MR, (F)1, A20 );

        // Y21 := L21 A11
        L21_VC_Star = L21_MC_Star;
        Y21_VC_Star.ResizeTo( A21.Height(), A21.Width() );
        blas::Hemm
        ( Right, Lower, 
          (F)1, A11_Star_Star.LocalMatrix(), L21_VC_Star.LocalMatrix(), 
          (F)0, Y21_VC_Star.LocalMatrix() );

        // A21 := A21 inv(L11)'
        A21_VC_Star = A21;
        blas::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit,
          (F)1, L11_Star_Star, A21_VC_Star );

        // A21 := A21 - 1/2 Y21
        blas::Axpy( (F)-0.5, Y21_VC_Star, A21_VC_Star );

        // A22 := A22 - (L21 A21' + A21 L21')
        A21Trans_Star_MC.TransposeFrom( A21_VC_Star );
        A21_VR_Star = A21_VC_Star;
        L21_VR_Star = L21_VC_Star;
        A21Herm_Star_MR.ConjugateTransposeFrom( A21_VR_Star );
        L21Herm_Star_MR.ConjugateTransposeFrom( L21_VR_Star );
        blas::internal::LocalTriangularRank2K
        ( Lower, Transpose,
          (F)-1, L21_MC_Star, A21Trans_Star_MC, 
                 L21Herm_Star_MR, A21Herm_Star_MR,
          (F)1, A22 );

        // A21 := A21 - 1/2 Y21
        blas::Axpy( (F)-0.5, Y21_VC_Star, A21_VC_Star );
        A21 = A21_VC_Star;
        //--------------------------------------------------------------------//
        A10_Star_MR.FreeAlignments();
        A10_Star_VR.FreeAlignments();
        A21_VC_Star.FreeAlignments();
        A21_VR_Star.FreeAlignments();
        A21Trans_Star_MC.FreeAlignments();
        A21Herm_Star_MR.FreeAlignments();
        L21_VC_Star.FreeAlignments();
        L21_VR_Star.FreeAlignments();
        L21_MC_Star.FreeAlignments();
        L21Herm_Star_MR.FreeAlignments();
        Y21_VC_Star.FreeAlignments();

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

template void elemental::lapack::internal::HegstRLVar4
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::HegstRLVar4
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& L );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstRLVar4
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::HegstRLVar4
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& L );
#endif

