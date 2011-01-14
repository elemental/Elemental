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
elemental::lapack::internal::HegstRLVar3
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
        YTL(g), YTR(g),  Y00(g), Y01(g), Y02(g),
        YBL(g), YBR(g),  Y10(g), Y11(g), Y12(g),
                         Y20(g), Y21(g), Y22(g);
    DistMatrix<F,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<F,MC,  MR  > Y(g);
    DistMatrix<F,Star,MR  > A11_Star_MR(g);
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,VC,  Star> A21_VC_Star(g);
    DistMatrix<F,Star,VR  > A10_Star_VR(g);
    DistMatrix<F,Star,MR  > A10_Star_MR(g);
    DistMatrix<F,Star,Star> L11_Star_Star(g);
    DistMatrix<F,Star,VR  > L10_Star_VR(g);
    DistMatrix<F,Star,MR  > L10_Star_MR(g);
    DistMatrix<F,MC,  Star> L21_MC_Star(g);
    DistMatrix<F,Star,Star> X11_Star_Star(g);
    DistMatrix<F,MC,  Star> X21_MC_Star(g);
    DistMatrix<F,MC,  Star> Z21_MC_Star(g);

    // We will use an entire extra matrix as temporary storage. If this is not
    // acceptable, use HegstRLVar4 instead.
    Y.AlignWith( A );
    Y.ResizeTo( A.Height(), A.Width() );
    Y.SetToZero();

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

        A11_Star_MR.AlignWith( Y21 );
        A21_VC_Star.AlignWith( A21 );
	A10_Star_VR.AlignWith( A10 );
	A10_Star_MR.AlignWith( A10 );
	L10_Star_VR.AlignWith( A10 );
	L10_Star_MR.AlignWith( A10 );
        L21_MC_Star.AlignWith( Y21 );
        X21_MC_Star.AlignWith( A20 );
        Z21_MC_Star.AlignWith( L20 );
        X11_Star_Star.ResizeTo( A11.Height(), A11.Width() );
        X21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        Z21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        //--------------------------------------------------------------------//
        blas::Axpy( (F)-0.5, Y10, A10 );

	A10_Star_VR = A10;
	L10_Star_VR = L10;
        blas::Her2k
        ( Lower, Normal, 
	  (F)1, A10_Star_VR.LocalMatrix(), L10_Star_VR.LocalMatrix(),
          (F)0, X11_Star_Star.LocalMatrix() );
        X11_Star_Star.MakeTrapezoidal( Left, Lower );
        A11.SumScatterUpdate( (F)-1, X11_Star_Star );

        A11_Star_Star = A11;
        L11_Star_Star = L11;
        lapack::internal::LocalHegst
        ( Right, Lower, A11_Star_Star, L11_Star_Star );
        A11 = A11_Star_Star;

	L10_Star_MR = L10_Star_VR;
        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
	  (F)1, A20, L10_Star_MR, (F)0, X21_MC_Star );
        A21.SumScatterUpdate( (F)-1, X21_MC_Star );

        A21_VC_Star = A21;
        blas::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit,
          (F)1, L11_Star_Star, A21_VC_Star );
        A21 = A21_VC_Star;

        blas::Axpy( (F)-0.5, Y10, A10 );
	A10_Star_VR = A10;
	blas::internal::LocalTrsm
        ( Left, Lower, Normal, NonUnit,
          (F)1, L11_Star_Star, A10_Star_VR );

	A10_Star_MR = A10_Star_VR;
	A10 = A10_Star_MR;
        L21_MC_Star = L21;
        blas::internal::LocalGemm
        ( Normal, Normal,
          (F)1, L21_MC_Star, A10_Star_MR, (F)1, Y20 );

        // Symmetrize A11[* ,* ] by copying the lower triangle into the upper
        // so that we can call a local gemm instead of worrying about
        // reproducing a hemm with nonsymmetric local matrices.
        {
            const int height = A11_Star_Star.LocalHeight();
            const int ldim = A11_Star_Star.LocalLDim();
            F* A11Buffer = A11_Star_Star.LocalBuffer();
            for( int i=1; i<height; ++i )
                for( int j=0; j<i; ++j )
                    A11Buffer[j+i*ldim] = Conj(A11Buffer[i+j*ldim]);
        }
        A11_Star_MR = A11_Star_Star;
        blas::internal::LocalGemm
        ( Normal, Normal, (F)1, L21_MC_Star, A11_Star_MR, (F)0, Y21 );

        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
	  (F)1, L20, A10_Star_MR, (F)0, Z21_MC_Star );
        Y21.SumScatterUpdate( (F)1, Z21_MC_Star );
        //--------------------------------------------------------------------//
        A11_Star_MR.FreeAlignments();
        A21_VC_Star.FreeAlignments();
	A10_Star_VR.FreeAlignments();
        A10_Star_MR.FreeAlignments();
	L10_Star_VR.FreeAlignments();
	L10_Star_MR.FreeAlignments();
        L21_MC_Star.FreeAlignments();
        X21_MC_Star.FreeAlignments();
        Z21_MC_Star.FreeAlignments();

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

template void elemental::lapack::internal::HegstRLVar3
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::HegstRLVar3
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& L );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstRLVar3
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::HegstRLVar3
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& L );
#endif
