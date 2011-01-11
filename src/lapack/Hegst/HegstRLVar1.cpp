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
elemental::lapack::internal::HegstRLVar1
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstRLVar1");
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
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,Star,VR  > A10_Star_VR(g);
    DistMatrix<F,Star,Star> L11_Star_Star(g);
    DistMatrix<F,Star,MC  > L10_Star_MC(g);
    DistMatrix<F,Star,VR  > L10_Star_VR(g);
    DistMatrix<F,MR,  Star> L10Herm_MR_Star(g);
    DistMatrix<F,VC,  Star> L10Herm_VC_Star(g);
    DistMatrix<F,MC,  MR  > Z10Herm(g);
    DistMatrix<F,MR,  MC  > Z10Herm_MR_MC(g);
    DistMatrix<F,MC,  Star> Z10Herm_MC_Star(g);
    DistMatrix<F,MR,  Star> Z10Herm_MR_Star(g);
    DistMatrix<F,MC,  MR  > Y10(g);
    DistMatrix<F,Star,Star> X11_Star_Star(g);

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

        A10_Star_VR.AlignWith( A10 );
        L10_Star_MC.AlignWith( A00 );
        L10_Star_VR.AlignWith( A00 );
        L10Herm_MR_Star.AlignWith( A00 );
        L10Herm_VC_Star.AlignWith( A00 );
        Z10Herm.AlignWith( A10 );
        Z10Herm_MR_MC.AlignWith( A10 );
        Z10Herm_MC_Star.AlignWith( A00 );
        Z10Herm_MR_Star.AlignWith( A00 );
        Y10.AlignWith( A10 );
        Z10Herm_MC_Star.ResizeTo( A10.Width(), A10.Height() );
        Z10Herm_MR_Star.ResizeTo( A10.Width(), A10.Height() );
        Y10.ResizeTo( A10.Height(), A10.Width() );
        X11_Star_Star.ResizeTo( A11.Height(), A11.Width() );
        //--------------------------------------------------------------------//
        L10Herm_MR_Star.ConjugateTransposeFrom( L10 );
        L10Herm_VC_Star = L10Herm_MR_Star;
        L10_Star_MC.ConjugateTransposeFrom( L10Herm_VC_Star );
        Z10Herm_MC_Star.SetToZero();
        Z10Herm_MR_Star.SetToZero();
        blas::internal::LocalHemmAccumulateRL
        ( (F)1, A00, L10_Star_MC, L10Herm_MR_Star, 
          Z10Herm_MC_Star, Z10Herm_MR_Star );
        Z10Herm.SumScatterFrom( Z10Herm_MC_Star );
        Z10Herm_MR_MC = Z10Herm;
        Z10Herm_MR_MC.SumScatterUpdate( (F)1, Z10Herm_MR_Star );
        blas::ConjTrans( Z10Herm_MR_MC.LocalMatrix(), Y10.LocalMatrix() );

        // This is the bottleneck because A10 only has blocksize rows
        blas::Trsm
        ( Right, Lower, ConjugateTranspose, NonUnit, (F)1, L00, A10 );

        blas::Axpy( (F)-0.5, Y10, A10 );

        A10_Star_VR = A10;
        L10_Star_VR = L10;
        blas::Her2k
        ( Lower, Normal,
          (F)-1, A10_Star_VR.LocalMatrix(), L10_Star_VR.LocalMatrix(), 
          (F)0, X11_Star_Star.LocalMatrix() );
        A11.SumScatterUpdate( (F)1, X11_Star_Star );

        A11_Star_Star = A11;
        L11_Star_Star = L11;
        lapack::internal::LocalHegst
        ( Right, Lower, A11_Star_Star, L11_Star_Star );
        A11 = A11_Star_Star;

        blas::Axpy( (F)-0.5, Y10, A10 );

        A10_Star_VR = A10;
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, NonUnit, (F)1, L11_Star_Star, A10_Star_VR );
        A10 = A10_Star_VR;
        //--------------------------------------------------------------------//
        A10_Star_VR.FreeAlignments();
        L10_Star_MC.FreeAlignments();
        L10_Star_VR.FreeAlignments();
        L10Herm_MR_Star.FreeAlignments();
        L10Herm_VC_Star.FreeAlignments();
        Z10Herm.FreeAlignments();
        Z10Herm_MR_MC.FreeAlignments();
        Z10Herm_MC_Star.FreeAlignments();
        Z10Herm_MR_Star.FreeAlignments();
        Y10.FreeAlignments();

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

template void elemental::lapack::internal::HegstRLVar1
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::HegstRLVar1
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& L );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstRLVar1
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::HegstRLVar1
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& L );
#endif

