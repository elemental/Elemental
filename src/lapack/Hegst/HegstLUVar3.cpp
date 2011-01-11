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
elemental::lapack::internal::HegstLUVar3
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstLUVar3");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( U.Height() != U.Width() )
        throw logic_error( "Triangular matrices must be square." );
    if( A.Height() != U.Height() )
        throw logic_error( "A and U must be the same size." );
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<F,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<F,MC,MR>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,Star,VR  > A12_Star_VR(g);
    DistMatrix<F,Star,Star> U11_Star_Star(g);
    DistMatrix<F,Star,MC  > U12_Star_MC(g);
    DistMatrix<F,Star,VR  > U12_Star_VR(g);
    DistMatrix<F,MR,  Star> U12Herm_MR_Star(g);
    DistMatrix<F,VC,  Star> U12Herm_VC_Star(g);
    DistMatrix<F,MC,  MR  > Z12Herm(g);
    DistMatrix<F,MR,  MC  > Z12Herm_MR_MC(g);
    DistMatrix<F,MC,  Star> Z12Herm_MC_Star(g);
    DistMatrix<F,MR,  Star> Z12Herm_MR_Star(g);
    DistMatrix<F,MC,  MR  > Y12(g);
    DistMatrix<F,Star,Star> X11_Star_Star(g);

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

        A12_Star_VR.AlignWith( A12 );
        U12_Star_MC.AlignWith( A22 );
        U12_Star_VR.AlignWith( A12 );
        U12Herm_MR_Star.AlignWith( A22 );
        U12Herm_VC_Star.AlignWith( A22 );
        Z12Herm.AlignWith( A12 );
        Z12Herm_MR_MC.AlignWith( A12 );
        Z12Herm_MC_Star.AlignWith( A22 );
        Z12Herm_MR_Star.AlignWith( A22 );
        Y12.AlignWith( A12 );
        Z12Herm_MC_Star.ResizeTo( A12.Width(), A12.Height() );
        Z12Herm_MR_Star.ResizeTo( A12.Width(), A12.Height() );
        Y12.ResizeTo( A12.Height(), A12.Width() );
        X11_Star_Star.ResizeTo( A11.Height(), A11.Width() );
        //--------------------------------------------------------------------//
        U12Herm_MR_Star.ConjugateTransposeFrom( U12 );
        U12Herm_VC_Star = U12Herm_MR_Star;
        U12_Star_MC.ConjugateTransposeFrom( U12Herm_VC_Star );
        Z12Herm_MC_Star.SetToZero();
        Z12Herm_MR_Star.SetToZero();
        blas::internal::LocalHemmAccumulateRU
        ( (F)1, A22, U12_Star_MC, U12Herm_MR_Star, 
          Z12Herm_MC_Star, Z12Herm_MR_Star );
        Z12Herm.SumScatterFrom( Z12Herm_MC_Star );
        Z12Herm_MR_MC = Z12Herm;
        Z12Herm_MR_MC.SumScatterUpdate( (F)1, Z12Herm_MR_Star );
        blas::ConjTrans( Z12Herm_MR_MC.LockedLocalMatrix(), Y12.LocalMatrix() );

        A12_Star_VR = A12;
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Left, Upper, Normal, NonUnit, (F)1, U11_Star_Star, A12_Star_VR );
        A12 = A12_Star_VR;

        blas::Axpy( (F)0.5, Y12, A12 );

        A11_Star_Star = A11;
        lapack::internal::LocalHegst
        ( Left, Upper, A11_Star_Star, U11_Star_Star );
        A11 = A11_Star_Star;

        A12_Star_VR = A12;
        U12_Star_VR = U12;
        blas::Her2k
        ( Upper, Normal,
          (F)1, A12_Star_VR.LocalMatrix(), U12_Star_VR.LocalMatrix(),
          (F)0, X11_Star_Star.LocalMatrix() );
        A11.SumScatterUpdate( (F)1, X11_Star_Star );

        blas::Axpy( (F)0.5, Y12, A12 );

        // This is the bottleneck because of A12 having only blocksize rows
        blas::Trmm
        ( Right, Upper, ConjugateTranspose, NonUnit, (F)1, U22, A12 );
        //--------------------------------------------------------------------//
        A12_Star_VR.FreeAlignments();
        U12_Star_MC.FreeAlignments();
        U12_Star_VR.FreeAlignments();
        U12Herm_MR_Star.FreeAlignments();
        U12Herm_VC_Star.FreeAlignments();
        Z12Herm.FreeAlignments();
        Z12Herm_MR_MC.FreeAlignments(); 
        Z12Herm_MC_Star.FreeAlignments();
        Z12Herm_MR_Star.FreeAlignments();
        Y12.FreeAlignments();

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

template void elemental::lapack::internal::HegstLUVar3
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::HegstLUVar3
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstLUVar3
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::HegstLUVar3
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& U );
#endif

