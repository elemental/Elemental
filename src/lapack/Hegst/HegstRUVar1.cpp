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
elemental::lapack::internal::HegstRUVar1
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstRUVar1");
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
    DistMatrix<F,VC,  Star> A01_VC_Star(g);
    DistMatrix<F,Star,Star> U11_Star_Star(g);
    DistMatrix<F,MC,  Star> U01_MC_Star(g);
    DistMatrix<F,VC,  Star> U01_VC_Star(g);
    DistMatrix<F,VR,  Star> U01_VR_Star(g);
    DistMatrix<F,Star,MR  > U01Herm_Star_MR(g);
    DistMatrix<F,Star,Star> X11_Star_Star(g);
    DistMatrix<F,MC,  MR  > Y01(g);
    DistMatrix<F,MR,  MC  > Z01_MR_MC(g);
    DistMatrix<F,MC,  Star> Z01_MC_Star(g);
    DistMatrix<F,MR,  Star> Z01_MR_Star(g);

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

        A01_VC_Star.AlignWith( A01 );
        U01_MC_Star.AlignWith( A00 );
        U01_VR_Star.AlignWith( A00 );
        U01_VC_Star.AlignWith( A00 );
        U01Herm_Star_MR.AlignWith( A00 );
        Y01.AlignWith( A01 );
        Z01_MR_MC.AlignWith( A01 );
        Z01_MC_Star.AlignWith( A00 );
        Z01_MR_Star.AlignWith( A00 );
        //--------------------------------------------------------------------//
        // Y01 := A00 U01
        U01_MC_Star = U01;
        U01_VR_Star = U01_MC_Star;
        U01Herm_Star_MR.ConjugateTransposeFrom( U01_VR_Star );
        Z01_MC_Star.ResizeTo( A01.Height(), A01.Width() );
        Z01_MR_Star.ResizeTo( A01.Height(), A01.Width() );
        Z01_MC_Star.SetToZero();
        Z01_MR_Star.SetToZero();
        blas::internal::LocalSymmetricAccumulateLU
        ( ConjugateTranspose, 
          (F)1, A00, U01_MC_Star, U01Herm_Star_MR, Z01_MC_Star, Z01_MR_Star );
        Z01_MR_MC.SumScatterFrom( Z01_MR_Star );
        Y01 = Z01_MR_MC;
        Y01.SumScatterUpdate( (F)1, Z01_MC_Star );

        // A01 := inv(U00)' A01
        //
        // This is the bottleneck because A01 only has blocksize columns
        blas::Trsm
        ( Left, Upper, ConjugateTranspose, NonUnit, (F)1, U00, A01 );

        // A01 := A01 - 1/2 Y01
        blas::Axpy( (F)-0.5, Y01, A01 );

        // A11 := A11 - (U01' A01 + A01' U01)
        A01_VC_Star = A01;
        U01_VC_Star = U01_MC_Star;
        X11_Star_Star.ResizeTo( A11.Height(), A11.Width() );
        blas::Her2k
        ( Upper, ConjugateTranspose,
          (F)-1, A01_VC_Star.LocalMatrix(), U01_VC_Star.LocalMatrix(),
          (F)0, X11_Star_Star.LocalMatrix() );
        A11.SumScatterUpdate( (F)1, X11_Star_Star );

        // A11 := inv(U11)' A11 inv(U11)
        A11_Star_Star = A11;
        U11_Star_Star = U11;
        lapack::internal::LocalHegst
        ( Right, Upper, A11_Star_Star, U11_Star_Star );
        A11 = A11_Star_Star;

        // A01 := A01 - 1/2 Y01
        blas::Axpy( (F)-0.5, Y01, A01 );

        // A01 := A01 inv(U11)
        A01_VC_Star = A01;
        blas::internal::LocalTrsm
        ( Right, Upper, Normal, NonUnit, (F)1, U11_Star_Star, A01_VC_Star );
        A01 = A01_VC_Star;
        //--------------------------------------------------------------------//
        A01_VC_Star.FreeAlignments();
        U01_MC_Star.FreeAlignments();
        U01_VR_Star.FreeAlignments();
        U01_VC_Star.FreeAlignments();
        U01Herm_Star_MR.FreeAlignments();
        Y01.FreeAlignments();
        Z01_MR_MC.FreeAlignments();
        Z01_MC_Star.FreeAlignments();
        Z01_MR_Star.FreeAlignments();

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

template void elemental::lapack::internal::HegstRUVar1
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::HegstRUVar1
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstRUVar1
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::HegstRUVar1
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& U );
#endif

