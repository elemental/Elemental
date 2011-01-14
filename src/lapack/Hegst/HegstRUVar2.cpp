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

// This routine has only partially been optimized. The ReduceScatter operations
// need to be (conjugate-)transposed in order to play nice with cache.
template<typename F> // F represents a real or complex field
void
elemental::lapack::internal::HegstRUVar2
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstRUVar2");
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
    DistMatrix<F,MC,  Star> A01_MC_Star(g);
    DistMatrix<F,VC,  Star> A01_VC_Star(g);
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,Star,VR  > A12_Star_VR(g);
    DistMatrix<F,MC,  Star> U01_MC_Star(g);
    DistMatrix<F,VR,  Star> U01_VR_Star(g);
    DistMatrix<F,Star,MR  > U01Herm_Star_MR(g);
    DistMatrix<F,Star,Star> U11_Star_Star(g);
    DistMatrix<F,MR,  Star> E01_MR_Star(g);
    DistMatrix<F,MC,  Star> F01_MC_Star(g);
    DistMatrix<F,MR,  MC  > E01_MR_MC(g);
    DistMatrix<F,MC,  MR  > E01(g);
    DistMatrix<F,Star,MR  > G11_Star_MR(g);
    DistMatrix<F,MC,  MR  > G11(g);
    DistMatrix<F,MR,  Star> H12Herm_MR_Star(g);
    DistMatrix<F,MR,  MC  > H12Herm_MR_MC(g);

    Matrix<F> H12Local;


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

        A01_MC_Star.AlignWith( U01 );
        U01_MC_Star.AlignWith( A00 );
        U01_VR_Star.AlignWith( A00 );
        U01Herm_Star_MR.AlignWith( A00 );
        E01_MR_Star.AlignWith( A00 );
        F01_MC_Star.AlignWith( A00 );
        E01.AlignWith( A01 );
        G11_Star_MR.AlignWith( U01 );
        G11.AlignWith( A11 );
        H12Herm_MR_Star.AlignWith( A02 );
        H12Herm_MR_MC.AlignWith( A12 );
        E01_MR_Star.ResizeTo( A01.Height(), A01.Width() ); 
        F01_MC_Star.ResizeTo( A01.Height(), A01.Width() );
        G11_Star_MR.ResizeTo( A11.Height(), A11.Width() );
        H12Herm_MR_Star.ResizeTo( A12.Width(), A12.Height() );
        E01_MR_Star.SetToZero();
        F01_MC_Star.SetToZero();
        G11_Star_MR.SetToZero();
        H12Herm_MR_Star.SetToZero();
        //--------------------------------------------------------------------//
        U01_MC_Star = U01;
        U01_VR_Star = U01_MC_Star;
        U01Herm_Star_MR.ConjugateTransposeFrom( U01_VR_Star );
        blas::internal::LocalSymmetricAccumulateLU
        ( ConjugateTranspose, 
          (F)1, A00, U01_MC_Star, U01Herm_Star_MR, F01_MC_Star, E01_MR_Star );
        E01_MR_MC.SumScatterFrom( E01_MR_Star );
        E01 = E01_MR_MC;
        E01.SumScatterUpdate( (F)1, F01_MC_Star );

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal,
          (F)1, U01_MC_Star, A01, (F)0, G11_Star_MR );

        blas::Axpy( (F)-1, E01, A01 );
        A01_MC_Star = A01;
        
        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal,
          (F)1, A01_MC_Star, U01, (F)1, G11_Star_MR );
        G11.SumScatterFrom( G11_Star_MR );
        G11.MakeTrapezoidal( Left, Upper );
        blas::Axpy( (F)-1, G11, A11 );

        U11_Star_Star = U11;
        A01_VC_Star = A01_MC_Star;
        blas::internal::LocalTrsm
        ( Right, Upper, Normal, NonUnit, (F)1, U11_Star_Star, A01_VC_Star );
        A01 = A01_VC_Star;

        A11_Star_Star = A11;
        lapack::internal::LocalHegst
        ( Right, Upper, A11_Star_Star, U11_Star_Star );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal,
          (F)1, A02, U01_MC_Star, (F)0, H12Herm_MR_Star );
        H12Herm_MR_MC.SumScatterFrom( H12Herm_MR_Star );
        blas::ConjTrans( H12Herm_MR_MC.LockedLocalMatrix(), H12Local );
        blas::Axpy( (F)-1, H12Local, A12.LocalMatrix() );

        A12_Star_VR = A12;
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (F)1, U11_Star_Star, A12_Star_VR );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_MC_Star.FreeAlignments();
        U01_MC_Star.FreeAlignments();
        U01_VR_Star.FreeAlignments();
        U01Herm_Star_MR.FreeAlignments();
        E01_MR_Star.FreeAlignments();
        F01_MC_Star.FreeAlignments();
        E01.FreeAlignments();
        G11_Star_MR.FreeAlignments();
        G11.FreeAlignments();
        H12Herm_MR_Star.FreeAlignments();
        H12Herm_MR_MC.FreeAlignments();

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

template void elemental::lapack::internal::HegstRUVar2
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::HegstRUVar2
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstRUVar2
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::HegstRUVar2
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& U );
#endif

