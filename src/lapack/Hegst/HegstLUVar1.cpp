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
elemental::lapack::internal::HegstLUVar1
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstLUVar1");
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
    DistMatrix<F,MC,  Star> A01_MC_Star(g);
    DistMatrix<F,MR,  Star> A01_MR_Star(g);
    DistMatrix<F,VC,  Star> A01_VC_Star(g);
    DistMatrix<F,Star,Star> U11_Star_Star(g);
    DistMatrix<F,MC,  Star> U01_MC_Star(g);
    DistMatrix<F,MR,  Star> U01_MR_Star(g);
    DistMatrix<F,VC,  Star> U01_VC_Star(g);
    DistMatrix<F,MC,  MR  > Y01(g);
    DistMatrix<F,VC,  Star> Y01_VC_Star(g);

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

        A01_MC_Star.AlignWith( A00 );
        A01_MR_Star.AlignWith( A00 );
        A01_VC_Star.AlignWith( A00 );
        U01_MC_Star.AlignWith( A00 );
        U01_MR_Star.AlignWith( A00 );
        U01_VC_Star.AlignWith( A00 );
        Y01.AlignWith( A01 );
        Y01_VC_Star.AlignWith( A01 );
        Y01_VC_Star.ResizeTo( A01.Height(), A01.Width() );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        U01_VC_Star = U01;
        blas::Hemm
        ( Right, Upper,
          (F)1, A11_Star_Star.LocalMatrix(), U01_VC_Star.LocalMatrix(),
          (F)0, Y01_VC_Star.LocalMatrix() );
        Y01 = Y01_VC_Star;

        blas::Trmm
        ( Left, Upper, Normal, NonUnit, (F)1, U00, A01 );

        blas::Axpy( (F)0.5, Y01, A01 );

        A01_MC_Star = A01;
        U01_MC_Star = U01;
        A01_VC_Star = A01_MC_Star;
        A01_MR_Star = A01_VC_Star;
        U01_MR_Star = U01_MC_Star;
        blas::internal::LocalTriangularRank2K
        ( Upper, ConjugateTranspose, ConjugateTranspose,
          (F)1, U01_MC_Star, A01_MC_Star, U01_MR_Star, A01_MR_Star,
          (F)1, A00 );

        blas::Axpy( (F)0.5, Y01_VC_Star, A01_VC_Star );
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Right, Upper, ConjugateTranspose, NonUnit,
          (F)1, U11_Star_Star, A01_VC_Star );
        A01 = A01_VC_Star;

        lapack::internal::LocalHegst
        ( Left, Upper, A11_Star_Star, U11_Star_Star );
        A11 = A11_Star_Star;
        //--------------------------------------------------------------------//
        A01_MC_Star.FreeAlignments();
        A01_MR_Star.FreeAlignments();
        A01_VC_Star.FreeAlignments();
        U01_MC_Star.FreeAlignments();
        U01_MR_Star.FreeAlignments();
        U01_VC_Star.FreeAlignments();
        Y01.FreeAlignments();
        Y01_VC_Star.FreeAlignments();

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

template void elemental::lapack::internal::HegstLUVar1
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::HegstLUVar1
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstLUVar1
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::HegstLUVar1
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& U );
#endif

