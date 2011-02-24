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
#include "elemental/basic_internal.hpp"
#include "elemental/advanced_internal.hpp"
using namespace std;
using namespace elemental;

template<typename F> // F represents a real or complex field
void
elemental::advanced::internal::HegstRUVar4
( DistMatrix<F,MC,MR>& A, const DistMatrix<F,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("advanced::internal::HegstRUVar4");
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
    DistMatrix<F,VC,  Star> A01_VC_Star(g);
    DistMatrix<F,Star,MC  > A01Trans_Star_MC(g);
    DistMatrix<F,Star,Star> A11_Star_Star(g);
    DistMatrix<F,Star,VR  > A12_Star_VR(g);
    DistMatrix<F,Star,VC  > A12_Star_VC(g);
    DistMatrix<F,Star,MC  > A12_Star_MC(g);
    DistMatrix<F,Star,MR  > A12_Star_MR(g);
    DistMatrix<F,Star,Star> U11_Star_Star(g);
    DistMatrix<F,MR,  Star> U12Trans_MR_Star(g);
    DistMatrix<F,VR,  Star> U12Trans_VR_Star(g);
    DistMatrix<F,Star,VR  > U12_Star_VR(g);
    DistMatrix<F,Star,VC  > U12_Star_VC(g);
    DistMatrix<F,Star,MC  > U12_Star_MC(g);
    DistMatrix<F,Star,VR  > Y12_Star_VR(g);

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

        A01_VC_Star.AlignWith( A02 );
        A01Trans_Star_MC.AlignWith( A02 );
        A12_Star_VR.AlignWith( A22 );
        A12_Star_VC.AlignWith( A22 );
        A12_Star_MC.AlignWith( A22 );
        A12_Star_MR.AlignWith( A22 );
        U12Trans_MR_Star.AlignWith( A02 );
        U12Trans_VR_Star.AlignWith( A02 );
        U12_Star_VR.AlignWith( A02 );
        U12_Star_VC.AlignWith( A22 );
        U12_Star_MC.AlignWith( A22 );
        Y12_Star_VR.AlignWith( A12 );
        //--------------------------------------------------------------------//
        // A01 := A01 inv(U11)
        A01_VC_Star = A01;
        U11_Star_Star = U11;
        basic::internal::LocalTrsm
        ( Right, Upper, Normal, NonUnit,
          (F)1, U11_Star_Star, A01_VC_Star );
        A01 = A01_VC_Star;

        // A11 := inv(U11)' A11 inv(U11)
        A11_Star_Star = A11;
        advanced::internal::LocalHegst
        ( Right, Upper, A11_Star_Star, U11_Star_Star );
        A11 = A11_Star_Star;

        // A02 := A02 - A01 U12
        A01Trans_Star_MC.TransposeFrom( A01_VC_Star );
        U12Trans_MR_Star.TransposeFrom( U12 );
        basic::internal::LocalGemm
        ( Transpose, Transpose, 
          (F)-1, A01Trans_Star_MC, U12Trans_MR_Star, (F)1, A02 );

        // Y12 := A11 U12
        U12Trans_VR_Star = U12Trans_MR_Star;
        U12_Star_VR.ResizeTo( A12.Height(), A12.Width() );
        basic::Trans
        ( U12Trans_VR_Star.LocalMatrix(), U12_Star_VR.LocalMatrix() );
        Y12_Star_VR.ResizeTo( A12.Height(), A12.Width() );
        basic::Hemm
        ( Left, Upper, 
          (F)1, A11_Star_Star.LocalMatrix(), U12_Star_VR.LocalMatrix(), 
          (F)0, Y12_Star_VR.LocalMatrix() );

        // A12 := inv(U11)' A12
        A12_Star_VR = A12;
        basic::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (F)1, U11_Star_Star, A12_Star_VR );

        // A12 := A12 - 1/2 Y12
        basic::Axpy( (F)-0.5, Y12_Star_VR, A12_Star_VR );

        // A22 := A22 - (A12 U12' + U12 A12')
        A12_Star_MR = A12_Star_VR;
        A12_Star_VC = A12_Star_VR;
        U12_Star_VC = U12_Star_VR;
        A12_Star_MC = A12_Star_VC;
        U12_Star_MC = U12_Star_VC;
        basic::internal::LocalTriangularRank2K
        ( Upper, ConjugateTranspose, ConjugateTranspose, Transpose,
          (F)-1, A12_Star_MC, U12_Star_MC, A12_Star_MR, U12Trans_MR_Star,
          (F)1, A22 );

        // A12 := A12 - 1/2 Y12
        basic::Axpy( (F)-0.5, Y12_Star_VR, A12_Star_VR );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_VC_Star.FreeAlignments();
        A01Trans_Star_MC.FreeAlignments();
        A12_Star_VR.FreeAlignments();
        A12_Star_VC.FreeAlignments();
        A12_Star_MC.FreeAlignments();
        A12_Star_MR.FreeAlignments();
        U12Trans_MR_Star.FreeAlignments();
        U12Trans_VR_Star.FreeAlignments();
        U12_Star_VR.FreeAlignments();
        U12_Star_VC.FreeAlignments();
        U12_Star_MC.FreeAlignments();
        Y12_Star_VR.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /**********************************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::advanced::internal::HegstRUVar4
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& U );

template void elemental::advanced::internal::HegstRUVar4
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::advanced::internal::HegstRUVar4
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& U );

template void elemental::advanced::internal::HegstRUVar4
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& U );
#endif

