/*
   Copyright (c) 2009-2010, Jack Poulson
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
template<typename T>
void
elemental::lapack::internal::HegstRL
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstRL");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( L.Height() != L.Width() )
        throw logic_error( "Triangular matrices must be square." );
    if( A.Height() != L.Height() )
        throw logic_error( "A and L must be the same size." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<T,MR,  Star> A10Herm_MR_Star(g);
    DistMatrix<T,Star,VR  > A10_Star_VR(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,VC,  Star> A21_VC_Star(g);
    DistMatrix<T,MR,  Star> L10Herm_MR_Star(g);
    DistMatrix<T,VC,  Star> L10Herm_VC_Star(g);
    DistMatrix<T,Star,MC  > L10_Star_MC(g);
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,MC,  Star> E10Herm_MC_Star(g);
    DistMatrix<T,MR,  Star> F10Herm_MR_Star(g);
    DistMatrix<T,MC,  MR  > E10Herm(g);
    DistMatrix<T,MR,  MC  > E10Herm_MR_MC(g);
    DistMatrix<T,MC,  Star> G11_MC_Star(g);
    DistMatrix<T,MC,  MR  > G11(g);
    DistMatrix<T,MC,  Star> H21_MC_Star(g);

    Matrix<T> E10Local;

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

        A10Herm_MR_Star.AlignWith( L10 );
        L10Herm_MR_Star.AlignWith( A00 );
        L10Herm_VC_Star.AlignWith( A00 );
        L10_Star_MC.AlignWith( A00 );
        E10Herm_MC_Star.AlignWith( A00 );
        F10Herm_MR_Star.AlignWith( A00 );
        E10Herm_MR_MC.AlignWith( A10 );
        G11_MC_Star.AlignWith( L10 );
        G11.AlignWith( A11 );
        H21_MC_Star.AlignWith( A20 );
        E10Herm_MC_Star.ResizeTo( A10.Width(), A10.Height() );
        F10Herm_MR_Star.ResizeTo( A10.Width(), A10.Height() );
        G11_MC_Star.ResizeTo( A11.Height(), A11.Width() );
        H21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        E10Herm_MC_Star.SetToZero();
        F10Herm_MR_Star.SetToZero();
        //--------------------------------------------------------------------//
        L10Herm_MR_Star.ConjugateTransposeFrom( L10 );
        L10Herm_VC_Star = L10Herm_MR_Star;
        L10_Star_MC.ConjugateTransposeFrom( L10Herm_VC_Star );
        blas::internal::LocalHemmAccumulateRL
        ( (T)1, A00, L10_Star_MC, L10Herm_MR_Star, 
          E10Herm_MC_Star, F10Herm_MR_Star );
        E10Herm.SumScatterFrom( E10Herm_MC_Star );
        E10Herm_MR_MC = E10Herm;
        E10Herm_MR_MC.SumScatterUpdate( (T)1, F10Herm_MR_Star );
        blas::ConjTrans( E10Herm_MR_MC.LockedLocalMatrix(), E10Local );

        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, A10, L10Herm_MR_Star, (T)0, G11_MC_Star );

        blas::Axpy( (T)-1, E10Local, A10.LocalMatrix() );
        A10Herm_MR_Star.ConjugateTransposeFrom( A10 );
        
        blas::internal::LocalGemm
        ( Normal, Normal,
          (T)1, L10, A10Herm_MR_Star, (T)1, G11_MC_Star );
        G11.SumScatterFrom( G11_MC_Star );
        G11.MakeTrapezoidal( Left, Lower );
        blas::Axpy( (T)-1, G11, A11 );

        L11_Star_Star = L11;
        A10_Star_VR.ConjugateTransposeFrom( A10Herm_MR_Star );
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, NonUnit, (T)1, L11_Star_Star, A10_Star_VR );
        A10 = A10_Star_VR;

        A11_Star_Star = A11;
        lapack::internal::LocalHegst
        ( Right, Lower, A11_Star_Star, L11_Star_Star );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( Normal, Normal,
          (T)1, A20, L10Herm_MR_Star, (T)0, H21_MC_Star );
        A21.SumScatterUpdate( (T)-1, H21_MC_Star );

        A21_VC_Star =  A21;
        blas::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (T)1, L11_Star_Star, A21_VC_Star );
        A21 = A21_VC_Star;
        //--------------------------------------------------------------------//
        A10Herm_MR_Star.FreeAlignments();
        L10Herm_MR_Star.FreeAlignments();
        L10Herm_VC_Star.FreeAlignments();
        L10_Star_MC.FreeAlignments();
        E10Herm_MC_Star.FreeAlignments();
        F10Herm_MR_Star.FreeAlignments();
        E10Herm_MR_MC.FreeAlignments();
        G11_MC_Star.FreeAlignments();
        G11.FreeAlignments();
        H21_MC_Star.FreeAlignments();

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

template<typename T>
void
elemental::lapack::internal::HegstRLNaive
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstRLNaive");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( L.Height() != L.Width() )
        throw logic_error( "Triangular matrices must be square." );
    if( A.Height() != L.Height() )
        throw logic_error( "A and L must be the same size." );
    if( A.GetGrid().VCRank() == 0 )
    { 
        cout << "HegstRLNaive exists solely for academic purposes. Please "
                "use HegstRL in real applications." << endl; 
    }
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR>
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    // Temporary distributions
    DistMatrix<T,Star,MR  > A10_Star_MR(g);
    DistMatrix<T,Star,VR  > A10_Star_VR(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,VC,  Star> A21_VC_Star(g);
    DistMatrix<T,Star,MR  > L10_Star_MR(g);
    DistMatrix<T,Star,MC  > L10_Star_MC(g);
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,Star,MC  > E10_Star_MC(g);
    DistMatrix<T,Star,MR  > F10_Star_MR(g);
    DistMatrix<T,MR,  MC  > E10_MR_MC(g);
    DistMatrix<T,MC,  MR  > E10(g);
    DistMatrix<T,MC,  Star> G11_MC_Star(g);
    DistMatrix<T,MC,  MR  > G11(g);
    DistMatrix<T,MC,  Star> H21_MC_Star(g);

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

        A10_Star_MR.AlignWith( L10 );
        L10_Star_MR.AlignWith( A00 );
        L10_Star_MC.AlignWith( A00 );
        E10_Star_MC.AlignWith( A00 );
        F10_Star_MR.AlignWith( A00 );
        E10.AlignWith( A10 );
        G11_MC_Star.AlignWith( L10 );
        G11.AlignWith( A11 );
        H21_MC_Star.AlignWith( A20 );
        E10_Star_MC.ResizeTo( A10.Height(), A10.Width() );
        F10_Star_MR.ResizeTo( A10.Height(), A10.Width() );
        G11_MC_Star.ResizeTo( A11.Height(), A11.Width() );
        H21_MC_Star.ResizeTo( A21.Height(), A21.Width() );
        E10_Star_MC.SetToZero();
        F10_Star_MR.SetToZero();
        //--------------------------------------------------------------------//
        L10_Star_MR = L10;
        L10_Star_MC = L10_Star_MR;
        blas::internal::LocalHemmAccumulateRL
        ( (T)1, A00, L10_Star_MC, L10_Star_MR, E10_Star_MC, F10_Star_MR );
        E10_MR_MC.SumScatterFrom( E10_Star_MC );
        E10 = E10_MR_MC;
        E10.SumScatterUpdate( (T)1, F10_Star_MR );
        
        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          (T)1, A10, L10_Star_MR, (T)0, G11_MC_Star );

        blas::Axpy( (T)-1, E10, A10 );
        A10_Star_MR = A10;
        
        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose,
          (T)1, L10, A10_Star_MR, (T)1, G11_MC_Star );
        G11.SumScatterFrom( G11_MC_Star );
        G11.MakeTrapezoidal( Left, Lower );
        blas::Axpy( (T)-1, G11, A11 );

        L11_Star_Star = L11;
        A10_Star_VR = A10_Star_MR;
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, NonUnit, (T)1, L11_Star_Star, A10_Star_VR );
        A10 = A10_Star_VR;

        A11_Star_Star = A11;
        lapack::internal::LocalHegst
        ( Right, Lower, A11_Star_Star, L11_Star_Star );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose,
          (T)1, A20, L10_Star_MR, (T)0, H21_MC_Star );
        A21.SumScatterUpdate( (T)-1, H21_MC_Star );

        A21_VC_Star =  A21;
        blas::internal::LocalTrsm
        ( Right, Lower, ConjugateTranspose, NonUnit, 
          (T)1, L11_Star_Star, A21_VC_Star );
        A21 = A21_VC_Star;
        //--------------------------------------------------------------------//
        A10_Star_MR.FreeAlignments();
        L10_Star_MR.FreeAlignments();
        L10_Star_MC.FreeAlignments();
        E10_Star_MC.FreeAlignments();
        F10_Star_MR.FreeAlignments();
        E10.FreeAlignments();
        G11_MC_Star.FreeAlignments();
        G11.FreeAlignments();
        H21_MC_Star.FreeAlignments();

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

template void elemental::lapack::internal::HegstRL
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::HegstRLNaive
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::HegstRL
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& L );

template void elemental::lapack::internal::HegstRLNaive
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& L );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstRL
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::HegstRLNaive
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::HegstRL
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& L );

template void elemental::lapack::internal::HegstRLNaive
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& L );
#endif

