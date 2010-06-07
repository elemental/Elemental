/*
   This file is part of elemental, a library for distributed-memory dense
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::lapack::internal::HegstFalseL
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstFalseL");
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
    DistMatrix<T,Star,MR  > A10_Star_MR(g);
    DistMatrix<T,Star,VR  > A10_Star_VR(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,VC,  Star> A21_VC_Star(g);
    DistMatrix<T,Star,MC  > L10_Star_MC(g);
    DistMatrix<T,Star,MR  > L10_Star_MR(g);
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
        L10_Star_MC.AlignWith( A00 );
        L10_Star_MR.AlignWith( A00 );
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
        ( false, Lower, A11_Star_Star, L11_Star_Star );
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
        L10_Star_MC.FreeAlignments();
        L10_Star_MR.FreeAlignments();
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

template void elemental::lapack::internal::HegstFalseL
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::HegstFalseL
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& L );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstFalseL
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::HegstFalseL
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& L );
#endif

