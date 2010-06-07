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
elemental::lapack::internal::HegstFalseU
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstFalseU");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( U.Height() != U.Width() )
        throw logic_error( "Triangular matrices must be square." );
    if( A.Height() != U.Height() )
        throw logic_error( "A and U must be the same size." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR>
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    // Temporary distributions
    DistMatrix<T,MC,  Star> A01_MC_Star(g);
    DistMatrix<T,VC,  Star> A01_VC_Star(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,Star,VR  > A12_Star_VR(g);
    DistMatrix<T,MR,  Star> U01_MR_Star(g);
    DistMatrix<T,MC,  Star> U01_MC_Star(g);
    DistMatrix<T,Star,Star> U11_Star_Star(g);
    DistMatrix<T,MR,  Star> E01_MR_Star(g);
    DistMatrix<T,MC,  Star> F01_MC_Star(g);
    DistMatrix<T,MR,  MC  > E01_MR_MC(g);
    DistMatrix<T,MC,  MR  > E01(g);
    DistMatrix<T,Star,MR  > G11_Star_MR(g);
    DistMatrix<T,MC,  MR  > G11(g);
    DistMatrix<T,Star,MR  > H12_Star_MR(g);

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
        U01_MR_Star.AlignWith( A00 );
        U01_MC_Star.AlignWith( A00 );
        E01_MR_Star.AlignWith( A00 );
        F01_MC_Star.AlignWith( A00 );
        E01.AlignWith( A01 );
        G11_Star_MR.AlignWith( U01 );
        G11.AlignWith( A11 );
        H12_Star_MR.AlignWith( A02 );
        E01_MR_Star.ResizeTo( A01.Height(), A01.Width() ); 
        F01_MC_Star.ResizeTo( A01.Height(), A01.Width() );
        G11_Star_MR.ResizeTo( A11.Height(), A11.Width() );
        H12_Star_MR.ResizeTo( A12.Height(), A12.Width() );
        E01_MR_Star.SetToZero();
        F01_MC_Star.SetToZero();
        //--------------------------------------------------------------------//
        U01_MC_Star = U01;
        U01_MR_Star = U01_MC_Star;
        blas::internal::LocalHemmAccumulateLU
        ( (T)1, A00, U01_MC_Star, U01_MR_Star, F01_MC_Star, E01_MR_Star );
        E01_MR_MC.SumScatterFrom( E01_MR_Star );
        E01 = E01_MR_MC;
        E01.SumScatterUpdate( (T)1, F01_MC_Star );

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal,
          (T)1, U01_MC_Star, A01, (T)0, G11_Star_MR );

        blas::Axpy( (T)-1, E01, A01 );
        A01_MC_Star = A01;
        
        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal,
          (T)1, A01_MC_Star, U01, (T)1, G11_Star_MR );
        G11.SumScatterFrom( G11_Star_MR );
        G11.MakeTrapezoidal( Left, Upper );
        blas::Axpy( (T)-1, G11, A11 );

        U11_Star_Star = U11;
        A01_VC_Star = A01_MC_Star;
        blas::internal::LocalTrsm
        ( Right, Upper, Normal, NonUnit, (T)1, U11_Star_Star, A01_VC_Star );
        A01 = A01_VC_Star;

        A11_Star_Star = A11;
        lapack::internal::LocalHegst
        ( false, Upper, A11_Star_Star, U11_Star_Star );
        A11 = A11_Star_Star;

        blas::internal::LocalGemm
        ( ConjugateTranspose, Normal,
          (T)1, U01_MC_Star, A02, (T)0, H12_Star_MR );
        A12.SumScatterUpdate( (T)-1, H12_Star_MR );

        A12_Star_VR = A12;
        blas::internal::LocalTrsm
        ( Left, Upper, ConjugateTranspose, NonUnit,
          (T)1, U11_Star_Star, A12_Star_VR );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_MC_Star.FreeAlignments();
        U01_MR_Star.FreeAlignments();
        U01_MC_Star.FreeAlignments();
        E01_MR_Star.FreeAlignments();
        F01_MC_Star.FreeAlignments();
        E01.FreeAlignments();
        G11_Star_MR.FreeAlignments();
        G11.FreeAlignments();
        H12_Star_MR.FreeAlignments();

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

template void elemental::lapack::internal::HegstFalseU
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::HegstFalseU
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstFalseU
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::HegstFalseU
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& U );
#endif

