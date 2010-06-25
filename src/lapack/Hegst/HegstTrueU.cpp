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
elemental::lapack::internal::HegstTrueU
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstTrueU");
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
    DistMatrix<T,VC,  Star> A01_VC_Star(g);
    DistMatrix<T,VR,  Star> A01_VR_Star(g);
    DistMatrix<T,Star,MC  > A01Herm_Star_MC(g);
    DistMatrix<T,Star,MR  > A01Herm_Star_MR(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,Star,VR  > A12_Star_VR(g);
    DistMatrix<T,MR,  Star> A12Herm_MR_Star(g);
    DistMatrix<T,VC,  Star> U01_VC_Star(g);
    DistMatrix<T,VR,  Star> U01_VR_Star(g);
    DistMatrix<T,Star,MC  > U01Herm_Star_MC(g);
    DistMatrix<T,Star,MR  > U01Herm_Star_MR(g);
    DistMatrix<T,Star,Star> U11_Star_Star(g);
    DistMatrix<T,VC,  Star> X01_VC_Star(g);

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

        A01_VC_Star.AlignWith( A00 );
        A01_VR_Star.AlignWith( A00 );
        A01Herm_Star_MC.AlignWith( A00 );
        A01Herm_Star_MR.AlignWith( A00 );
        A12Herm_MR_Star.AlignWith( A02 );
        U01_VC_Star.AlignWith( A00 );
        U01_VR_Star.AlignWith( A00 );
        U01Herm_Star_MC.AlignWith( A00 );
        U01Herm_Star_MR.AlignWith( A00 );
        X01_VC_Star.AlignWith( A01 );
        X01_VC_Star.ResizeTo( A01.Height(), A01.Width() );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        U01_VC_Star = U01;
        blas::Hemm
        ( Right, Upper, 
          (T)0.5, A11_Star_Star.LockedLocalMatrix(), 
                  U01_VC_Star.LockedLocalMatrix(), 
          (T)0, X01_VC_Star.LocalMatrix() );

        A01_VC_Star = A01;
        blas::Axpy( (T)1, X01_VC_Star, A01_VC_Star );


        A01Herm_Star_MC.ConjugateTransposeFrom( A01_VC_Star );
        A01_VR_Star = A01_VC_Star;
        A01Herm_Star_MR.ConjugateTransposeFrom( A01_VR_Star );

        U01Herm_Star_MC.ConjugateTransposeFrom( U01_VC_Star );
        U01_VR_Star = U01_VC_Star;
        U01Herm_Star_MR.ConjugateTransposeFrom( U01_VR_Star );

        blas::internal::LocalTriangularRank2K
        ( Upper, ConjugateTranspose, ConjugateTranspose, 
          (T)1, A01Herm_Star_MC, U01Herm_Star_MC, 
                A01Herm_Star_MR, U01Herm_Star_MR,
          (T)1, A00 );

        blas::Axpy( (T)1, X01_VC_Star, A01_VC_Star );
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Right, Upper, ConjugateTranspose, NonUnit,
          (T)1, U11_Star_Star, A01_VC_Star );
        A01 = A01_VC_Star;

        A12Herm_MR_Star.ConjugateTransposeFrom( A12 );
        blas::internal::LocalGemm
        ( ConjugateTranspose, ConjugateTranspose, 
          (T)1, U01Herm_Star_MC, A12Herm_MR_Star, (T)1, A02 );

        lapack::internal::LocalHegst
        ( true, Upper, A11_Star_Star, U11_Star_Star );
        A11 = A11_Star_Star;

        A12_Star_VR.ConjugateTransposeFrom( A12Herm_MR_Star );
        blas::internal::LocalTrmm
        ( Left, Upper, Normal, NonUnit, (T)1, U11_Star_Star, A12_Star_VR );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_VC_Star.FreeAlignments();
        A01_VR_Star.FreeAlignments();
        A01Herm_Star_MC.FreeAlignments();
        A01Herm_Star_MR.FreeAlignments();
        A12Herm_MR_Star.FreeAlignments();
        U01_VC_Star.FreeAlignments();
        U01_VR_Star.FreeAlignments();
        U01Herm_Star_MC.FreeAlignments();
        U01Herm_Star_MR.FreeAlignments();
        X01_VC_Star.FreeAlignments();

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

template<typename T>
void
elemental::lapack::internal::HegstTrueUNaive
( DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::HegstTrueUNaive");
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( U.Height() != U.Width() )
        throw logic_error( "Triangular matrices must be square." );
    if( A.Height() != U.Height() )
        throw logic_error( "A and U must be the same size." );
    if( A.GetGrid().VCRank() == 0 )
    {
        cout << "HegstTrueUNaive exists solely for academic purposes. Please "
                "use HegstTrueU for real applications." << endl;
    }
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
    DistMatrix<T,VC,  Star> A01_VC_Star(g);
    DistMatrix<T,MC,  Star> A01_MC_Star(g);
    DistMatrix<T,MR,  Star> A01_MR_Star(g);
    DistMatrix<T,Star,Star> A11_Star_Star(g);
    DistMatrix<T,Star,VR  > A12_Star_VR(g);
    DistMatrix<T,Star,MR  > A12_Star_MR(g);
    DistMatrix<T,VC,  Star> U01_VC_Star(g);
    DistMatrix<T,MC,  Star> U01_MC_Star(g);
    DistMatrix<T,MR,  Star> U01_MR_Star(g);
    DistMatrix<T,Star,Star> U11_Star_Star(g);
    DistMatrix<T,VC,  Star> X01_VC_Star(g);

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

        A01_VC_Star.AlignWith( A00 );
        A01_MC_Star.AlignWith( A00 );
        A01_MR_Star.AlignWith( A00 );
        A12_Star_MR.AlignWith( A02 );
        U01_VC_Star.AlignWith( A00 );
        U01_MC_Star.AlignWith( A00 );
        U01_MR_Star.AlignWith( A00 );
        X01_VC_Star.AlignWith( A01 );
        X01_VC_Star.ResizeTo( A01.Height(), A01.Width() );
        //--------------------------------------------------------------------//
        A11_Star_Star = A11;
        U01_VC_Star = U01;
        blas::Hemm
        ( Right, Upper, 
          (T)0.5, A11_Star_Star.LockedLocalMatrix(), 
                  U01_VC_Star.LockedLocalMatrix(), 
          (T)0, X01_VC_Star.LocalMatrix() );

        A01_VC_Star = A01;
        blas::Axpy( (T)1, X01_VC_Star, A01_VC_Star );

        A01_MC_Star = A01_VC_Star;
        A01_MR_Star = A01_VC_Star;
        U01_MC_Star = U01_VC_Star;
        U01_MR_Star = U01_VC_Star;
        blas::internal::LocalTriangularRank2K
        ( Upper, ConjugateTranspose, ConjugateTranspose,
          (T)1, A01_MC_Star, U01_MC_Star, A01_MR_Star, U01_MR_Star, (T)1, A00 );

        blas::Axpy( (T)1, X01_VC_Star, A01_VC_Star );
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Right, Upper, ConjugateTranspose, NonUnit,
          (T)1, U11_Star_Star, A01_VC_Star );
        A01 = A01_VC_Star;

        A12_Star_MR = A12;
        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, U01_MC_Star, A12_Star_MR, (T)1, A02 );

        lapack::internal::LocalHegst
        ( true, Upper, A11_Star_Star, U11_Star_Star );
        A11 = A11_Star_Star;

        A12_Star_VR = A12_Star_MR;
        blas::internal::LocalTrmm
        ( Left, Upper, Normal, NonUnit, (T)1, U11_Star_Star, A12_Star_VR );
        A12 = A12_Star_VR;
        //--------------------------------------------------------------------//
        A01_VC_Star.FreeAlignments();
        A01_MC_Star.FreeAlignments();
        A01_MR_Star.FreeAlignments();
        A12_Star_MR.FreeAlignments();
        U01_VC_Star.FreeAlignments();
        U01_MC_Star.FreeAlignments();
        U01_MR_Star.FreeAlignments();
        X01_VC_Star.FreeAlignments();

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

template void elemental::lapack::internal::HegstTrueU
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::HegstTrueUNaive
( DistMatrix<float,MC,MR>& A, const DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::HegstTrueU
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& U );

template void elemental::lapack::internal::HegstTrueUNaive
( DistMatrix<double,MC,MR>& A, const DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::HegstTrueU
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::HegstTrueUNaive
( DistMatrix<scomplex,MC,MR>& A, const DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::HegstTrueU
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& U );

template void elemental::lapack::internal::HegstTrueUNaive
( DistMatrix<dcomplex,MC,MR>& A, const DistMatrix<dcomplex,MC,MR>& U );
#endif

