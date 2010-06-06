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
elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvU");
#endif
    TrinvUVar3( diagonal, U );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvUVar3");
    if( U.Height() != U.Width() )
        throw logic_error( "Nonsquare matrices cannot be triangular." );
#endif
    const Grid& grid = U.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        UTL(grid), UTR(grid),  U00(grid), U01(grid), U02(grid),
        UBL(grid), UBR(grid),  U10(grid), U11(grid), U12(grid),
                               U20(grid), U21(grid), U22(grid);

    // Temporary distributions

    DistMatrix<T,VC,  Star> U01_VC_Star(grid);
    DistMatrix<T,Star,Star> U11_Star_Star(grid);
    DistMatrix<T,Star,VR  > U12_Star_VR(grid);
    DistMatrix<T,Star,MC  > U01Trans_Star_MC(grid);
    DistMatrix<T,MR,  Star> U12Trans_MR_Star(grid);

    // Start the algorithm
    PartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    while( UBR.Height() < U.Height() )
    {
        RepartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        U01Trans_Star_MC.AlignWith( U02 );
        U12Trans_MR_Star.AlignWith( U02 );
        //--------------------------------------------------------------------//
        U11_Star_Star = U11;
        lapack::Trinv( Upper, diagonal, U11_Star_Star.LocalMatrix() );
        U11 = U11_Star_Star;

        U01_VC_Star = U01;
        blas::internal::LocalTrmm
        ( Right, Upper, Normal, diagonal, (T)-1, U11_Star_Star, U01_VC_Star );

        // We transpose before the communication to avoid cache-thrashing
        // in the unpacking stage.
        U12Trans_MR_Star.TransposeFrom( U12 );
        U01Trans_Star_MC.TransposeFrom( U01_VC_Star );

        blas::internal::LocalGemm
        ( Transpose, Transpose, 
          (T)1, U01Trans_Star_MC, U12Trans_MR_Star, (T)1, U02 );
        U01.TransposeFrom( U01Trans_Star_MC );

        U12_Star_VR.TransposeFrom( U12Trans_MR_Star );
        blas::internal::LocalTrmm
        ( Left, Upper, Normal, diagonal, (T)1, U11_Star_Star, U12_Star_VR );
        U12 = U12_Star_VR;
        //--------------------------------------------------------------------//
        U01Trans_Star_MC.FreeAlignments();
        U12Trans_MR_Star.FreeAlignments();

        SlidePartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<float,MC,MR>& U );

template void elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<double,MC,MR>& U );

template void elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<scomplex,MC,MR>& U );

template void elemental::lapack::internal::TrinvU
( Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& U );

template void elemental::lapack::internal::TrinvUVar3
( Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& U );
#endif

