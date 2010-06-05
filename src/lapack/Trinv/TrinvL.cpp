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
elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvL");
#endif
    TrinvLVar3( diagonal, L );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TrinvLVar3");
    if( L.Height() != L.Width() )
        throw "Nonsquare matrices cannot be triangular.";
#endif
    const Grid& grid = L.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(grid), LTR(grid),  L00(grid), L01(grid), L02(grid),
        LBL(grid), LBR(grid),  L10(grid), L11(grid), L12(grid),
                               L20(grid), L21(grid), L22(grid);

    // Temporary distributions
    DistMatrix<T,Star,MR  > L10_Star_MR(grid);
    DistMatrix<T,Star,VR  > L10_Star_VR(grid);
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,MC,  Star> L21_MC_Star(grid);
    DistMatrix<T,VC,  Star> L21_VC_Star(grid);

    // Start the algorithm
    PartitionDownDiagonal( L, LTL, LTR,
                              LBL, LBR );
    while( LTL.Height() < L.Height() )
    {
        RepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        L10_Star_MR.AlignWith( L20 );
        L21_MC_Star.AlignWith( L20 );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11;
        lapack::Trinv( Lower, diagonal, L11_Star_Star.LocalMatrix() );
        L11 = L11_Star_Star;

        L10_Star_VR = L10;
        blas::Trmm
        ( Left, Lower, Normal, diagonal,
          (T)-1, L11_Star_Star.LockedLocalMatrix(),
                 L10_Star_VR.LocalMatrix()         );

        L21_MC_Star = L21;
        L10_Star_MR = L10_Star_VR;
        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, L21_MC_Star, L10_Star_MR, (T)1, L20 );
        L10 = L10_Star_MR;

        L21_VC_Star = L21_MC_Star;
        blas::Trmm
        ( Right, Lower, Normal, diagonal,
          (T)1, L11_Star_Star.LockedLocalMatrix(),
                L21_VC_Star.LocalMatrix()         );
        L21 = L21_VC_Star;
        //--------------------------------------------------------------------//
        L10_Star_MR.FreeAlignments();
        L21_MC_Star.FreeAlignments();

        SlidePartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<float,MC,MR>& L );

template void elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<double,MC,MR>& L );

template void elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<double,MC,MR>& L );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<scomplex,MC,MR>& L );

template void elemental::lapack::internal::TrinvL
( Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& L );

template void elemental::lapack::internal::TrinvLVar3
( Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& L );
#endif

