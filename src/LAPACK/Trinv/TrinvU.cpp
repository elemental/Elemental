/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "Elemental/LAPACKInternal.hpp"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::LAPACK::Internal::TrinvU
( const Diagonal diagonal, DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::TrinvU");
#endif
    TrinvUVar3( diagonal, U );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::LAPACK::Internal::TrinvUVar3
( const Diagonal diagonal, DistMatrix<T,MC,MR>& U )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::TrinvUVar3");
    if( U.Height() != U.Width() )
        throw "Nonsquare matrices cannot be triangular.";
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
    PartitionUpDiagonal( U, UTL, UTR,
                            UBL, UBR );
    while( UBR.Height() < U.Height() )
    {
        RepartitionUpDiagonal( UTL, /**/ UTR,  U00, U01, /**/ U02,
                                    /**/       U10, U11, /**/ U12,
                              /*************/ /******************/
                               UBL, /**/ UBR,  U20, U21, /**/ U22 );

        U01Trans_Star_MC.AlignWith( U02 );
        U12Trans_MR_Star.AlignWith( U02 );
        //--------------------------------------------------------------------//
        U11_Star_Star = U11;
        LAPACK::Trinv( Upper, diagonal, U11_Star_Star.LocalMatrix() );
        U11 = U11_Star_Star;

        U01_VC_Star = U01;
        BLAS::Trmm
        ( Right, Upper, Normal, diagonal,
          (T)-1, U11_Star_Star.LockedLocalMatrix(),
                 U01_VC_Star.LocalMatrix()         );

        // We transpose before the communication to avoid cache-thrashing
        // in the unpacking stage.
        U12Trans_MR_Star.TransposeFrom( U12 );
        U01Trans_Star_MC.TransposeFrom( U01_VC_Star );

        BLAS::Gemm
        ( Transpose, Transpose,
          (T)1, U01Trans_Star_MC.LockedLocalMatrix(),
                U12Trans_MR_Star.LockedLocalMatrix(),
          (T)1, U02.LocalMatrix()                    );
        U01.TransposeFrom( U01Trans_Star_MC );

        U12_Star_VR.TransposeFrom( U12Trans_MR_Star );
        BLAS::Trmm
        ( Left, Upper, Normal, diagonal,
          (T)1, U11_Star_Star.LockedLocalMatrix(),
                U12_Star_VR.LocalMatrix()         );
        U12 = U12_Star_VR;
        //--------------------------------------------------------------------//
        U01Trans_Star_MC.FreeConstraints();
        U12Trans_MR_Star.FreeConstraints();

        SlidePartitionUpDiagonal( UTL, /**/ UTR,  U00, /**/ U01, U02,
                                 /*************/ /******************/
                                       /**/       U10, /**/ U11, U12,
                                  UBL, /**/ UBR,  U20, /**/ U21, U22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::LAPACK::Internal::TrinvU
( const Diagonal diagonal, DistMatrix<float,MC,MR>& U );

template void Elemental::LAPACK::Internal::TrinvUVar3
( const Diagonal diagonal, DistMatrix<float,MC,MR>& U );

template void Elemental::LAPACK::Internal::TrinvU
( const Diagonal diagonal, DistMatrix<double,MC,MR>& U );

template void Elemental::LAPACK::Internal::TrinvUVar3
( const Diagonal diagonal, DistMatrix<double,MC,MR>& U );

#ifndef WITHOUT_COMPLEX
template void Elemental::LAPACK::Internal::TrinvU
( const Diagonal diagonal, DistMatrix<scomplex,MC,MR>& U );

template void Elemental::LAPACK::Internal::TrinvUVar3
( const Diagonal diagonal, DistMatrix<scomplex,MC,MR>& U );

template void Elemental::LAPACK::Internal::TrinvU
( const Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& U );

template void Elemental::LAPACK::Internal::TrinvUVar3
( const Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& U );
#endif

