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
#include "ElementalLAPACK_Internal.h"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::LAPACK::Internal::TrinvL
( const Diagonal diagonal, DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::TrinvL");
#endif
    TrinvL_Var3( diagonal, L );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
Elemental::LAPACK::Internal::TrinvL_Var3
( const Diagonal diagonal, DistMatrix<T,MC,MR>& L )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::TrinvL_Var3");
    if( L.Height() != L.Width() )
    {
        if( L.GetGrid().VCRank() == 0 )
            cerr << "Nonsquare matrices cannot be triangular." << endl;
        DumpCallStack();
        throw exception();
    }
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
        RepartitionDownDiagonal( LTL, /**/ LTR,  L00, /**/ L01, L02,
                                /*************/ /******************/
                                      /**/       L10, /**/ L11, L12,
                                 LBL, /**/ LBR,  L20, /**/ L21, L22 );

        L10_Star_MR.AlignWith( L20 );
        L21_MC_Star.AlignWith( L20 );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11;
        LAPACK::Trinv( Lower, diagonal, L11_Star_Star.LocalMatrix() );
        L11 = L11_Star_Star;

        L10_Star_VR = L10;
        BLAS::Trmm( Left, Lower, Normal, diagonal,
                    (T)-1, L11_Star_Star.LockedLocalMatrix(),
                           L10_Star_VR.LocalMatrix()         );

        L21_MC_Star = L21;
        L10_Star_MR = L10_Star_VR;
        BLAS::Gemm( Normal, Normal,
                    (T)1, L21_MC_Star.LockedLocalMatrix(),
                          L10_Star_MR.LockedLocalMatrix(),
                    (T)1, L20.LocalMatrix()               );
        L10 = L10_Star_MR;

        L21_VC_Star = L21_MC_Star;
        BLAS::Trmm( Right, Lower, Normal, diagonal,
                    (T)1, L11_Star_Star.LockedLocalMatrix(),
                          L21_VC_Star.LocalMatrix()         );
        L21 = L21_VC_Star;
        //--------------------------------------------------------------------//
        L10_Star_MR.FreeConstraints();
        L21_MC_Star.FreeConstraints();

        SlidePartitionDownDiagonal( LTL, /**/ LTR,  L00, L01, /**/ L02,
                                         /**/       L10, L11, /**/ L12,
                                   /*************/ /******************/
                                    LBL, /**/ LBR,  L20, L21, /**/ L22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::LAPACK::Internal::TrinvL
( const Diagonal diagonal, DistMatrix<float,MC,MR>& L );

template void Elemental::LAPACK::Internal::TrinvL_Var3
( const Diagonal diagonal, DistMatrix<float,MC,MR>& L );

template void Elemental::LAPACK::Internal::TrinvL
( const Diagonal diagonal, DistMatrix<double,MC,MR>& L );

template void Elemental::LAPACK::Internal::TrinvL_Var3
( const Diagonal diagonal, DistMatrix<double,MC,MR>& L );

#ifndef WITHOUT_COMPLEX
template void Elemental::LAPACK::Internal::TrinvL
( const Diagonal diagonal, DistMatrix<scomplex,MC,MR>& L );

template void Elemental::LAPACK::Internal::TrinvL_Var3
( const Diagonal diagonal, DistMatrix<scomplex,MC,MR>& L );

template void Elemental::LAPACK::Internal::TrinvL
( const Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& L );

template void Elemental::LAPACK::Internal::TrinvL_Var3
( const Diagonal diagonal, DistMatrix<dcomplex,MC,MR>& L );
#endif

