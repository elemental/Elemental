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
#include "ElementalBLAS_Internal.h"
using namespace std;
using namespace Elemental;

// Left Upper Normal (Non)Unit Trsm
//   X := triu(U)^-1  X, or
//   X := triuu(U)^-1 X
template<typename T>
void
BLAS::Internal::TrsmLUN
( const Diagonal diagonal,
  const T alpha, const DistMatrix<T,MC,MR>& U,
                       DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrsmLUN");
#endif
    const Grid& grid = U.GetGrid();
#ifndef RELEASE
    if( U.GetGrid() != X.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "U and X must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( U.Height() != U.Width() || U.Width() != X.Height() )
    {
        if( grid.VCRank() == 0 )
        {
            cerr << "Nonconformal TrsmLUN: " <<
            endl << "  U ~ " << U.Height() << " x " << U.Width() <<
            endl << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        }
        DumpCallStack();
        throw exception();
    }
#endif
    // Matrix views
    DistMatrix<T,MC,MR> 
        UTL(grid), UTR(grid),  U00(grid), U01(grid), U02(grid),
        UBL(grid), UBR(grid),  U10(grid), U11(grid), U12(grid),
                               U20(grid), U21(grid), U22(grid);

    DistMatrix<T,MC,MR> XT(grid),  X0(grid),
                        XB(grid),  X1(grid),
                                   X2(grid);

    // Temporary distributions
    DistMatrix<T,MC,  Star> U01_MC_Star(grid);
    DistMatrix<T,Star,Star> U11_Star_Star(grid);
    DistMatrix<T,Star,MR  > X1_Star_MR(grid);
    DistMatrix<T,Star,VR  > X1_Star_VR(grid);

    // Start the algorithm
    BLAS::Scal( alpha, X );
    LockedPartitionUpDiagonal( U, UTL, UTR,
                                  UBL, UBR );
    PartitionUp( X, XT,
                    XB );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal( UTL, /**/ UTR,   U00, U01, /**/ U02,
                                          /**/        U10, U11, /**/ U12,
                                    /*************/  /******************/
                                     UBL, /**/ UBR,   U20, U21, /**/ U22 );

        RepartitionUp( XT,  X0,
                            X1,
                      /**/ /**/
                       XB,  X2 );

        U01_MC_Star.AlignWith( X0 );
        X1_Star_MR.AlignWith( X0 );
        //--------------------------------------------------------------------//
        U11_Star_Star = U11; // U11[*,*] <- U11[MC,MR]
        X1_Star_VR    = X1;  // X1[*,VR] <- X1[MC,MR]
        
        // X1[*,VR] := (U11[*,*])^-1 X1[*,VR]
        BLAS::Trsm( Left, Upper, Normal, diagonal,
                    (T)1, U11_Star_Star.LockedLocalMatrix(),
                          X1_Star_VR.LocalMatrix()          );

        X1_Star_MR  = X1_Star_VR; // X1[*,MR]  <- X1[*,VR]
        X1          = X1_Star_MR; // X1[MC,MR] <- X1[*,MR]
        U01_MC_Star = U01;        // U01[MC,*] <- U01[MC,MR]

        // X0[MC,MR] -= U01[MC,*] X1[*,MR]
        BLAS::Gemm( Normal, Normal, 
                    (T)-1, U01_MC_Star.LockedLocalMatrix(),
                           X1_Star_MR.LockedLocalMatrix(),
                    (T) 1, X0.LocalMatrix()                );
        //--------------------------------------------------------------------//
        U01_MC_Star.FreeConstraints();
        X1_Star_MR.FreeConstraints();

        SlideLockedPartitionUpDiagonal( UTL, /**/ UTR,  U00, /**/ U01, U02,
                                       /*************/ /******************/
                                             /**/       U10, /**/ U11, U12,
                                        UBL, /**/ UBR,  U20, /**/ U21, U22 );

        SlidePartitionUp( XT,  X0,
                         /**/ /**/
                               X1,
                          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::TrsmLUN
( const Diagonal diagonal,
  const float alpha, 
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& X );

template void Elemental::BLAS::Internal::TrsmLUN
( const Diagonal diagonal,
  const double alpha, 
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrsmLUN
( const Diagonal diagonal,
  const scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& X );

template void Elemental::BLAS::Internal::TrsmLUN
( const Diagonal diagonal,
  const dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

