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
#include "ElementalBLASInternal.h"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Internal::TrsvUN
( const Diagonal diagonal, 
  const DistMatrix<T,MC,MR>& U, 
        DistMatrix<T,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrsvUN");
#endif
    const Grid& grid = U.GetGrid();
#ifndef RELEASE
    if( U.GetGrid() != x.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "U and x must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( U.Height() != U.Width() )
    {
        if( grid.VCRank() == 0 )
            cerr << "U must be square." << endl;
        DumpCallStack();
        throw exception();
    }
    if( x.Width() != 1 && x.Height() != 1 )
    {
        if( grid.VCRank() == 0 )
            cerr << "x must be a vector." << endl;
        DumpCallStack();
        throw exception();
    }
    {
        const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( U.Width() != xLength )
        {
            if( grid.VCRank() == 0 )
                cerr << "Nonconformal TrsvUN." << endl;
            DumpCallStack();
            throw exception();
        }
    }
#endif
    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<T,MC,MR> 
            UTL(grid), UTR(grid),  U00(grid), U01(grid), U02(grid),
            UBL(grid), UBR(grid),  U10(grid), U11(grid), U12(grid),
                                   U20(grid), U21(grid), U22(grid);

        DistMatrix<T,MC,MR> xT(grid),  x0(grid),
                            xB(grid),  x1(grid),
                                       x2(grid);

        // Temporary distributions
        DistMatrix<T,Star,Star> U11_Star_Star(grid);
        DistMatrix<T,Star,Star> x1_Star_Star(grid);
        DistMatrix<T,MR,  Star> x1_MR_Star(grid);
        DistMatrix<T,MC,  Star> z0_MC_Star(grid);

        // Start the algorithm
        LockedPartitionUpDiagonal( U, UTL, UTR,
                                      UBL, UBR );
        PartitionUp( x, xT,
                        xB );
        while( xT.Height() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, U01, /**/ U02,
                   /**/       U10, U11, /**/ U12,
             /*************/ /******************/
              UBL, /**/ UBR,  U20, U21, /**/ U22 );

            RepartitionUp( xT,  x0,
                                x1,
                          /**/ /**/
                           xB,  x2 );

            x1_MR_Star.ConformWith( U01 );
            z0_MC_Star.AlignWith( U01 );
            z0_MC_Star.ResizeTo( x0.Height(), 1 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            U11_Star_Star = U11;
            BLAS::Trsv( Upper, Normal, diagonal,
                        U11_Star_Star.LockedLocalMatrix(),
                        x1_Star_Star.LocalMatrix()        );
            x1 = x1_Star_Star;

            x1_MR_Star = x1_Star_Star;
            BLAS::Gemv( Normal, (T)-1, 
                        U01.LockedLocalMatrix(), 
                        x1_MR_Star.LockedLocalMatrix(),
                        (T)0, z0_MC_Star.LocalMatrix() );
            x0.ReduceScatterUpdate( (T)1, z0_MC_Star );
            //----------------------------------------------------------------//
            x1_MR_Star.FreeConstraints();
            z0_MC_Star.FreeConstraints();

            SlideLockedPartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, /**/ U01, U02,
             /*************/ /******************/
                   /**/       U10, /**/ U11, U12,
              UBL, /**/ UBR,  U20, /**/ U21, U22 );

            SlidePartitionUp( xT,  x0,
                             /**/ /**/
                                   x1,
                              xB,  x2 );
        }
    }
    else
    {
        // Matrix views 
        DistMatrix<T,MC,MR> 
            UTL(grid), UTR(grid),  U00(grid), U01(grid), U02(grid),
            UBL(grid), UBR(grid),  U10(grid), U11(grid), U12(grid),
                                   U20(grid), U21(grid), U22(grid);

        DistMatrix<T,MC,MR> xL(grid), xR(grid),
                            x0(grid), x1(grid), x2(grid);

        // Temporary distributions
        DistMatrix<T,Star,Star> U11_Star_Star(grid);
        DistMatrix<T,Star,Star> x1_Star_Star(grid);
        DistMatrix<T,Star,MR  > x1_Star_MR(grid);
        DistMatrix<T,Star,MC  > z0_Star_MC(grid);
        DistMatrix<T,MR,  MC  > z0_MR_MC(grid);
        DistMatrix<T,MC,  MR  > z0(grid);

        // Start the algorithm
        LockedPartitionUpDiagonal( U, UTL, UTR,
                                      UBL, UBR );
        PartitionLeft( x,  xL, xR );
        while( xL.Width() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, U01, /**/ U02,
                   /**/       U10, U11, /**/ U12,
             /*************/ /******************/
              UBL, /**/ UBR,  U20, U21, /**/ U22 );

            RepartitionLeft( xL,     /**/ xR,
                             x0, x1, /**/ x2 );

            x1_Star_MR.ConformWith( U01 );
            z0_Star_MC.AlignWith( U01 );
            z0.AlignWith( x0 );
            z0_Star_MC.ResizeTo( 1, x0.Width() );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            U11_Star_Star = U11;
            BLAS::Trsv( Upper, Normal, diagonal,
                        U11_Star_Star.LockedLocalMatrix(),
                        x1_Star_Star.LocalMatrix()        );
            x1 = x1_Star_Star;

            x1_Star_MR = x1_Star_Star;
            BLAS::Gemv( Normal, (T)-1, 
                        U01.LockedLocalMatrix(), 
                        x1_Star_MR.LockedLocalMatrix(),
                        (T)0, z0_Star_MC.LocalMatrix() );
            z0_MR_MC.ReduceScatterFrom( z0_Star_MC );
            z0 = z0_MR_MC;
            BLAS::Axpy( (T)1, z0, x0 );
            //----------------------------------------------------------------//
            x1_Star_MR.FreeConstraints();
            z0_Star_MC.FreeConstraints();
            z0.FreeConstraints(); 

            SlideLockedPartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, /**/ U01, U02,
             /*************/ /******************/
                   /**/       U10, /**/ U11, U12,
              UBL, /**/ UBR,  U20, /**/ U21, U22 );

            SlidePartitionLeft( xL, /**/ xR,
                                x0, /**/ x1, x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::TrsvUN
( const Diagonal diagonal,
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& x );

template void Elemental::BLAS::Internal::TrsvUN
( const Diagonal diagonal,
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrsvUN
( const Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& x );

template void Elemental::BLAS::Internal::TrsvUN
( const Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

