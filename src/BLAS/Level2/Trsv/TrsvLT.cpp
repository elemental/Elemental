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
#include "Elemental/BLASInternal.hpp"
using namespace std;
using namespace Elemental;

template<typename T>
void
Elemental::BLAS::Internal::TrsvLT
( const Orientation orientation,
  const Diagonal diagonal, 
  const DistMatrix<T,MC,MR>& L, 
        DistMatrix<T,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrsvLT");
    if( L.GetGrid() != x.GetGrid() )
        throw "L and x must be distributed over the same grid.";
    if( orientation == Normal )
        throw "TrsvLT expects a (conjugate-)transpose option.";
    if( L.Height() != L.Width() )
        throw "L must be square.";
    if( x.Width() != 1 && x.Height() != 1 )
        throw "x must be a vector.";
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( L.Width() != xLength )
        throw "Nonconformal TrsvLT.";
#endif
    const Grid& grid = L.GetGrid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<T,MC,MR> 
            LTL(grid), LTR(grid),  L00(grid), L01(grid), L02(grid),
            LBL(grid), LBR(grid),  L10(grid), L11(grid), L12(grid),
                                   L20(grid), L21(grid), L22(grid);

        DistMatrix<T,MC,MR> xT(grid),  x0(grid),
                            xB(grid),  x1(grid),
                                       x2(grid);

        // Temporary distributions
        DistMatrix<T,Star,Star> L11_Star_Star(grid);
        DistMatrix<T,Star,Star> x1_Star_Star(grid);
        DistMatrix<T,MC,  Star> x1_MC_Star(grid);
        DistMatrix<T,MR,  Star> z0_MR_Star(grid);
        DistMatrix<T,MR,  MC  > z0_MR_MC(grid);
        DistMatrix<T,MC,  MR  > z0(grid);

        // Start the algorithm
        LockedPartitionUpDiagonal( L, LTL, LTR,
                                      LBL, LBR );
        PartitionUp( x, xT,
                        xB );
        while( xT.Height() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            RepartitionUp( xT,  x0,
                                x1,
                          /**/ /**/
                           xB,  x2 );

            x1_MC_Star.ConformWith( L10 );
            z0_MR_Star.AlignWith( L10 );
            z0_MR_Star.ResizeTo( x0.Height(), 1 );
            z0.AlignWith( x0 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            BLAS::Trsv( Lower, orientation, diagonal,
                        L11_Star_Star.LockedLocalMatrix(),
                        x1_Star_Star.LocalMatrix()        );
            x1 = x1_Star_Star;

            x1_MC_Star = x1_Star_Star;
            BLAS::Gemv( orientation, (T)-1, 
                        L10.LockedLocalMatrix(), 
                        x1_MC_Star.LockedLocalMatrix(),
                        (T)0, z0_MR_Star.LocalMatrix() );
            z0_MR_MC.ReduceScatterFrom( z0_MR_Star );
            z0 = z0_MR_MC;
            BLAS::Axpy( (T)1, z0, x0 );
            //----------------------------------------------------------------//
            x1_MC_Star.FreeConstraints();
            z0_MR_Star.FreeConstraints();
            z0.FreeConstraints();

            SlideLockedPartitionUpDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

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
            LTL(grid), LTR(grid),  L00(grid), L01(grid), L02(grid),
            LBL(grid), LBR(grid),  L10(grid), L11(grid), L12(grid),
                                   L20(grid), L21(grid), L22(grid);

        DistMatrix<T,MC,MR> xL(grid), xR(grid),
                            x0(grid), x1(grid), x2(grid);

        // Temporary distributions
        DistMatrix<T,Star,Star> L11_Star_Star(grid);
        DistMatrix<T,Star,Star> x1_Star_Star(grid);
        DistMatrix<T,Star,MC  > x1_Star_MC(grid);
        DistMatrix<T,Star,MR  > z0_Star_MR(grid);

        // Start the algorithm
        LockedPartitionUpDiagonal( L, LTL, LTR,
                                      LBL, LBR );
        PartitionLeft( x,  xL, xR );
        while( xL.Width() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            RepartitionLeft( xL,     /**/ xR,
                             x0, x1, /**/ x2 );

            x1_Star_MC.ConformWith( L10 );
            z0_Star_MR.AlignWith( L10 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            BLAS::Trsv( Lower, orientation, diagonal,
                        L11_Star_Star.LockedLocalMatrix(),
                        x1_Star_Star.LocalMatrix()        );
            x1 = x1_Star_Star;

            x1_Star_MC = x1_Star_Star;
            BLAS::Gemv( orientation, (T)-1, 
                        L10.LockedLocalMatrix(), 
                        x1_Star_MC.LockedLocalMatrix(),
                        (T)0, z0_Star_MR.LocalMatrix() );
            x0.ReduceScatterUpdate( (T)1, z0_Star_MR );
            //----------------------------------------------------------------//
            x1_Star_MC.FreeConstraints();
            z0_Star_MR.FreeConstraints();

            SlideLockedPartitionUpDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            SlidePartitionLeft( xL, /**/ xR,
                                x0, /**/ x1, x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::TrsvLT
( const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& x );

template void Elemental::BLAS::Internal::TrsvLT
( const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrsvLT
( const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& x );

template void Elemental::BLAS::Internal::TrsvLT
( const Orientation orientation,
  const Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

