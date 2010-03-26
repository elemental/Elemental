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
Elemental::BLAS::Internal::TrsvLN
( const Diagonal diagonal, 
  const DistMatrix<T,MC,MR>& L, 
        DistMatrix<T,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("BLAS::Internal::TrsvLN");
#endif
    const Grid& grid = L.GetGrid();
#ifndef RELEASE
    if( L.GetGrid() != x.GetGrid() )
    {
        if( grid.VCRank() == 0 )
            cerr << "L and x must be distributed over the same grid." << endl;
        DumpCallStack();
        throw exception();
    }
    if( L.Height() != L.Width() )
    {
        if( grid.VCRank() == 0 )
            cerr << "L must be square." << endl;
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
        if( L.Width() != xLength )
        {
            if( grid.VCRank() == 0 )
                cerr << "Nonconformal TrsvLN." << endl;
            DumpCallStack();
            throw exception();
        }
    }
#endif
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
        DistMatrix<T,MR,  Star> x1_MR_Star(grid);
        DistMatrix<T,MC,  Star> z2_MC_Star(grid);

        // Start the algorithm
        LockedPartitionDownDiagonal( L, LTL, LTR,
                                        LBL, LBR );
        PartitionDown( x, xT,
                          xB );
        while( xB.Height() > 0 )
        {
            LockedRepartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            RepartitionDown( xT,  x0,
                            /**/ /**/
                                  x1,
                             xB,  x2 );

            x1_MR_Star.ConformWith( L21 );
            z2_MC_Star.AlignWith( L21 );
            z2_MC_Star.ResizeTo( x2.Height(), 1 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            BLAS::Trsv( Lower, Normal, diagonal,
                        L11_Star_Star.LockedLocalMatrix(),
                        x1_Star_Star.LocalMatrix()        );
            x1 = x1_Star_Star;

            x1_MR_Star = x1_Star_Star;
            BLAS::Gemv( Normal, (T)-1, 
                        L21.LockedLocalMatrix(), 
                        x1_MR_Star.LockedLocalMatrix(),
                        (T)0, z2_MC_Star.LocalMatrix() );
            x2.ReduceScatterUpdate( (T)1, z2_MC_Star );
            //----------------------------------------------------------------//
            x1_MR_Star.FreeConstraints();
            z2_MC_Star.FreeConstraints();

            SlideLockedPartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            SlidePartitionDown( xT,  x0,
                                     x1,
                               /**/ /**/
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
        DistMatrix<T,Star,MR  > x1_Star_MR(grid);
        DistMatrix<T,Star,MC  > z2_Star_MC(grid);
        DistMatrix<T,MR,  MC  > z2_MR_MC(grid);
        DistMatrix<T,MC,  MR  > z2(grid);

        // Start the algorithm
        LockedPartitionDownDiagonal( L, LTL, LTR,
                                        LBL, LBR );
        PartitionRight( x,  xL, xR );
        while( xR.Width() > 0 )
        {
            LockedRepartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            RepartitionRight( xL, /**/ xR,
                              x0, /**/ x1, x2 );

            x1_Star_MR.ConformWith( L21 );
            z2_Star_MC.AlignWith( L21 );
            z2.AlignWith( x2 );
            z2_Star_MC.ResizeTo( 1, x2.Width() );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            BLAS::Trsv( Lower, Normal, diagonal,
                        L11_Star_Star.LockedLocalMatrix(),
                        x1_Star_Star.LocalMatrix()        );
            x1 = x1_Star_Star;

            x1_Star_MR = x1_Star_Star;
            BLAS::Gemv( Normal, (T)-1, 
                        L21.LockedLocalMatrix(), 
                        x1_Star_MR.LockedLocalMatrix(),
                        (T)0, z2_Star_MC.LocalMatrix() );
            z2_MR_MC.ReduceScatterFrom( z2_Star_MC );
            z2 = z2_MR_MC;
            BLAS::Axpy( (T)1, z2, x2 );
            //----------------------------------------------------------------//
            x1_Star_MR.FreeConstraints();
            z2_Star_MC.FreeConstraints();
            z2.FreeConstraints(); 

            SlideLockedPartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            SlidePartitionRight( xL,     /**/ xR,
                                 x0, x1, /**/ x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::BLAS::Internal::TrsvLN
( const Diagonal diagonal,
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& x );

template void Elemental::BLAS::Internal::TrsvLN
( const Diagonal diagonal,
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void Elemental::BLAS::Internal::TrsvLN
( const Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& x );

template void Elemental::BLAS::Internal::TrsvLN
( const Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

