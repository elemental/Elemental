/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::TrsvLN
( Diagonal diagonal, 
  const DistMatrix<T,MC,MR>& L, 
        DistMatrix<T,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrsvLN");
    if( L.GetGrid() != x.GetGrid() )
        throw "L and x must be distributed over the same grid.";
    if( L.Height() != L.Width() )
        throw "L must be square.";
    if( x.Width() != 1 && x.Height() != 1 )
        throw "x must be a vector.";
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( L.Width() != xLength )
        throw "Nonconformal TrsvLN.";
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
        DistMatrix<T,MR,  Star> x1_MR_Star(grid);
        DistMatrix<T,MC,  Star> z2_MC_Star(grid);

        // Start the algorithm
        LockedPartitionDownDiagonal
        ( L, LTL, LTR,
             LBL, LBR );
        PartitionDown
        ( x, xT,
             xB );
        while( xB.Height() > 0 )
        {
            LockedRepartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            RepartitionDown
            ( xT,  x0,
             /**/ /**/
                   x1,
              xB,  x2 );

            x1_MR_Star.AlignWith( L21 );
            z2_MC_Star.AlignWith( L21 );
            z2_MC_Star.ResizeTo( x2.Height(), 1 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            blas::Trsv
            ( Lower, Normal, diagonal,
              L11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_MR_Star = x1_Star_Star;
            blas::Gemv
            ( Normal, (T)-1, 
              L21.LockedLocalMatrix(), 
              x1_MR_Star.LockedLocalMatrix(),
              (T)0, z2_MC_Star.LocalMatrix() );
            x2.SumScatterUpdate( (T)1, z2_MC_Star );
            //----------------------------------------------------------------//
            x1_MR_Star.FreeAlignments();
            z2_MC_Star.FreeAlignments();

            SlideLockedPartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            SlidePartitionDown
            ( xT,  x0,
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
        LockedPartitionDownDiagonal
        ( L, LTL, LTR,
             LBL, LBR );
        PartitionRight( x,  xL, xR );
        while( xR.Width() > 0 )
        {
            LockedRepartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            RepartitionRight
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );

            x1_Star_MR.AlignWith( L21 );
            z2_Star_MC.AlignWith( L21 );
            z2.AlignWith( x2 );
            z2_Star_MC.ResizeTo( 1, x2.Width() );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            L11_Star_Star = L11;
            blas::Trsv
            ( Lower, Normal, diagonal,
              L11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_Star_MR = x1_Star_Star;
            blas::Gemv
            ( Normal, (T)-1, 
              L21.LockedLocalMatrix(), 
              x1_Star_MR.LockedLocalMatrix(),
              (T)0, z2_Star_MC.LocalMatrix() );
            z2_MR_MC.SumScatterFrom( z2_Star_MC );
            z2 = z2_MR_MC;
            blas::Axpy( (T)1, z2, x2 );
            //----------------------------------------------------------------//
            x1_Star_MR.FreeAlignments();
            z2_Star_MC.FreeAlignments();
            z2.FreeAlignments(); 

            SlideLockedPartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

            SlidePartitionRight
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrsvLN
( Diagonal diagonal,
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& x );

template void elemental::blas::internal::TrsvLN
( Diagonal diagonal,
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrsvLN
( Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& x );

template void elemental::blas::internal::TrsvLN
( Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

