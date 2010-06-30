/*
   This file is part of Elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (c) 2009-2010 Jack Poulson <jack.poulson@gmail.com>.
   All rights reserved.

   This file is released under the terms of the license contained in the file
   LICENSE-PURE.
*/
#include "elemental/blas_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::blas::internal::TrsvUN
( Diagonal diagonal, 
  const DistMatrix<T,MC,MR>& U, 
        DistMatrix<T,MC,MR>& x )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrsvUN");
    if( U.GetGrid() != x.GetGrid() )
        throw logic_error( "U and x must be distributed over the same grid." );
    if( U.Height() != U.Width() )
        throw logic_error( "U must be square." );
    if( x.Width() != 1 && x.Height() != 1 )
        throw logic_error( "x must be a vector." );
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( U.Width() != xLength )
        throw logic_error( "Nonconformal TrsvUN." );
#endif
    const Grid& g = U.GetGrid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<T,MC,MR> 
            UTL(g), UTR(g),  U00(g), U01(g), U02(g),
            UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                             U20(g), U21(g), U22(g);

        DistMatrix<T,MC,MR> 
            xT(g),  x0(g),
            xB(g),  x1(g),
                    x2(g);

        // Temporary distributions
        DistMatrix<T,Star,Star> U11_Star_Star(g);
        DistMatrix<T,Star,Star> x1_Star_Star(g);
        DistMatrix<T,MR,  Star> x1_MR_Star(g);
        DistMatrix<T,MC,  Star> z0_MC_Star(g);

        // Start the algorithm
        LockedPartitionUpDiagonal
        ( U, UTL, UTR,
             UBL, UBR, 0 );
        PartitionUp
        ( x, xT,
             xB, 0 );
        while( xT.Height() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, U01, /**/ U02,
                   /**/       U10, U11, /**/ U12,
             /*************/ /******************/
              UBL, /**/ UBR,  U20, U21, /**/ U22 );

            RepartitionUp
            ( xT,  x0,
                   x1,
             /**/ /**/
              xB,  x2 );

            x1_MR_Star.AlignWith( U01 );
            z0_MC_Star.AlignWith( U01 );
            z0_MC_Star.ResizeTo( x0.Height(), 1 );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            U11_Star_Star = U11;
            blas::Trsv
            ( Upper, Normal, diagonal,
              U11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_MR_Star = x1_Star_Star;
            blas::Gemv
            ( Normal, (T)-1, 
              U01.LockedLocalMatrix(), 
              x1_MR_Star.LockedLocalMatrix(),
              (T)0, z0_MC_Star.LocalMatrix() );
            x0.SumScatterUpdate( (T)1, z0_MC_Star );
            //----------------------------------------------------------------//
            x1_MR_Star.FreeAlignments();
            z0_MC_Star.FreeAlignments();

            SlideLockedPartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, /**/ U01, U02,
             /*************/ /******************/
                   /**/       U10, /**/ U11, U12,
              UBL, /**/ UBR,  U20, /**/ U21, U22 );

            SlidePartitionUp
            ( xT,  x0,
             /**/ /**/
                   x1,
              xB,  x2 );
        }
    }
    else
    {
        // Matrix views 
        DistMatrix<T,MC,MR> 
            UTL(g), UTR(g),  U00(g), U01(g), U02(g),
            UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                             U20(g), U21(g), U22(g);

        DistMatrix<T,MC,MR> 
            xL(g), xR(g),
            x0(g), x1(g), x2(g);

        // Temporary distributions
        DistMatrix<T,Star,Star> U11_Star_Star(g);
        DistMatrix<T,Star,Star> x1_Star_Star(g);
        DistMatrix<T,Star,MR  > x1_Star_MR(g);
        DistMatrix<T,Star,MC  > z0_Star_MC(g);
        DistMatrix<T,MR,  MC  > z0_MR_MC(g);
        DistMatrix<T,MC,  MR  > z0(g);

        // Start the algorithm
        LockedPartitionUpDiagonal
        ( U, UTL, UTR,
             UBL, UBR, 0 );
        PartitionLeft( x,  xL, xR, 0 );
        while( xL.Width() > 0 )
        {
            LockedRepartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, U01, /**/ U02,
                   /**/       U10, U11, /**/ U12,
             /*************/ /******************/
              UBL, /**/ UBR,  U20, U21, /**/ U22 );

            RepartitionLeft
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );

            x1_Star_MR.AlignWith( U01 );
            z0_Star_MC.AlignWith( U01 );
            z0.AlignWith( x0 );
            z0_Star_MC.ResizeTo( 1, x0.Width() );
            //----------------------------------------------------------------//
            x1_Star_Star = x1;
            U11_Star_Star = U11;
            blas::Trsv
            ( Upper, Normal, diagonal,
              U11_Star_Star.LockedLocalMatrix(),
              x1_Star_Star.LocalMatrix() );
            x1 = x1_Star_Star;

            x1_Star_MR = x1_Star_Star;
            blas::Gemv
            ( Normal, (T)-1, 
              U01.LockedLocalMatrix(), 
              x1_Star_MR.LockedLocalMatrix(),
              (T)0, z0_Star_MC.LocalMatrix() );
            z0_MR_MC.SumScatterFrom( z0_Star_MC );
            z0 = z0_MR_MC;
            blas::Axpy( (T)1, z0, x0 );
            //----------------------------------------------------------------//
            x1_Star_MR.FreeAlignments();
            z0_Star_MC.FreeAlignments();
            z0.FreeAlignments(); 

            SlideLockedPartitionUpDiagonal
            ( UTL, /**/ UTR,  U00, /**/ U01, U02,
             /*************/ /******************/
                   /**/       U10, /**/ U11, U12,
              UBL, /**/ UBR,  U20, /**/ U21, U22 );

            SlidePartitionLeft
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrsvUN
( Diagonal diagonal,
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& x );

template void elemental::blas::internal::TrsvUN
( Diagonal diagonal,
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& x );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrsvUN
( Diagonal diagonal,
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& x );

template void elemental::blas::internal::TrsvUN
( Diagonal diagonal,
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& x );
#endif

