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

// Left Upper (Conjugate)Transpose (Non)Unit Trmm
//   X := triu(U)^T  X, 
//   X := triu(U)^H  X,
//   X := triuu(U)^T X, or
//   X := triuu(U)^H X
template<typename T>
void
elemental::blas::internal::TrmmLUT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmLUT");
    if( U.GetGrid() != X.GetGrid() )
        throw "U and X must be distributed over the same grid.";
    if( orientation == Normal )
        throw "TrmmLUT expects a (Conjugate)Transpose option.";
    if( U.Height() != U.Width() || U.Height() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLUT: " << endl
            << "  U ~ " << U.Height() << " x " << U.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
    }
#endif
    const Grid& grid = U.GetGrid();

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
    DistMatrix<T,Star,VR  > X1_Star_VR(grid);
    DistMatrix<T,Star,MR  > D1_Star_MR(grid);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR );
    PartitionUp
    ( X, XT,
         XB );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12,
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        U01_MC_Star.AlignWith( X0 );
        D1_Star_MR.AlignWith( X1 );
        D1_Star_MR.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_Star_VR = X1;
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Left, Upper, orientation, diagonal, (T)1, U11_Star_Star, X1_Star_VR );
        X1 = X1_Star_VR;
        
        U01_MC_Star = U01;
        blas::internal::LocalGemm
        ( orientation, Normal, (T)1, U01_MC_Star, X0, (T)0, D1_Star_MR );
        X1.SumScatterUpdate( (T)1, D1_Star_MR );
        //--------------------------------------------------------------------//
        U01_MC_Star.FreeAlignments();
        D1_Star_MR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( UTL, /**/ UTR,   U00, /**/ U01, U02,
         /*************/  /******************/
               /**/        U10, /**/ U11, U12,
          UBL, /**/ UBR,   U20, /**/ U21, U22 );

        SlidePartitionUp
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrmmLUT
( Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmLUT
( Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmLUT
( Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmLUT
( Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

