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

// Right Upper Normal (Non)Unit Trmm
//   X := X triu(U), and
//   X := X triuu(U)
template<typename T>
void
elemental::blas::internal::TrmmRUN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRUN");
    if( U.GetGrid() != X.GetGrid() )
        throw "U and X must be distributed over the same grid.";
    if( U.Height() != U.Width() || X.Width() != U.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmRUN: " << endl
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

    DistMatrix<T,MC,MR> XL(grid), XR(grid),
                        X0(grid), X1(grid), X2(grid);

    // Temporary distributions
    DistMatrix<T,MR,  Star> U01_MR_Star(grid);
    DistMatrix<T,Star,Star> U11_Star_Star(grid); 
    DistMatrix<T,VC,  Star> X1_VC_Star(grid);    
    DistMatrix<T,MC,  Star> D1_MC_Star(grid);
    
    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR );
    PartitionLeft( X, XL, XR );
    while( XL.Width() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12, 
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        RepartitionLeft
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 ); 

        U01_MR_Star.AlignWith( X0 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Right, Upper, Normal, diagonal, (T)1, U11_Star_Star, X1_VC_Star );
        X1 = X1_VC_Star;
 
        U01_MR_Star = U01;
        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, X0, U01_MR_Star, (T)0, D1_MC_Star );
        X1.SumScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        U01_MR_Star.FreeAlignments();
        D1_MC_Star.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        SlidePartitionLeft
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrmmRUN
( Diagonal diagonal,
  float alpha, const DistMatrix<float,MC,MR>& U,
                     DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmRUN
( Diagonal diagonal,
  double alpha, const DistMatrix<double,MC,MR>& U,
                      DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmRUN
( Diagonal diagonal,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& U,
                        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmRUN
( Diagonal diagonal,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& U,
                        DistMatrix<dcomplex,MC,MR>& X );
#endif

