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

// Right Upper (Conjugate)Transpose (Non)Unit Trmm
//   X := X triu(U)^T, 
//   X := X triu(U)^H,
//   X := X triuu(U)^T, or
//   X := X triuu(U)^H
template<typename T>
void
elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRUT");
    if( U.GetGrid() != X.GetGrid() )
        throw "U and X must be distributed over the same grid.";
    if( orientation == Normal )
        throw "TrmmRUT expects a (Conjugate)Transpose option.";
    if( U.Height() != U.Width() || X.Width() != U.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmRUT: " << endl
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
    DistMatrix<T,Star,Star> U11_Star_Star(grid);
    DistMatrix<T,Star,MR  > U12_Star_MR(grid);
    DistMatrix<T,VC,  Star> X1_VC_Star(grid);
    DistMatrix<T,MC,  Star> D1_MC_Star(grid);
    
    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR );
    PartitionRight( X, XL, XR );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        U12_Star_MR.AlignWith( X2 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        U11_Star_Star = U11;
        blas::internal::LocalTrmm
        ( Right, Upper, orientation, diagonal, 
          (T)1, U11_Star_Star, X1_VC_Star );
        X1 = X1_VC_Star;
 
        U12_Star_MR = U12;
        blas::internal::LocalGemm
        ( Normal, orientation, (T)1, X2, U12_Star_MR, (T)0, D1_MC_Star );
        X1.SumScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        U12_Star_MR.FreeAlignments();
        D1_MC_Star.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, U01, /**/ U02,
               /**/       U10, U11, /**/ U12, 
         /*************/ /******************/
          UBL, /**/ UBR,  U20, U21, /**/ U22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmRUT
( Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

