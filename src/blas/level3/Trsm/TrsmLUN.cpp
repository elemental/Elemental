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

// Left Upper Normal (Non)Unit Trsm
//   X := triu(U)^-1  X, or
//   X := triuu(U)^-1 X
template<typename T>
void
blas::internal::TrsmLUN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& U,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrsmLUN");
    if( U.GetGrid() != X.GetGrid() )
        throw logic_error( "U and X must be distributed over the same grid." );
    if( U.Height() != U.Width() || U.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmLUN: " << endl
            << "  U ~ " << U.Height() << " x " << U.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = U.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<T,MC,MR> XT(g),  X0(g),
                        XB(g),  X1(g),
                                X2(g);

    // Temporary distributions
    DistMatrix<T,MC,  Star> U01_MC_Star(g);
    DistMatrix<T,Star,Star> U11_Star_Star(g);
    DistMatrix<T,Star,MR  > X1_Star_MR(g);
    DistMatrix<T,Star,VR  > X1_Star_VR(g);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionUp
    ( X, XT,
         XB, 0 );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        U01_MC_Star.AlignWith( X0 );
        X1_Star_MR.AlignWith( X0 );
        //--------------------------------------------------------------------//
        U11_Star_Star = U11; // U11[*,*] <- U11[MC,MR]
        X1_Star_VR    = X1;  // X1[*,VR] <- X1[MC,MR]
        
        // X1[*,VR] := (U11[*,*])^-1 X1[*,VR]
        blas::internal::LocalTrsm
        ( Left, Upper, Normal, diagonal, (T)1, U11_Star_Star, X1_Star_VR );

        X1_Star_MR  = X1_Star_VR; // X1[*,MR]  <- X1[*,VR]
        X1          = X1_Star_MR; // X1[MC,MR] <- X1[*,MR]
        U01_MC_Star = U01;        // U01[MC,*] <- U01[MC,MR]

        // X0[MC,MR] -= U01[MC,*] X1[*,MR]
        blas::internal::LocalGemm
        ( Normal, Normal, (T)-1, U01_MC_Star, X1_Star_MR, (T)1, X0 );
        //--------------------------------------------------------------------//
        U01_MC_Star.FreeAlignments();
        X1_Star_MR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

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

template void elemental::blas::internal::TrsmLUN
( Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrsmLUN
( Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrsmLUN
( Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrsmLUN
( Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

