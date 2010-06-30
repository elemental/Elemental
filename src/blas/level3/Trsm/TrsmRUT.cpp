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

// Right Upper (Conjugate)Transpose (Non)Unit Trsm
//   X := X triu(U)^-T, 
//   X := X triu(U)^-H,
//   X := X triuu(U)^-T, or
//   X := X triuu(U)^-H
template<typename T>
void
elemental::blas::internal::TrsmRUT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& U,
        DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrsmRUT");
    if( U.GetGrid() != X.GetGrid() )
        throw logic_error( "U and X must be distributed over the same grid." );
    if( orientation == Normal )
        throw logic_error( "TrsmRUT expects a (Conjugate)Transpose option." );
    if( U.Height() != U.Width() || X.Width() != U.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmRUT: " << endl
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

    DistMatrix<T,MC,MR> XL(g), XR(g),
                        X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,MR,  Star> U01_MR_Star(g);
    DistMatrix<T,Star,Star> U11_Star_Star(g);
    DistMatrix<T,MC,  Star> X1_MC_Star(g);
    DistMatrix<T,VC,  Star> X1_VC_Star(g);
    
    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionLeft( X, XL, XR, 0 );
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

        X1_MC_Star.AlignWith( X0 );
        U01_MR_Star.AlignWith( X0 );
        //--------------------------------------------------------------------//
        U11_Star_Star = U11; // U11[*,*] <- U11[MC,MR]
        X1_VC_Star    = X1;  // X1[VC,*] <- X1[MC,MR]

        // X1[VC,*] := X1[VC,*] (U11[*,*])^-(T/H)
        blas::internal::LocalTrsm
        ( Right, Upper, orientation, diagonal, 
          (T)1, U11_Star_Star, X1_VC_Star );

        X1_MC_Star  = X1_VC_Star; // X1[MC,*]  <- X1[VC,*]
        X1          = X1_MC_Star; // X1[MC,MR] <- X1[MC,*]
        U01_MR_Star = U01;        // U01[MR,*] <- U01[MC,MR]

        // X0[MC,MR] -= X1[MC,*] (U01[MR,*])^(T/H)
        //            = X1[MC,*] (U01^(T/H))[*,MR]
        blas::internal::LocalGemm
        ( Normal, orientation, (T)-1, X1_MC_Star, U01_MR_Star, (T)1, X0 );
        //--------------------------------------------------------------------//
        X1_MC_Star.FreeAlignments();
        U01_MR_Star.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12, 
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        SlidePartitionLeft
        ( XL, /**/     XR,
          X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrsmRUT
( Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& U,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrsmRUT
( Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& U,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrsmRUT
( Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& U,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrsmRUT
( Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& U,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

