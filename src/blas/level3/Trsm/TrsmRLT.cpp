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

// Right Lower (Conjugate)Transpose (Non)Unit Trsm
//   X := X tril(L)^-T,
//   X := X tril(L)^-H,
//   X := X trilu(L)^-T, or
//   X := X trilu(L)^-H
template<typename T>
void
elemental::blas::internal::TrsmRLT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& L,
        DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrsmRLT");
    if( L.GetGrid() != X.GetGrid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( orientation == Normal )
        throw logic_error( "TrsmRLT expects a (Conjugate)Transpose option." );
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmRLT: " << endl
            << "  L ~ " << L.Height() << " x " << L.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = L.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<T,MC,MR> XL(g), XR(g),
                        X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,MR,  Star> L21_MR_Star(g);
    DistMatrix<T,MC,  Star> X1_MC_Star(g);
    DistMatrix<T,VC,  Star> X1_VC_Star(g);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionRight
        ( XL, /**/     XR,
          X0, /**/ X1, X2 );

        X1_MC_Star.AlignWith( X2 );
        L21_MR_Star.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11; // L11[*,*] <- L11[MC,MR]
        X1_VC_Star    = X1;  // X1[VC,*] <- X1[MC,MR]
        
        // X1[VC,*] := X1[VC,*] (L11[*,*])^-(T/H)
        blas::internal::LocalTrsm
        ( Right, Lower, orientation, diagonal, 
          (T)1, L11_Star_Star, X1_VC_Star );

        X1_MC_Star  = X1_VC_Star; // X1[MC,*]  <- X1[VC,*]
        X1          = X1_MC_Star; // X1[MC,MR] <- X1[MC,*]
        L21_MR_Star = L21;        // L21[MR,*] <- L21[MC,MR]

        // X2[MC,MR] -= X1[MC,*] (L21[MR,*])^(T/H)
        //            = X1[MC,*] (L21^(T/H))[*,MR]
        blas::internal::LocalGemm
        ( Normal, orientation, (T)-1, X1_MC_Star, L21_MR_Star, (T)1, X2 );
        //--------------------------------------------------------------------//
        X1_MC_Star.FreeAlignments();
        L21_MR_Star.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrsmRLT
( Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrsmRLT
( Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrsmRLT
( Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrsmRLT
( Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

