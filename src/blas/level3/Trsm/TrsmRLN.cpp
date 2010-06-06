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

// Right Lower Normal (Non)Unit Trsm
//   X := X tril(L)^-1, and
//   X := X trilu(L)^-1
template<typename T>
void
elemental::blas::internal::TrsmRLN
( Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& L,
        DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrsmRLN");
    if( L.GetGrid() != X.GetGrid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmRLN: " << endl
            << "  L ~ " << L.Height() << " x " << L.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& grid = L.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        LTL(grid), LTR(grid),  L00(grid), L01(grid), L02(grid),
        LBL(grid), LBR(grid),  L10(grid), L11(grid), L12(grid),
                               L20(grid), L21(grid), L22(grid);

    DistMatrix<T,MC,MR> XL(grid), XR(grid),
                        X0(grid), X1(grid), X2(grid);

    // Temporary distributions
    DistMatrix<T,Star,MR  > L10_Star_MR(grid);
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,MC,  Star> X1_MC_Star(grid);
    DistMatrix<T,VC,  Star> X1_VC_Star(grid);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionLeft( X, XL, XR, 0 );
    while( XL.Width() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );
 
        RepartitionLeft
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );

        X1_MC_Star.AlignWith( X0 );
        L10_Star_MR.AlignWith( X0 );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11; // L11[*,*] <- L11[MC,MR]
        X1_VC_Star    = X1;  // X1[VC,*] <- X1[MC,MR]

        // X1[VC,*] := X1[VC,*] (L11[*,*])^-1
        blas::internal::LocalTrsm
        ( Right, Lower, Normal, diagonal, (T)1, L11_Star_Star, X1_VC_Star );

        X1_MC_Star  = X1_VC_Star; // X1[MC,*]  <- X1[VC,*]
        X1          = X1_MC_Star; // X1[MC,MR] <- X1[MC,*]
        L10_Star_MR = L10;        // L10[*,MR] <- L10[MC,MR]

        // X0[MC,MR] -= X1[MC,*] L10[*,MR]
        blas::internal::LocalGemm
        ( Normal, Normal, (T)-1, X1_MC_Star, L10_Star_MR, (T)1, X0 );
        //--------------------------------------------------------------------//
        X1_MC_Star.FreeAlignments();
        L10_Star_MR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        SlidePartitionLeft
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrsmRLN
( Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrsmRLN
( Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrsmRLN
( Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrsmRLN
( Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

