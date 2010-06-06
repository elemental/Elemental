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

// Left Lower Normal (Non)Unit Trsm 
//   X := tril(L)^-1  X, or
//   X := trilu(L)^-1 X
template<typename T>
void
elemental::blas::internal::TrsmLLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrsmLLN");
    if( L.GetGrid() != X.GetGrid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmLLN: " << endl
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

    DistMatrix<T,MC,MR> XT(grid),  X0(grid), 
                        XB(grid),  X1(grid),
                                   X2(grid);

    // Temporary distributions
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,MC,  Star> L21_MC_Star(grid);
    DistMatrix<T,Star,MR  > X1_Star_MR(grid);
    DistMatrix<T,Star,VR  > X1_Star_VR(grid);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        L21_MC_Star.AlignWith( X2 );
        X1_Star_MR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11; // L11[*,*] <- L11[MC,MR]
        X1_Star_VR    = X1;  // X1[*,VR] <- X1[MC,MR]

        // X1[*,VR] := (L11[*,*])^-1 X1[*,VR]
        blas::internal::LocalTrsm
        ( Left, Lower, Normal, diagonal, (T)1, L11_Star_Star, X1_Star_VR );

        X1_Star_MR  = X1_Star_VR; // X1[*,MR]  <- X1[*,VR]
        X1          = X1_Star_MR; // X1[MC,MR] <- X1[*,MR]
        L21_MC_Star = L21;        // L21[MC,*] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,*] X1[*,MR]
        blas::internal::LocalGemm
        ( Normal, Normal, (T)-1, L21_MC_Star, X1_Star_MR, (T)1, X2 );
        //--------------------------------------------------------------------//
        L21_MC_Star.FreeAlignments();
        X1_Star_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::TrsmLLN
( Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrsmLLN
( Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrsmLLN
( Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrsmLLN
( Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

