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

// Left Lower Normal (Non)Unit Trmm 
//   X := tril(L)  X, or
//   X := trilu(L) X
template<typename T>
void
elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmLLN");
    if( L.GetGrid() != X.GetGrid() )
        throw "L and X must be distributed over the same grid.";
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLLN: " << endl
            << "  L ~ " << L.Height() << " x " << L.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        const string& s = msg.str();
        throw s.c_str();
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
    DistMatrix<T,Star,MC  > L10_Star_MC(grid);
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,Star,VR  > X1_Star_VR(grid);
    DistMatrix<T,Star,MR  > D1_Star_MR(grid);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionUpDiagonal
    ( L, LTL, LTR,
         LBL, LBR );
    PartitionUp
    ( X, XT,
         XB );
    while( XT.Height() > 0 )
    {
        LockedRepartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        RepartitionUp
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );

        L10_Star_MC.AlignWith( X0 );
        D1_Star_MR.AlignWith( X1 );
        D1_Star_MR.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        L11_Star_Star = L11;
        X1_Star_VR = X1;
        blas::internal::LocalTrmm
        ( Left, Lower, Normal, diagonal, (T)1, L11_Star_Star, X1_Star_VR );
        X1 = X1_Star_VR;

        L10_Star_MC = L10;
        blas::internal::LocalGemm
        ( Normal, Normal, (T)1, L10_Star_MC, X0, (T)0, D1_Star_MR );
        X1.SumScatterUpdate( (T)1, D1_Star_MR );
        //--------------------------------------------------------------------//
        L10_Star_MC.FreeAlignments();
        D1_Star_MR.FreeAlignments();

        SlideLockedPartitionUpDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12, 
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

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

template void elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  float alpha, const DistMatrix<float,MC,MR>& L,
                     DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  double alpha, const DistMatrix<double,MC,MR>& L,
                      DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& L,
                        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmLLN
( Diagonal diagonal,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& L,
                        DistMatrix<dcomplex,MC,MR>& X );
#endif

