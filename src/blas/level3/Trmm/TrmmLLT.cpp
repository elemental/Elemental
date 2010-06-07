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

// Left Lower (Conjugate)Transpose (Non)Unit Trmm
//   X := tril(L)^T,
//   X := tril(L)^H,
//   X := trilu(L)^T, or
//   X := trilu(L)^H
template<typename T>
void
elemental::blas::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  T alpha, 
  const DistMatrix<T,MC,MR>& L,
        DistMatrix<T,MC,MR>& X   )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmLLT");
    if( L.GetGrid() != X.GetGrid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( orientation == Normal )
        throw logic_error( "TrmmLLT expects a (Conjugate)Transpose option." );
    if( L.Height() != L.Width() || L.Height() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmLLT: " << endl
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

    DistMatrix<T,MC,MR> XT(g),  X0(g),
                        XB(g),  X1(g),
                                X2(g);

    // Temporary distributions
    DistMatrix<T,Star,Star> L11_Star_Star(g);
    DistMatrix<T,MC,  Star> L21_MC_Star(g);
    DistMatrix<T,Star,VR  > X1_Star_VR(g);
    DistMatrix<T,Star,MR  > D1_Star_MR(g);

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
        D1_Star_MR.AlignWith( X1 );
        D1_Star_MR.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_Star_VR = X1;
        L11_Star_Star = L11;
        blas::internal::LocalTrmm
        ( Left, Lower, orientation, diagonal, (T)1, L11_Star_Star, X1_Star_VR );
        X1 = X1_Star_VR;
 
        L21_MC_Star = L21;
        blas::internal::LocalGemm
        ( orientation, Normal, (T)1, L21_MC_Star, X2, (T)0, D1_Star_MR );
        X1.SumScatterUpdate( (T)1, D1_Star_MR );
       //--------------------------------------------------------------------//
        L21_MC_Star.FreeAlignments();
        D1_Star_MR.FreeAlignments();

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

template void elemental::blas::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmLLT
( Orientation orientation, 
  Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X );
#endif

