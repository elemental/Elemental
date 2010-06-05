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

// Right Lower Normal (Non)Unit Trmm
//   X := X tril(L), and
//   X := X trilu(L)
template<typename T>
void
elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  T alpha, const DistMatrix<T,MC,MR>& L,
                 DistMatrix<T,MC,MR>& X )
{
#ifndef RELEASE
    PushCallStack("blas::internal::TrmmRLN");
    if( L.GetGrid() != X.GetGrid() )
        throw "L and X must be distributed over the same grid.";
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrmmRLN: " << endl
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

    DistMatrix<T,MC,MR> XL(grid), XR(grid),
                        X0(grid), X1(grid), X2(grid);

    // Temporary distributions
    DistMatrix<T,Star,Star> L11_Star_Star(grid);
    DistMatrix<T,MR,  Star> L21_MR_Star(grid);
    DistMatrix<T,VC,  Star> X1_VC_Star(grid);
    DistMatrix<T,MC,  Star> D1_MC_Star(grid);

    // Start the algorithm
    blas::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR );
    PartitionRight( X, XL, XR );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );
 
        RepartitionRight
        ( XL, /**/ XR,
          X0, /**/ X1, X2 );

        L21_MR_Star.AlignWith( X2 );
        D1_MC_Star.AlignWith( X1 );
        D1_MC_Star.ResizeTo( X1.Height(), X1.Width() );
        //--------------------------------------------------------------------//
        X1_VC_Star = X1;
        L11_Star_Star = L11;
        blas::Trmm
        ( Right, Lower, Normal, diagonal,
          (T)1, L11_Star_Star.LockedLocalMatrix(),
                X1_VC_Star.LocalMatrix() );
        X1 = X1_VC_Star;
 
        L21_MR_Star = L21;
        blas::Gemm
        ( Normal, Normal, 
          (T)1, X2.LockedLocalMatrix(),
                L21_MR_Star.LockedLocalMatrix(),
          (T)0, D1_MC_Star.LocalMatrix() );
        X1.SumScatterUpdate( (T)1, D1_MC_Star );
       //--------------------------------------------------------------------//
        L21_MR_Star.FreeAlignments();
        D1_MC_Star.FreeAlignments();

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

template void elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  float alpha, const DistMatrix<float,MC,MR>& L,
                     DistMatrix<float,MC,MR>& X );

template void elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  double alpha, const DistMatrix<double,MC,MR>& L,
                      DistMatrix<double,MC,MR>& X );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  scomplex alpha, const DistMatrix<scomplex,MC,MR>& L,
                        DistMatrix<scomplex,MC,MR>& X );

template void elemental::blas::internal::TrmmRLN
( Diagonal diagonal,
  dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& L,
                        DistMatrix<dcomplex,MC,MR>& X );
#endif

