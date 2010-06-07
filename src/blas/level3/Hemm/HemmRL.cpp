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

template<typename T>
void
elemental::blas::internal::HemmRL
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmRL");
#endif
    blas::internal::HemmRLC( alpha, A, B, beta, C );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
void
elemental::blas::internal::HemmRLC
( T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("blas::internal::HemmRLC");
    if( A.GetGrid() != B.GetGrid() || B.GetGrid() != C.GetGrid() )
        throw logic_error( "{A,B,C} must be distributed over the same grid." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  AColPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),  ARowPan(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    DistMatrix<T,MC,MR> CL(g), CR(g),
                        C0(g), C1(g), C2(g),
                        CLeft(g), CRight(g);

    // Temporary distributions
    DistMatrix<T,MC,Star> B1_MC_Star(g);
    DistMatrix<T,MR,Star> AColPan_MR_Star(g);
    DistMatrix<T,Star,MR> ARowPan_Star_MR(g);

    // Start the algorithm
    blas::Scal( beta, C );
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    PartitionRight( C, CL, CR, 0 );
    while( CR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        RepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );

        ARowPan.LockedView1x2( A10, A11 );

        AColPan.LockedView2x1
        ( A11,
          A21 );

        CLeft.View1x2( C0, C1 );

        CRight.View1x2( C1, C2 );

        B1_MC_Star.AlignWith( C );
        AColPan_MR_Star.AlignWith( CRight );
        ARowPan_Star_MR.AlignWith( CLeft );
        //--------------------------------------------------------------------//
        B1_MC_Star = B1;

        ARowPan_Star_MR = ARowPan;
        AColPan_MR_Star = AColPan;
        ARowPan_Star_MR.MakeTrapezoidal( Right, Lower );
        AColPan_MR_Star.MakeTrapezoidal( Left, Lower, -1 );

        blas::internal::LocalGemm
        ( Normal, Normal, alpha, B1_MC_Star, ARowPan_Star_MR, (T)1, CLeft );

        blas::internal::LocalGemm
        ( Normal, ConjugateTranspose, 
          alpha, B1_MC_Star, AColPan_MR_Star, (T)1, CRight );
        //--------------------------------------------------------------------//
        B1_MC_Star.FreeAlignments();
        AColPan_MR_Star.FreeAlignments();
        ARowPan_Star_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );

        SlidePartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::blas::internal::HemmRL
( float alpha, const DistMatrix<float,MC,MR>& A,
               const DistMatrix<float,MC,MR>& B,
  float beta,        DistMatrix<float,MC,MR>& C );

template void elemental::blas::internal::HemmRL
( double alpha, const DistMatrix<double,MC,MR>& A,
                const DistMatrix<double,MC,MR>& B,
  double beta,        DistMatrix<double,MC,MR>& C );

#ifndef WITHOUT_COMPLEX
template void elemental::blas::internal::HemmRL
( scomplex alpha, const DistMatrix<scomplex,MC,MR>& A,
                  const DistMatrix<scomplex,MC,MR>& B,
  scomplex beta,        DistMatrix<scomplex,MC,MR>& C );

template void elemental::blas::internal::HemmRL
( dcomplex alpha, const DistMatrix<dcomplex,MC,MR>& A,
                  const DistMatrix<dcomplex,MC,MR>& B,
  dcomplex beta,        DistMatrix<dcomplex,MC,MR>& C );
#endif

