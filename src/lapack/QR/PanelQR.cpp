/*
   This file is part of elemental, a library for distributed-memory dense 
   linear algebra.

   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This program is released under the terms of the license contained in the 
   file LICENSE.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

template<typename T>
void
elemental::lapack::internal::PanelQR
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelQR");
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ALeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);

    // Temporary distributions
    DistMatrix<T,MC,Star> ALeftCol_MC_Star(g);
    DistMatrix<T,MR,Star> Z_MR_Star(g);

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        ALeftCol.View2x1
        ( alpha11,
          a21 );

        ARightPan.View2x1
        ( a12,
          A22 );

        ALeftCol_MC_Star.AlignWith( ARightPan );
        Z_MR_Star.AlignWith( ARightPan );
        Z_MR_Star.ResizeTo( ARightPan.Width(), 1 );
        //--------------------------------------------------------------------//
        T tau = lapack::internal::Reflector( alpha11, a21 );
        ALeftCol_MC_Star = ALeftCol;
        blas::Gemv
        ( ConjugateTranspose, 
          (T)1, ARightPan.LockedLocalMatrix(), ALeftCol.LockedLocalMatrix(),
          (T)0, Z_MR_Star.LocalMatrix() );
        Z_MR_Star.SumOverCol(); 

        blas::Ger
        ( -tau, 
          ALeftCol_MC_Star.LockedLocalMatrix(), 
          Z_MR_Star.LockedLocalMatrix(),
          ARightPan.LocalMatrix() );
        //--------------------------------------------------------------------//
        ALeftCol_MC_Star.FreeAlignments();
        Z_MR_Star.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::internal::PanelQR
( DistMatrix<float,MC,MR>& A );

template void
elemental::lapack::internal::PanelQR
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::internal::PanelQR
( DistMatrix<scomplex,MC,MR>& A );

template void
elemental::lapack::internal::PanelQR
( DistMatrix<dcomplex,MC,MR>& A );
#endif

