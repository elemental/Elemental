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

// On exit, the upper triangle of A is overwritten by R, and the Householder
// transforms that determine Q are stored below the diagonal of A with an 
// implicit one on the diagonal. On exit, the column-vector t stores the 
// coefficients of the Householder transforms.
template<typename T>
void
elemental::lapack::QR
( DistMatrix<T,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::QR");
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ALeftPan(g), ARightPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ( ATL.Height() < A.Height() && ATL.Width() < A.Width() ) )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        ALeftPan.View2x1
        ( A11,
          A21 );

        ARightPan.View2x1
        ( A12,
          A22 );

        //--------------------------------------------------------------------//
        lapack::internal::PanelQR( ALeftPan );
        lapack::UT( Left, Lower, Normal, 0, ALeftPan, ARightPan );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::QR
( DistMatrix<float,MC,MR>& A );

template void
elemental::lapack::QR
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::QR
( DistMatrix<scomplex,MC,MR>& A );

template void
elemental::lapack::QR
( DistMatrix<dcomplex,MC,MR>& A );
#endif

