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
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::QR");
    if( A.GetGrid() != t.GetGrid() )
        throw logic_error( "A and t must be distributed over the same grid." );
    if( t.Viewing() && 
        ( t.Height() != min(A.Height(),A.Width()) || t.Width() != 1 ) )
        throw logic_error
              ( "t must be a column vector of the same length as the minimum "
                "dimension of A." );
    if( ( t.Viewing() || t.ConstrainedColAlignment() ) && 
        ( t.ColAlignment() != A.ColAlignment() + 
                              A.RowAlignment()*A.GetGrid().Height() ) )
        throw logic_error( "t is not aligned with A." );
#endif
    if( !t.Viewing() ) 
    {
        if( !t.ConstrainedColAlignment() )
            t.AlignWithDiag( A );
        t.ResizeTo( min(A.Height(),A.Width()), 1 );
    }
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),  ALeftPan(g), ARightPan(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);
    DistMatrix<T,MD,Star>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ( ATL.Height() < A.Height() && ATL.Width() < A.Width() ) )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        ALeftPan.View2x1
        ( A11,
          A21 );

        ARightPan.View2x1
        ( A12,
          A22 );

        //--------------------------------------------------------------------//
        lapack::internal::PanelQR( ALeftPan, t1 );
        lapack::UT( Lower, Normal, 0, ALeftPan, ARightPan );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::QR
( DistMatrix<float,MC,MR>& A, DistMatrix<float,MD,Star>& t );

template void
elemental::lapack::QR
( DistMatrix<double,MC,MR>& A, DistMatrix<double,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::QR
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<scomplex,MD,Star>& t );

template void
elemental::lapack::QR
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<dcomplex,MD,Star>& t );
#endif

