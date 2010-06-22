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
( DistMatrix<T,MC,MR>& A, DistMatrix<T,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelQR");
    if( A.GetGrid() != t.GetGrid() )
        throw logic_error( "A and t must be distributed over the same grid." );
    if( t.Height() != min(A.Height(),A.Width()) || t.Width() != 1 )
        throw logic_error
              ( "t must be a column vector of the same length as the minimum "
                "dimension of A." );
    if( t.ColAlignment() != A.ColAlignment() + 
                            A.RowAlignment()*A.GetGrid().Height() )
        throw logic_error( "t is not aligned with A." );
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<T,MC,MR>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ALeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);

    DistMatrix<T,MD,Star>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    DistMatrix<T,MC,Star> ALeftCol_MC_Star(g);
    DistMatrix<T,MR,Star> Z_MR_Star(g);

    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionDown
        ( tT,  t0,
         /**/ /****/
               tau1,
          tB,  t2 );

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

        tau1.Set( 0, 0, tau );
        //--------------------------------------------------------------------//
        ALeftCol_MC_Star.FreeAlignments();
        Z_MR_Star.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void
elemental::lapack::internal::PanelQR
( DistMatrix<float,MC,MR>& A, DistMatrix<float,MD,Star>& t );

template void
elemental::lapack::internal::PanelQR
( DistMatrix<double,MC,MR>& A, DistMatrix<double,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::internal::PanelQR
( DistMatrix<scomplex,MC,MR>& A, DistMatrix<scomplex,MD,Star>& t );

template void
elemental::lapack::internal::PanelQR
( DistMatrix<dcomplex,MC,MR>& A, DistMatrix<dcomplex,MD,Star>& t );
#endif

