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

template<typename R>
void
elemental::lapack::internal::TridiagU
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MD,Star>& d,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<R,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TridiagU");
#endif
    const Grid& grid = A.GetGrid();
#ifndef RELEASE
    if( A.GetGrid() != d.GetGrid() ||
        d.GetGrid() != e.GetGrid() ||
        e.GetGrid() != t.GetGrid() )
        throw logic_error
        ( "A, d, e, and t must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( d.Height() != A.Height() || d.Width() != 1 )
        throw logic_error
        ( "d must be a column vector of the same length as A's width." );
    if( e.Height() != A.Height()-1 || e.Width() != 1 )
        throw logic_error
        ( "e must be a column vector of length one less than the width of A." );
    if( t.Height() != A.Height()-1 || t.Width() != 1 )
        throw logic_error
        ( "t must be a column vector of length one less than the width of A." );
    if( d.ColAlignment() != A.ColAlignment() + A.RowAlignment()*grid.Height() )
        throw logic_error( "d is not aligned with A." );
    if( e.ColAlignment() != A.ColAlignment() +
                            ((A.RowAlignment()+1)%grid.Width())*grid.Height() )
        throw logic_error( "e is not aligned with A." );
    if( t.ColAlignment() != ((A.ColAlignment()+1)%grid.Height()) + 
                            ((A.RowAlignment()+1)%grid.Width())*grid.Height() )
        throw logic_error( "t is not aligned with A." );
#endif

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),  
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),
                               A20(grid), A21(grid), A22(grid),
        A11Expanded(grid);
    DistMatrix<R,MD,Star> dT(grid),  d0(grid), 
                          dB(grid),  d1(grid),
                                     d2(grid);
    DistMatrix<R,MD,Star> eT(grid),  e0(grid),
                          eB(grid),  e1(grid), 
                                     e2(grid);
    DistMatrix<R,MD,Star> tT(grid),  t0(grid), 
                          tB(grid),  t1(grid),
                                     t2(grid);

    // Temporary distributions
    DistMatrix<R,Star,Star> A11_Star_Star(grid);
    DistMatrix<R,Star,Star> d1_Star_Star(grid);
    DistMatrix<R,Star,Star> e1_Star_Star(grid);
    DistMatrix<R,Star,Star> t1_Star_Star(grid);
    DistMatrix<R,MC,  MR  > W01(grid), W11(grid),  WPan(grid);

    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( d, dT,
         dB, 0 );
    PartitionUp
    ( e, eT,
         eB, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    while( ABR.Height() < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        RepartitionUp
        ( dT,  d0,
               d1,
         /**/ /**/
          dB,  d2 );

        RepartitionUp
        ( eT,  e0,
               e1,
         /**/ /**/
          eB,  e2 );

        RepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        if( A00.Height() > 0 )
        {
            A11Expanded.View
            ( ATL, 
              A00.Height()-1, A00.Width()-1, A11.Height()+1, A11.Width()+1 );
            WPan.AlignWith( A01 );
            WPan.ResizeTo( ATL.Height(), A11.Width() );
            PartitionUp
            ( WPan, W01,
                    W11, A11.Height() );
            //----------------------------------------------------------------//
            lapack::internal::PanelTridiagU( ATL, WPan, e1, t1 );
            blas::Syr2k( Upper, Normal, (R)-1, A01, W01, (R)1, A00 );
            A11Expanded.SetDiagonal( e1, 1 );
            A11.GetDiagonal( d1 );
            //----------------------------------------------------------------//
            WPan.FreeAlignments();
        }
        else
        {
            A11_Star_Star = A11;
            d1_Star_Star = d1;
            e1_Star_Star = e1;
            t1_Star_Star = t1;

            lapack::Tridiag
            ( Upper, 
              A11_Star_Star.LocalMatrix(),         
              d1_Star_Star.LocalMatrix(),
              e1_Star_Star.LocalMatrix(),
              t1_Star_Star.LocalMatrix() );
            
            A11 = A11_Star_Star;
            d1 = d1_Star_Star;
            e1 = e1_Star_Star;
            t1 = t1_Star_Star;
        }

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        SlidePartitionUp
        ( dT,  d0,
         /**/ /**/
               d1,
          dB,  d2 );

        SlidePartitionUp
        ( eT,  e0,
         /**/ /**/
               e1,
          eB,  e2 );

        SlidePartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );
    }
        
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::TridiagU
( DistMatrix<complex<R>,MC,MR  >& A,
  DistMatrix<R,         MD,Star>& d,
  DistMatrix<R,         MD,Star>& e,
  DistMatrix<complex<R>,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TridiagU");
#endif
    const Grid& grid = A.GetGrid();
#ifndef RELEASE
    if( A.GetGrid() != d.GetGrid() ||
        d.GetGrid() != e.GetGrid() ||
        e.GetGrid() != t.GetGrid() )
        throw logic_error
        ( "A, d, e, and t must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( d.Height() != A.Height() || d.Width() != 1 )
        throw logic_error
        ( "d must be a column vector of the same length as A's width." );
    if( e.Height() != A.Height()-1 || e.Width() != 1 )
        throw logic_error
        ( "e must be a column vector of length one less than the width of A." );
    if( t.Height() != A.Height()-1 || t.Width() != 1 )
        throw logic_error
        ( "t must be a column vector of length one less than the width of A." );
    if( d.ColAlignment() != A.ColAlignment() + A.RowAlignment()*grid.Height() )
        throw logic_error( "d is not aligned with A." );
    if( e.ColAlignment() != A.ColAlignment() +
                            ((A.RowAlignment()+1)%grid.Width())*grid.Height() )
        throw logic_error( "e is not aligned with A." );
    if( t.ColAlignment() != ((A.ColAlignment()+1)%grid.Height()) + 
                            ((A.RowAlignment()+1)%grid.Width())*grid.Height() )
        throw logic_error( "t is not aligned with A." );
#endif
    typedef complex<R> C;

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), A01(grid), A02(grid),  
        ABL(grid), ABR(grid),  A10(grid), A11(grid), A12(grid),
                               A20(grid), A21(grid), A22(grid),
        A11Expanded(grid);
    DistMatrix<R,MD,Star> dT(grid),  d0(grid), 
                          dB(grid),  d1(grid),
                                     d2(grid);
    DistMatrix<R,MD,Star> eT(grid),  e0(grid),
                          eB(grid),  e1(grid), 
                                     e2(grid);
    DistMatrix<C,MD,Star> tT(grid),  t0(grid), 
                          tB(grid),  t1(grid),
                                     t2(grid);

    // Temporary distributions
    Matrix<C> A11_Herm;
    DistMatrix<C,Star,Star> A11_Star_Star(grid);
    DistMatrix<R,Star,Star> d1_Star_Star(grid);
    DistMatrix<R,Star,Star> e1_Star_Star(grid);
    DistMatrix<C,Star,Star> t1_Star_Star(grid);
    DistMatrix<C,MC,  MR  > W01(grid), W11(grid),  WPan(grid);

    PartitionUpDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUp
    ( d, dT,
         dB, 0 );
    PartitionUp
    ( e, eT,
         eB, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    while( ABR.Height() < A.Height() )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        RepartitionUp
        ( dT,  d0,
               d1,
         /**/ /**/
          dB,  d2 );

        RepartitionUp
        ( eT,  e0,
               e1,
         /**/ /**/
          eB,  e2 );

        RepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        if( A00.Height() > 0 )
        {
            // HERE
            A11Expanded.View
            ( ATL, 
              A00.Height()-1, A00.Width()-1, A11.Height()+1, A11.Width()+1 );
            WPan.AlignWith( A01 );
            WPan.ResizeTo( ATL.Height(), A11.Width() );
            PartitionUp
            ( WPan, W01,
                    W11, A11.Height() );
            //----------------------------------------------------------------//
            lapack::internal::PanelTridiagU( ATL, WPan, e1, t1 );
            blas::Her2k( Upper, Normal, (C)-1, A01, W01, (C)1, A00 );
            A11Expanded.SetDiagonal( e1, 1 );
            A11.GetRealDiagonal( d1 );
            //----------------------------------------------------------------//
            WPan.FreeAlignments();
        }
        else
        {
            A11_Star_Star = A11;
            d1_Star_Star = d1;
            e1_Star_Star = e1;
            t1_Star_Star = t1;

            lapack::Tridiag
            ( Upper, 
              A11_Star_Star.LocalMatrix(),        
              d1_Star_Star.LocalMatrix(),
              e1_Star_Star.LocalMatrix(),
              t1_Star_Star.LocalMatrix() );
            
            A11 = A11_Star_Star;
            d1 = d1_Star_Star;
            e1 = e1_Star_Star;
            t1 = t1_Star_Star;
        }

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        SlidePartitionUp
        ( dT,  d0,
         /**/ /**/
               d1,
          dB,  d2 );

        SlidePartitionUp
        ( eT,  e0,
         /**/ /**/
               e1,
          eB,  e2 );

        SlidePartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );
    }
        
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::internal::TridiagU
( DistMatrix<float,MC,MR  >& A,
  DistMatrix<float,MD,Star>& d,
  DistMatrix<float,MD,Star>& e,
  DistMatrix<float,MD,Star>& t );

template void elemental::lapack::internal::TridiagU
( DistMatrix<double,MC,MR  >& A,
  DistMatrix<double,MD,Star>& d,
  DistMatrix<double,MD,Star>& e,
  DistMatrix<double,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::TridiagU
( DistMatrix<scomplex,MC,MR  >& A,
  DistMatrix<float,   MD,Star>& d,
  DistMatrix<float,   MD,Star>& e,
  DistMatrix<scomplex,MD,Star>& t );

template void elemental::lapack::internal::TridiagU
( DistMatrix<dcomplex,MC,MR  >& A,
  DistMatrix<double,  MD,Star>& d,
  DistMatrix<double,  MD,Star>& e,
  DistMatrix<dcomplex,MD,Star>& t );
#endif

