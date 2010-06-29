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
elemental::lapack::internal::TridiagL
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MD,Star>& d,
  DistMatrix<R,MD,Star>& e )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TridiagL");
#endif
    const Grid& g = A.GetGrid();
#ifndef RELEASE
    if( A.GetGrid() != d.GetGrid() || d.GetGrid() != e.GetGrid() )
        throw logic_error
        ( "A, d, and e must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( d.Viewing() && ( d.Height() != A.Height() || d.Width() != 1 ) )
        throw logic_error
        ( "d must be a column vector of the same length as A's width." );
    if( e.Viewing() && ( e.Height() != A.Height()-1 || e.Width() != 1 ) )
        throw logic_error
        ( "e must be a column vector of length one less than the width of A." );
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) && 
        ( d.ColAlignment() != A.ColAlignment() + A.RowAlignment()*g.Height() ) )
        throw logic_error( "d is not aligned with A." );
    if( ( e.Viewing() || e.ConstrainedColAlignment() ) && 
        ( e.ColAlignment() != ((A.ColAlignment()+1) % g.Height())
                              + A.RowAlignment() * g.Height() ) )
        throw logic_error( "e is not aligned with A." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiag( A );
        d.ResizeTo( A.Height(), 1 );
    }
    if( !e.Viewing() )
    {
        if( !e.ConstrainedColAlignment() )
            e.AlignWithDiag( A, -1 );
        e.ResizeTo( A.Height()-1, 1 );
    }

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g),
        A11Expanded(g);
    DistMatrix<R,MD,Star> dT(g),  d0(g), 
                          dB(g),  d1(g),
                                  d2(g);
    DistMatrix<R,MD,Star> eT(g),  e0(g), 
                          eB(g),  e1(g), 
                                  e2(g);

    // Temporary distributions
    DistMatrix<R,Star,Star> A11_Star_Star(g);
    DistMatrix<R,Star,Star> d1_Star_Star(g);
    DistMatrix<R,Star,Star> e1_Star_Star(g);
    DistMatrix<R,MC,  MR  > W11(g),  WPan(g),
                            W21(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( d, dT,
         dB, 0 );
    PartitionDown
    ( e, eT,
         eB, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( dT,  d0,
         /**/ /**/
               d1,
          dB,  d2 );

        RepartitionDown
        ( eT,  e0,
         /**/ /**/
               e1,
          eB,  e2 );

        if( A22.Height() > 0 )
        {
            A11Expanded.View( ABR, 0, 0, A11.Height()+1, A11.Width()+1 );
            WPan.AlignWith( A11 );
            WPan.ResizeTo( ABR.Height(), A11.Width() );
            PartitionDown
            ( WPan, W11,
                    W21, A11.Height() );
            //----------------------------------------------------------------//
            lapack::internal::PanelTridiagL( ABR, WPan, e1 );
            blas::Syr2k( Lower, Normal, (R)-1, A21, W21, (R)1, A22 );
            A11Expanded.SetDiagonal( e1, -1 );
            A11.GetDiagonal( d1 );
            //----------------------------------------------------------------//
            WPan.FreeAlignments();
        }
        else
        {
            A11_Star_Star = A11;
            d1_Star_Star.ResizeTo( d1.Height(), 1 );
            e1_Star_Star.ResizeTo( e1.Height(), 1 );

            lapack::Tridiag
            ( Lower, 
              A11_Star_Star.LocalMatrix(),
              d1_Star_Star.LocalMatrix(), 
              e1_Star_Star.LocalMatrix() );

            A11 = A11_Star_Star;
            d1 = d1_Star_Star;
            e1 = e1_Star_Star;
        }

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDown
        ( dT,  d0,
               d1,
         /**/ /**/
          dB,  d2 );

        SlidePartitionDown
        ( eT,  e0,
               e1,
         /**/ /**/
          eB,  e2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::TridiagL
( DistMatrix<complex<R>,MC,MR  >& A,
  DistMatrix<R,MD,Star>& d,
  DistMatrix<R,MD,Star>& e )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TridiagL");
#endif
    const Grid& g = A.GetGrid();
#ifndef RELEASE
    if( A.GetGrid() != d.GetGrid() || d.GetGrid() != e.GetGrid() )
        throw logic_error
        ( "A, d, and e must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( d.Viewing() && ( d.Height() != A.Height() || d.Width() != 1 ) )
        throw logic_error
        ( "d must be a column vector of the same length as A's width." );
    if( e.Viewing() && ( e.Height() != A.Height()-1 || e.Width() != 1 ) )
        throw logic_error
        ( "e must be a column vector of length one less than the width of A." );
    if( ( d.Viewing() || d.ConstrainedColAlignment() ) && 
        ( d.ColAlignment() != A.ColAlignment() + A.RowAlignment()*g.Height() ) )
        throw logic_error( "d is not aligned with A." );
    if( ( e.Viewing() || e.ConstrainedColAlignment() ) && 
        ( e.ColAlignment() != ((A.ColAlignment()+1) % g.Height())
                              + A.RowAlignment() * g.Height() ) )
        throw logic_error( "e is not aligned with A." );
#endif
    if( !d.Viewing() )
    {
        if( !d.ConstrainedColAlignment() )
            d.AlignWithDiag( A );
        d.ResizeTo( A.Height(), 1 );
    }
    if( !e.Viewing() )
    {
        if( !e.ConstrainedColAlignment() )
            e.AlignWithDiag( A, -1 );
        e.ResizeTo( A.Height()-1, 1 );
    }
    typedef complex<R> C;

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g), 
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g),
        A11Expanded(g);
    DistMatrix<R,MD,Star> dT(g),  d0(g), 
                          dB(g),  d1(g),
                                  d2(g);
    DistMatrix<R,MD,Star> eT(g),  e0(g), 
                          eB(g),  e1(g), 
                                  e2(g);

    // Temporary distributions
    DistMatrix<C,Star,Star> A11_Star_Star(g);
    DistMatrix<R,Star,Star> d1_Star_Star(g);
    DistMatrix<R,Star,Star> e1_Star_Star(g);
    DistMatrix<C,MC,  MR  > W11(g),  WPan(g),
                            W21(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDown
    ( d, dT,
         dB, 0 );
    PartitionDown
    ( e, eT,
         eB, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        RepartitionDown
        ( dT,  d0,
         /**/ /**/
               d1,
          dB,  d2 );

        RepartitionDown
        ( eT,  e0,
         /**/ /**/
               e1,
          eB,  e2 );

        if( A22.Height() > 0 )
        {
            A11Expanded.View( ABR, 0, 0, A11.Height()+1, A11.Width()+1 );
            WPan.AlignWith( A11 );
            WPan.ResizeTo( ABR.Height(), A11.Width() );
            PartitionDown
            ( WPan, W11,
                    W21, A11.Height() );
            //----------------------------------------------------------------//
            lapack::internal::PanelTridiagL( ABR, WPan, e1 );
            blas::Her2k( Lower, Normal, (C)-1, A21, W21, (C)1, A22 );
            A11Expanded.SetDiagonal( e1, -1 );
            A11.GetRealDiagonal( d1 );
            //----------------------------------------------------------------//
            WPan.FreeAlignments();
        }
        else
        {
            A11_Star_Star = A11;
            d1_Star_Star.ResizeTo( d1.Height(), 1 );
            e1_Star_Star.ResizeTo( e1.Height(), 1 );

            lapack::Tridiag
            ( Lower, 
              A11_Star_Star.LocalMatrix(),
              d1_Star_Star.LocalMatrix(), 
              e1_Star_Star.LocalMatrix() );

            A11 = A11_Star_Star;
            d1 = d1_Star_Star;
            e1 = e1_Star_Star;
        }

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlidePartitionDown
        ( dT,  d0,
               d1,
         /**/ /**/
          dB,  d2 );

        SlidePartitionDown
        ( eT,  e0,
               e1,
         /**/ /**/
          eB,  e2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::internal::TridiagL
( DistMatrix<float,MC,MR  >& A,
  DistMatrix<float,MD,Star>& d,
  DistMatrix<float,MD,Star>& e );

template void elemental::lapack::internal::TridiagL
( DistMatrix<double,MC,MR  >& A,
  DistMatrix<double,MD,Star>& d,
  DistMatrix<double,MD,Star>& e );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::TridiagL
( DistMatrix<scomplex,MC,MR  >& A,
  DistMatrix<float,   MD,Star>& d,
  DistMatrix<float,   MD,Star>& e );

template void elemental::lapack::internal::TridiagL
( DistMatrix<dcomplex,MC,MR  >& A,
  DistMatrix<double,  MD,Star>& d,
  DistMatrix<double,  MD,Star>& e );
#endif

