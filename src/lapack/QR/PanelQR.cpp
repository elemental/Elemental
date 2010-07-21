/*
   Copyright (C) 2009-2010 Jack Poulson <jack.poulson@gmail.com>

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with Elemental.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "elemental/blas_internal.hpp"
#include "elemental/lapack_internal.hpp"
using namespace std;
using namespace elemental;

template<typename R>
void
elemental::lapack::internal::PanelQR
( DistMatrix<R,MC,MR>& A )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelQR");
#endif
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<R,MC,MR>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ALeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);

    // Temporary distributions
    DistMatrix<R,MC,Star> ALeftCol_MC_Star(g);
    DistMatrix<R,MR,Star> Z_MR_Star(g);

    PushBlocksizeStack( 1 );
    PartitionDownLeftDiagonal
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
        R tau = lapack::Reflector( alpha11, a21 );

        bool myDiagonalEntry = ( g.MCRank() == alpha11.ColAlignment() && 
                                 g.MRRank() == alpha11.RowAlignment() );
        R alpha = (R)0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.LocalEntry(0,0);
            alpha11.LocalEntry(0,0) = (R)1;
        }

        ALeftCol_MC_Star = ALeftCol;

        blas::Gemv
        ( Transpose, 
          (R)1, ARightPan.LockedLocalMatrix(), 
                ALeftCol_MC_Star.LockedLocalMatrix(),
          (R)0, Z_MR_Star.LocalMatrix() );
        Z_MR_Star.SumOverCol(); 

        blas::Ger
        ( -tau, 
          ALeftCol_MC_Star.LockedLocalMatrix(), 
          Z_MR_Star.LockedLocalMatrix(),
          ARightPan.LocalMatrix() );

        if( myDiagonalEntry )
            alpha11.LocalEntry(0,0) = alpha;
        //--------------------------------------------------------------------//
        ALeftCol_MC_Star.FreeAlignments();
        Z_MR_Star.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::PanelQR
( DistMatrix<complex<R>,MC,MR  >& A,
  DistMatrix<complex<R>,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelQR");
    if( A.GetGrid() != t.GetGrid() )
        throw logic_error( "A and t must be distributed over the same grid." );
    if( t.Height() != min(A.Height(),A.Width()) || t.Width() != 1 )
        throw logic_error
              ( "t must be a column vector of height equal to the minimum "
                "dimension of A." );
    if( !t.AlignedWithDiag( A, 0 ) )
        throw logic_error( "t must be aligned with A's main diagonal." );
#endif
    typedef complex<R> C;
    const Grid& g = A.GetGrid();

    // Matrix views
    DistMatrix<C,MC,MR>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ALeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);
    DistMatrix<C,MD,Star>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    DistMatrix<C,MC,Star> ALeftCol_MC_Star(g);
    DistMatrix<C,MR,Star> Z_MR_Star(g);

    PushBlocksizeStack( 1 );
    PartitionDownLeftDiagonal
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
        C tau = lapack::Reflector( alpha11, a21 );
        tau1.Set( 0, 0, tau );

        bool myDiagonalEntry = ( g.MCRank() == alpha11.ColAlignment() && 
                                 g.MRRank() == alpha11.RowAlignment() );
        C alpha = (C)0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.LocalEntry(0,0);
            alpha11.LocalEntry(0,0) = (C)1;
        }

        ALeftCol_MC_Star = ALeftCol;

        blas::Gemv
        ( ConjugateTranspose, 
          (C)1, ARightPan.LockedLocalMatrix(), 
                ALeftCol_MC_Star.LockedLocalMatrix(),
          (C)0, Z_MR_Star.LocalMatrix() );
        Z_MR_Star.SumOverCol(); 

        blas::Ger
        ( -conj(tau), 
          ALeftCol_MC_Star.LockedLocalMatrix(), 
          Z_MR_Star.LockedLocalMatrix(),
          ARightPan.LocalMatrix() );

        if( myDiagonalEntry )
            alpha11.LocalEntry(0,0) = alpha;
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
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif

template void
elemental::lapack::internal::PanelQR
( DistMatrix<float,MC,MR>& A );

template void
elemental::lapack::internal::PanelQR
( DistMatrix<double,MC,MR>& A );

#ifndef WITHOUT_COMPLEX
template void
elemental::lapack::internal::PanelQR
( DistMatrix<scomplex,MC,MR  >& A,
  DistMatrix<scomplex,MD,Star>& t );

template void
elemental::lapack::internal::PanelQR
( DistMatrix<dcomplex,MC,MR  >& A,
  DistMatrix<dcomplex,MD,Star>& t );
#endif

