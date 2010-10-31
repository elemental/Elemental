/*
   Copyright (c) 2009-2010, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
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
            alpha = alpha11.GetLocalEntry(0,0);
            alpha11.SetLocalEntry(0,0,1);
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
            alpha11.SetLocalEntry(0,0,alpha);
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
            alpha = alpha11.GetLocalEntry(0,0);
            alpha11.SetLocalEntry(0,0,1);
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
            alpha11.SetLocalEntry(0,0,alpha);
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

