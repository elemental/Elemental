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
elemental::lapack::internal::TridiagUSquare
( DistMatrix<R,MC,MR>& paddedA, int padding )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TridiagUSquare");
    if( paddedA.Height() != paddedA.Width() )
        throw logic_error("A must be square.");
#endif
    const Grid& g = paddedA.Grid();

    // Separate off the padded parts of A and ensure that they are zero
    DistMatrix<R,MC,MR>
        paddedATL(g), paddedATR(g),
        paddedABL(g), A(g);
    PartitionDownLeftDiagonal
    ( paddedA, paddedATL, paddedATR,
               paddedABL, A,         padding );
    paddedATL.SetToZero();
    paddedATR.SetToZero();
    paddedABL.SetToZero();

#ifndef RELEASE
    if( g.Height() != g.Width() ||
        paddedA.ColAlignment() != paddedA.RowAlignment() )
        throw logic_error("Square Tridiag is for square, diag aligned grids.");
    if( paddedA.Height() % g.Height() != 0 )
        throw logic_error("Square Tridiag requires a padded matrix.");
#endif

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ABR(g), paddedA00(g), paddedA01(g), paddedA02(g),  
                paddedA10(g), A11(g),       A12(g),
                paddedA20(g), A21(g),       A22(g),
        ATL(g), A11Expanded(g), A01(g), A00(g);

    // Temporary distributions
    DistMatrix<R,MC,  Star> A01_MC_Star(g);
    DistMatrix<R,MR,  Star> A01_MR_Star(g);
    DistMatrix<R,MC,  Star> W01_MC_Star(g);
    DistMatrix<R,MR,  Star> W01_MR_Star(g);
    DistMatrix<R,Star,Star> A11_Star_Star(g);
    DistMatrix<R,MD,  Star> e1(g);
    DistMatrix<R,MC,  MR  > W01(g), W11(g),  WPan(g);

    PartitionUpDiagonal
    ( paddedA, paddedATL, paddedATR,
               paddedABL, ABR,       0 );
    while( ABR.Height()+padding < paddedA.Height() )
    {
        int bsize = min( Blocksize(), paddedA.Height()-(ABR.Height()+padding) );
        RepartitionUpDiagonal
        ( paddedATL, /**/ paddedATR,  paddedA00, paddedA01, /**/ paddedA02,
                     /**/             paddedA10, A11,       /**/ A12,
         /*************************/ /************************************/
          paddedABL, /**/ ABR,        paddedA20, A21,       /**/ A22,
          bsize );

        ATL.View
        ( paddedATL, padding, padding,
          paddedATL.Height()-padding, paddedATL.Width()-padding );
        A01.View
        ( paddedA01, padding, 0,
          paddedA01.Height()-padding, paddedA01.Width() );
        A00.View
        ( paddedA00, padding, padding,
          paddedA00.Height()-padding, paddedA00.Width()-padding );

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
            e1.AlignWithDiag( A11Expanded, 1 );
            e1.ResizeTo( WPan.Width(), 1 );
            A01_MC_Star.AlignWith( A00 );
            A01_MR_Star.AlignWith( A00 );
            W01_MC_Star.AlignWith( A00 );
            W01_MR_Star.AlignWith( A00 );
            //----------------------------------------------------------------//
            lapack::internal::PanelTridiagUSquare
            ( paddedATL, WPan, e1, padding );

            // Perform the Syr2k for square process grids on A00
            A01_MC_Star = A01; // Allgather row 
            A01_MR_Star = A01; // Pairwise exchange then Allgather column
            W01_MC_Star = W01; // Allgather row
            W01_MR_Star = W01; // Pairwise exchange then Allgather column
            blas::internal::LocalTriangularRank2K
            ( Upper, Transpose, Transpose, (R)-1,
              A01_MC_Star, W01_MC_Star, A01_MR_Star, W01_MR_Star, (R)1, A00 );

            A11Expanded.SetDiagonal( e1, 1 );
            //----------------------------------------------------------------//
            WPan.FreeAlignments();
            e1.FreeAlignments();
            A01_MC_Star.FreeAlignments();
            A01_MR_Star.FreeAlignments();
            W01_MC_Star.FreeAlignments();
            W01_MR_Star.FreeAlignments();
        }
        else
        {
            A11_Star_Star = A11;
            lapack::Tridiag( Upper, A11_Star_Star.LocalMatrix() );
            A11 = A11_Star_Star;
        }

        SlidePartitionUpDiagonal
        ( paddedATL, /**/ paddedATR,  paddedA00, /**/ paddedA01, paddedA02,
         /*************************/ /************************************/
                     /**/             paddedA10, /**/ A11,       A12,
          paddedABL, /**/ ABR,        paddedA20, /**/ A21,       A22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::TridiagUSquare
( DistMatrix<complex<R>,MC,  MR  >& paddedA,
  DistMatrix<complex<R>,Star,Star>& t,
  int padding )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TridiagUSquare");
    if( A.Height() != A.Width() )
        throw logic_error("A must be square.");
#endif
    const Grid& g = paddedA.Grid();
    typedef complex<R> C;

    // Separate off the padded parts of A and ensure that they are zero
    DistMatrix<C,MC,MR>
        paddedATL(g), paddedATR(g),
        paddedABL(g), A(g);
    PartitionDownLeftDiagonal
    ( paddedA, paddedATL, paddedATR,
               paddedABL, A,         padding );
    paddedATL.SetToZero();
    paddedATR.SetToZero();
    paddedABL.SetToZero();

#ifndef RELEASE
    if( paddedA.Grid() != t.Grid() )
        throw logic_error("A and t must be distributed over the same grid.");
    if( t.Viewing() )
        throw logic_error("t must not be a view.");
    if( g.Height() != g.Width() ||
        paddedA.ColAlignment() != paddedA.RowAlignment() )
        throw logic_error("Square Tridiag is for square, diag aligned grids.");
    if( paddedA.Height() % g.Height() != 0 )
        throw logic_error("Square Tridiag requires a padded matrix.");
#endif

    DistMatrix<C,MD,Star> tDiag(g);
    tDiag.AlignWithDiag( A, 1 );
    tDiag.ResizeTo( A.DiagonalLength( 1 ), 1 );

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ABR(g),  paddedA00(g), paddedA01(g), paddedA02(g),  
                 paddedA10(g), A11(g),       A12(g),
                 paddedA20(g), A21(g),       A22(g),
        ATL(g), A11Expanded(g), A01(g), A00(g);
    DistMatrix<C,MD,Star> tT(g),  t0(g),
                          tB(g),  t1(g), 
                                  t2(g);

    // Temporary distributions
    Matrix<C> A11_Herm;
    DistMatrix<C,MC,  Star> A01_MC_Star(g);
    DistMatrix<C,MR,  Star> A01_MR_Star(g);
    DistMatrix<C,MC,  Star> W01_MC_Star(g);
    DistMatrix<C,MR,  Star> W01_MR_Star(g);
    DistMatrix<C,Star,Star> A11_Star_Star(g);
    DistMatrix<R,MD,  Star> e1(g);
    DistMatrix<C,Star,Star> t1_Star_Star(g);
    DistMatrix<C,MC,  MR  > W01(g), W11(g),  WPan(g);

    PartitionUpDiagonal
    ( paddedA, paddedATL, paddedATR,
               paddedABL, ABR,       0 );
    PartitionUp
    ( tDiag, tT,
             tB, 0 );
    while( ABR.Height()+padding < paddedA.Height() )
    {
        int bsize = min( Blocksize(), paddedA.Height()-(ABR.Height()+padding) );
        RepartitionUpDiagonal
        ( paddedATL, /**/ paddedATR,  paddedA00, paddedA01, /**/ paddedA02,
                     /**/             paddedA10, A11,       /**/ A12,
         /*************************/ /************************************/
          paddedABL, /**/ ABR,        paddedA20, A21,       /**/ A22,
          bsize );

        RepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        ATL.View
        ( paddedATL, padding, padding,
          paddedATL.Height()-padding, paddedATL.Width()-padding );
        A01.View
        ( paddedA01, padding, 0,
          paddedA01.Height()-padding, paddedA01.Width() );
        A00.View
        ( paddedA00, padding, padding,
          paddedA00.Height()-padding, paddedA00.Width()-padding );

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
            e1.AlignWithDiag( A11Expanded, 1 );
            e1.ResizeTo( WPan.Width(), 1 );
            A01_MC_Star.AlignWith( A00 );
            A01_MR_Star.AlignWith( A00 );
            W01_MC_Star.AlignWith( A00 );
            W01_MR_Star.AlignWith( A00 );
            //----------------------------------------------------------------//
            lapack::internal::PanelTridiagUSquare
            ( paddedATL, WPan, e1, t1, padding );

            // Perform the Her2k for square process grids on A00
            A01_MC_Star = A01; // Allgather row
            A01_MR_Star = A01; // Pairwise exchange then Allgather column
            W01_MC_Star = W01; // Allgather row
            W01_MR_Star = W01; // Pairwise exchange then Allgather column
            blas::internal::LocalTriangularRank2K
            ( Upper, ConjugateTranspose, ConjugateTranspose, (C)-1,
              A01_MC_Star, W01_MC_Star, A01_MR_Star, W01_MR_Star, (C)1, A00 );

            A11Expanded.SetDiagonal( e1, 1 );
            //----------------------------------------------------------------//
            WPan.FreeAlignments();
            e1.FreeAlignments();
            A01_MC_Star.FreeAlignments();
            A01_MR_Star.FreeAlignments();
            W01_MC_Star.FreeAlignments();
            W01_MR_Star.FreeAlignments();
        }
        else
        {
            A11_Star_Star = A11;
            t1_Star_Star.ResizeTo( t1.Height(), 1 );

            lapack::Tridiag
            ( Upper, 
              A11_Star_Star.LocalMatrix(),        
              t1_Star_Star.LocalMatrix() );
            
            A11 = A11_Star_Star;
            t1 = t1_Star_Star;
        }

        SlidePartitionUpDiagonal
        ( paddedATL, /**/ paddedATR,  paddedA00, /**/ paddedA01, paddedA02,
         /*************************/ /************************************/
                     /**/             paddedA10, /**/ A11,       A12,
          paddedABL, /**/ ABR,        paddedA20, /**/ A21,       A22 );

        SlidePartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );
    }
    // Redistribute from matrix-diag form to fully replicated
    t = tDiag;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::internal::TridiagUSquare
( DistMatrix<float,MC,MR>& paddedA, int padding );

template void elemental::lapack::internal::TridiagUSquare
( DistMatrix<double,MC,MR>& paddedA, int padding );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::TridiagUSquare
( DistMatrix<scomplex,MC,  MR  >& paddedA,
  DistMatrix<scomplex,Star,Star>& t,
  int padding );

template void elemental::lapack::internal::TridiagUSquare
( DistMatrix<dcomplex,MC,  MR  >& paddedA,
  DistMatrix<dcomplex,Star,Star>& t,
  int padding );
#endif

