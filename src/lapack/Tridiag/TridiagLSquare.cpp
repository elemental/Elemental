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
elemental::lapack::internal::TridiagLSquare
( DistMatrix<R,MC,MR>& paddedA )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TridiagLSquare");
    if( paddedA.Height() != paddedA.Width() )
        throw logic_error( "paddedA must be square." );
#endif
    const Grid& g = paddedA.Grid();
    const int padding = g.Height();

    // Separate off the padded parts of A and ensure that they are zero
    DistMatrix<R,MC,MR>
        A(g),         paddedATR(g),
        paddedABL(g), paddedABR(g);
    PartitionUpLeftDiagonal
    ( paddedA, A,         paddedATR,
               paddedABL, paddedABR, padding );
    paddedATR.SetToZero();
    paddedABL.SetToZero();
    paddedABR.SetToZero();

#ifndef RELEASE
    if( g.Height() != g.Width() || 
        paddedA.ColAlignment() != paddedA.RowAlignment() )
        throw logic_error("Square Tridiag is for square, diag aligned grids.");
#endif

    if( g.InGrid() )
    {
        // Matrix views 
        DistMatrix<R,MC,MR> 
            ATL(g),  A00(g),       A01(g),       paddedA02(g), 
                     A10(g),       A11(g),       paddedA12(g),
                     paddedA20(g), paddedA21(g), paddedA22(g),
            ABR(g), A11Expanded(g), A21(g), A22(g);

        // Temporary distributions
        DistMatrix<R,MC,  Star> A21_MC_Star(g);
        DistMatrix<R,MR,  Star> A21_MR_Star(g);
        DistMatrix<R,MC,  Star> W21_MC_Star(g);
        DistMatrix<R,MR,  Star> W21_MR_Star(g);
        DistMatrix<R,Star,Star> A11_Star_Star(g);
        DistMatrix<R,MD,  Star> e1(g);
        DistMatrix<R,MC,  MR  > W11(g),  WPan(g),
                                W21(g);

        PartitionDownDiagonal
        ( paddedA, ATL,       paddedATR,
                   paddedABL, paddedABR, 0 );
        while( ATL.Height()+padding < paddedA.Height() )
        {
            int bsize = 
                min( Blocksize(), paddedA.Height()-(ATL.Height()+padding) );
            RepartitionDownDiagonal
            ( ATL,       /**/ paddedATR,  A00,       /**/ A01,       paddedA02,
             /*************************/ /************************************/
                         /**/             A10,       /**/ A11,       paddedA12,
              paddedABL, /**/ paddedABR,  paddedA20, /**/ paddedA21, paddedA22,
              bsize );

            ABR.View
            ( paddedABR, 0, 0,
              paddedABR.Height()-padding, paddedABR.Width()-padding );
            A21.View
            ( paddedA21, 0, 0, 
              paddedA21.Height()-padding, paddedA21.Width() );
            A22.View
            ( paddedA22, 0, 0, 
              paddedA22.Height()-padding, paddedA22.Width()-padding );

            if( A22.Height() > 0 )
            {
                A11Expanded.View( ABR, 0, 0, A11.Height()+1, A11.Width()+1 );
                WPan.AlignWith( A11 );
                WPan.ResizeTo( ABR.Height(), A11.Width() );
                PartitionDown
                ( WPan, W11,
                        W21, A11.Height() );
                e1.AlignWithDiag( ABR, -1 );
                e1.ResizeTo( WPan.Width(), 1 );
                A21_MC_Star.AlignWith( A22 );
                A21_MR_Star.AlignWith( A22 );
                W21_MC_Star.AlignWith( A22 );
                W21_MR_Star.AlignWith( A22 );
                //------------------------------------------------------------//
                lapack::internal::PanelTridiagLSquare
                ( paddedABR, WPan, e1 );

                // Perform the Syr2k for square process grids on A22
                A21_MC_Star = A21; // Allgather row
                A21_MR_Star = A21; // Pairwise exchange then Allgather column
                W21_MC_Star = W21; // Allgather row
                W21_MR_Star = W21; // Pairwise exchange then Allgather column
                blas::internal::LocalTriangularRank2K
                ( Lower, Transpose, Transpose, (R)-1, 
                  A21_MC_Star, W21_MC_Star, A21_MR_Star, W21_MR_Star, (R)1, 
                  A22 );

                A11Expanded.SetDiagonal( e1, -1 );
                //------------------------------------------------------------//
                WPan.FreeAlignments();
                e1.FreeAlignments();
                A21_MC_Star.FreeAlignments();
                A21_MR_Star.FreeAlignments();
                W21_MC_Star.FreeAlignments();
                W21_MR_Star.FreeAlignments();
            }
            else
            {
                A11_Star_Star = A11;
                lapack::Tridiag( Lower, A11_Star_Star.LocalMatrix() );
                A11 = A11_Star_Star;
            }

            SlidePartitionDownDiagonal
            ( ATL,       /**/ paddedATR,  A00,       A01,       /**/ paddedA02,
                         /**/             A10,       A11,       /**/ paddedA12,
             /*************************/ /************************************/
              paddedABL, /**/ paddedABR,  paddedA20, paddedA21, /**/ paddedA22);
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::TridiagLSquare
( DistMatrix<complex<R>,MC,  MR  >& paddedA,
  DistMatrix<complex<R>,Star,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::TridiagLSquare");
    if( paddedA.Height() != paddedA.Width() )
        throw logic_error( "A must be square." );
#endif
    const Grid& g = paddedA.Grid();
    const int padding = g.Height();
    typedef complex<R> C;

    // Separate off the padded parts of A and ensure that they are zero
    DistMatrix<C,MC,MR> 
        A(g),         paddedATR(g),
        paddedABL(g), paddedABR(g);
    PartitionUpLeftDiagonal
    ( paddedA, A,         paddedATR,
               paddedABL, paddedABR, padding );
    paddedATR.SetToZero();
    paddedABL.SetToZero();
    paddedABR.SetToZero();

#ifndef RELEASE
    if( A.Grid() != t.Grid() )
        throw logic_error("A and t must be distributed over the same grid.");
    if( t.Viewing() )
        throw logic_error("t must not be a view.");
    if( g.Height() != g.Width() ||
        A.ColAlignment() != A.RowAlignment() )
        throw logic_error("Square Tridiag is for square, diag aligned grids.");
#endif
    typedef complex<R> C;

    DistMatrix<C,MD,Star> tDiag(g);
    tDiag.AlignWithDiag( A, -1 );
    tDiag.ResizeTo( A.Height()-1, 1 );

    if( g.InGrid() )
    {
        // Matrix views 
        DistMatrix<C,MC,MR> 
            ATL(g), A00(g),       A01(g),       paddedA02(g), 
                    A10(g),       A11(g),       paddedA12(g),
                    paddedA20(g), paddedA21(g), paddedA22(g),
            ABR(g), A11Expanded(g), A21(g), A22(g);
        DistMatrix<C,MD,Star> tT(g),  t0(g), 
                              tB(g),  t1(g),
                                      t2(g);

        // Temporary distributions
        DistMatrix<C,MC,  Star> A21_MC_Star(g);
        DistMatrix<C,MR,  Star> A21_MR_Star(g);
        DistMatrix<C,MC,  Star> W21_MC_Star(g);
        DistMatrix<C,MR,  Star> W21_MR_Star(g);
        DistMatrix<C,Star,Star> A11_Star_Star(g);
        DistMatrix<R,MD,  Star> e1(g);
        DistMatrix<C,Star,Star> t1_Star_Star(g);
        DistMatrix<C,MC,  MR  > W11(g),  WPan(g),
                                W21(g);

        PartitionDownDiagonal
        ( paddedA, ATL,       paddedATR,
                   paddedABL, paddedABR, 0 );
        PartitionDown
        ( tDiag, tT,
                 tB, 0 );
        while( ATL.Height()+padding < paddedA.Height() )
        {
            int bsize = 
                min( Blocksize(), paddedA.Height()-(ATL.Height()+padding) );
            RepartitionDownDiagonal
            ( ATL,       /**/ paddedATR,  A00,       /**/ A01,       paddedA02,
             /*************************/ /************************************/
                         /**/             A10,       /**/ A11,       paddedA12,
              paddedABL, /**/ paddedABR,  paddedA20, /**/ paddedA21, paddedA22,
              bsize );

            RepartitionDown
            ( tT,  t0,
             /**/ /**/
                   t1,
              tB,  t2 );

            ABR.View
            ( paddedABR, 0, 0,
              paddedABR.Height()-padding, paddedABR.Width()-padding );
            A21.View
            ( paddedA21, 0, 0, 
              paddedA21.Height()-padding, paddedA21.Width() );
            A22.View
            ( paddedA22, 0, 0, 
              paddedA22.Height()-padding, paddedA22.Width()-padding );
            
            if( A22.Height() > 0 )
            {
                A11Expanded.View( ABR, 0, 0, A11.Height()+1, A11.Width()+1 );
                WPan.AlignWith( A11 );
                WPan.ResizeTo( ABR.Height(), A11.Width() );
                PartitionDown
                ( WPan, W11,
                        W21, A11.Height() );
                e1.AlignWithDiag( ABR, -1 );
                e1.ResizeTo( WPan.Width(), 1 );
                A21_MC_Star.AlignWith( A22 );
                A21_MR_Star.AlignWith( A22 );
                W21_MC_Star.AlignWith( A22 );
                W21_MR_Star.AlignWith( A22 );
                //------------------------------------------------------------//
                lapack::internal::PanelTridiagLSquare
                ( paddedABR, WPan, e1, t1 );

                // Perform the Her2k for square process grids on A22
                A21_MC_Star = A21; // Allgather row
                A21_MR_Star = A21; // Pairwise exchange then Allgather column
                W21_MC_Star = W21; // Allgather row
                W21_MR_Star = W21; // Pairwise exchange then Allgather column
                blas::internal::LocalTriangularRank2K
                ( Lower, ConjugateTranspose, ConjugateTranspose, (C)-1,
                  A21_MC_Star, W21_MC_Star, A21_MR_Star, W21_MR_Star, (C)1, 
                  A22 );

                A11Expanded.SetDiagonal( e1, -1 );
                //-------------------------------------------------------------//
                WPan.FreeAlignments();
                e1.FreeAlignments();
                A21_MC_Star.FreeAlignments();
                A21_MR_Star.FreeAlignments();
                W21_MC_Star.FreeAlignments();
                W21_MR_Star.FreeAlignments();
            }
            else
            {
                A11_Star_Star = A11;
                t1_Star_Star.ResizeTo( t1.Height(), 1 );

                lapack::Tridiag
                ( Lower, A11_Star_Star.LocalMatrix(), 
                  t1_Star_Star.LocalMatrix() );

                A11 = A11_Star_Star;
                t1 = t1_Star_Star;
            }

            SlidePartitionDownDiagonal
            ( ATL,       /**/ paddedATR,  A00,       A01,       /**/ paddedA02,
                         /**/             A10,       A11,       /**/ paddedA12,
             /*************************/ /************************************/
              paddedABL, /**/ paddedABR,  paddedA20, paddedA21, /**/ paddedA22);

            SlidePartitionDown
            ( tT,  t0,
                   t1,
             /**/ /**/
              tB,  t2 );
        }
    } 

    // Redistribute from matrix-diagonal form to fully replicated
    t = tDiag;
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::internal::TridiagLSquare
( DistMatrix<float,MC,MR>& paddedA );

template void elemental::lapack::internal::TridiagLSquare
( DistMatrix<double,MC,MR>& paddedA );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::TridiagLSquare
( DistMatrix<scomplex,MC,  MR  >& paddedA, 
  DistMatrix<scomplex,Star,Star>& t );

template void elemental::lapack::internal::TridiagLSquare
( DistMatrix<dcomplex,MC,  MR  >& paddedA, 
  DistMatrix<dcomplex,Star,Star>& t );
#endif

