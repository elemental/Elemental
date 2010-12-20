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
using namespace elemental::wrappers::mpi;

template<typename R>
void
elemental::lapack::internal::PanelTridiagLSquare
( DistMatrix<R,MC,MR  >& paddedA,
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MD,Star>& e,
  int padding )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelTridiagLSquare");
#endif 
    const Grid& g = paddedA.Grid();

    // Separate out the padded and original A (the padded parts should be zero)
    DistMatrix<R,MC,MR> 
        A(g),         paddedATR(g),
        paddedABL(g), paddedABR(g);
    PartitionUpLeftDiagonal
    ( paddedA, A,         paddedATR,
               paddedABL, paddedABR, padding );

#ifndef RELEASE
    if( A.Grid() != W.Grid() || W.Grid() != e.Grid() )
        throw logic_error
        ( "A, d, and e must be distributed over the same grid." );
    if( paddedA.Height() != paddedA.Width() )
        throw logic_error( "A must be square." );
    if( A.Height() != W.Height() )
        throw logic_error( "A and W must be the same height." );
    if( W.Height() < W.Width() )
        throw logic_error( "W must be a column panel." );
    if( W.ColAlignment() != A.ColAlignment() || 
        W.RowAlignment() != A.RowAlignment() )
        throw logic_error( "W and A must be aligned." );
    if( e.Height() != W.Width() || e.Width() != 1 )
        throw logic_error
        ( "e must be a column vector of the same length as W's width." );
    if( !e.AlignedWithDiag( A, -1 ) )
        throw logic_error( "e is not aligned with A." );
    if( g.Height() != g.Width() || 
        A.ColAlignment() != A.RowAlignment() )
        throw logic_error("Square Tridiag is for square, diag aligned grids.");
    if( paddedA.Height() % g.Height() != 0 )
        throw logic_error("Square Tridiag requires a padded matrix.");
#endif

    DistMatrix<R,MC,MR> 
        ATL(g), A00(g),       a01(g),       paddedA02(g),  
                a10(g),       alpha11(g),   paddeda12(g),  
                paddedA20(g), paddeda21(g), paddedA22(g),
        ABL(g), A20(g), a21(g), A22(g), ACol(g), alpha21T(g), a21B(g);
    DistMatrix<R,MC,MR> 
        WTL(g), WTR(g),  W00(g), w01(g),     W02(g),  WCol(g),
        WBL(g), WBR(g),  w10(g), omega11(g), w12(g),
                         W20(g), w21(g),     W22(g);
    DistMatrix<R,MD,Star> 
        eT(g),  e0(g),
        eB(g),  epsilon1(g),
                e2(g);

    // Temporary distributions
    DistMatrix<R,MC,  Star> a21_MC_Star(g);
    DistMatrix<R,MC,  Star> paddeda21_MC_Star(g);
    DistMatrix<R,MR,  Star> paddeda21_MR_Star(g);
    DistMatrix<R,Star,MR  > z10_Star_MR(g);
    DistMatrix<R,MC,  MR  > z21(g);
    DistMatrix<R,MC,  Star> z21_MC_Star(g);
    DistMatrix<R,MC,  Star> paddedz21_MC_Star(g);
    DistMatrix<R,MR,  Star> z21_MR_Star(g);
    DistMatrix<R,MR,  Star> paddedz21_MR_Star(g);
    DistMatrix<R,MR,  MC  > z21_MR_MC(g);
    vector<R> u;

    // Push to the blocksize of 1, then pop at the end of the routine
    PushBlocksizeStack( 1 );

    PartitionDownLeftDiagonal
    ( paddedA, ATL,       paddedATR,
               paddedABL, paddedABR, 0 );
    PartitionDownLeftDiagonal
    ( W, WTL, WTR,
         WBL, WBR, 0 );
    PartitionDown
    ( e,  eT,
          eB, 0 );
    while( WTL.Width() < W.Width() )
    {
        RepartitionDownDiagonal
        ( ATL,       /**/ paddedATR,  A00,       /**/ a01,       paddedA02,
         /*************************/ /*************************************/
                     /**/             a10,       /**/ alpha11,   paddeda12, 
          paddedABL, /**/ paddedABR,  paddedA20, /**/ paddeda21, paddedA22 );
        
        RepartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, /**/ w01,     W02,
         /*************/ /**********************/
               /**/       w10, /**/ omega11, w12,
          WBL, /**/ WBR,  W20, /**/ w21,     W22 );

        RepartitionDown
        ( eT,  e0,
         /**/ /********/
               epsilon1,
          eB,  e2 );

        ABL.View
        ( paddedABL, 0, 0, ABL.Height()-padding, ABL.Width() );
        A20.View
        ( paddedA20, 0, 0, A20.Height()-padding, A20.Width() );
        a21.View
        ( paddeda21, 0, 0, a21.Height()-padding, a21.Width() );
        A22.View
        ( paddedA22, 0, 0, A22.Height()-padding, A22.Width()-padding );

        ACol.View2x1
        ( alpha11,
          a21 );

        WCol.View2x1
        ( omega11,
          w21 );
            
        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        paddeda21_MC_Star.AlignWith( paddedA22 );
        paddeda21_MR_Star.AlignWith( paddedA22 );
        z10_Star_MR.AlignWith( W20 );
        z21.AlignWith( w21 );
        paddedz21_MC_Star.AlignWith( paddedA22 );
        paddedz21_MR_Star.AlignWith( paddedA22 );
        paddedz21_MC_Star.ResizeTo( paddeda21.Height(), 1 );
        paddedz21_MR_Star.ResizeTo( paddeda21.Height(), 1 );
        z21_MR_MC.AlignColsWith( A22 );
        z10_Star_MR.ResizeTo( 1, w10.Width() );
        //--------------------------------------------------------------------//
        blas::Gemv( Normal, (R)-1, ABL, w10, (R)1, ACol );
        blas::Gemv( Normal, (R)-1, WBL, a10, (R)1, ACol );

        R tau = 0; // Initializing avoids false compiler warnings
        const bool thisIsMyColumn = ( g.MRRank() == a21.RowAlignment() );
        if( thisIsMyColumn )
            tau = lapack::internal::ColReflector( alpha21T, a21B );
            
        alpha21T.GetDiagonal( epsilon1 );
        alpha21T.Set( 0, 0, (R)1 );

        paddeda21_MC_Star = paddeda21; // Braodcast in row
        paddeda21_MR_Star = paddeda21; // pairwise exchange, Bcast in column

        // Perform the local portion of the Symv
        paddedz21_MC_Star.LocalMatrix() = paddeda21_MC_Star.LocalMatrix();
        if( g.MCRank() < g.MRRank() )
        {
            // We're above the diagonal, so we must zero our local diagonal
            u.resize( paddedA22.LocalWidth() );
            for( int j=0; j<paddedA22.LocalWidth(); ++j )
            {
                u[j] = paddedA22.GetLocalEntry(j,j);
                paddedA22.SetLocalEntry(j,j,0);
            }
            blas::Trmv
            ( Lower, Normal, NonUnit, 
              paddedA22.LockedLocalMatrix(), paddedz21_MC_Star.LocalMatrix() );
            for( int j=0; j<paddedA22.LocalWidth(); ++j )
                paddedA22.SetLocalEntry(j,j,u[j]);
        }
        else
        {
            blas::Trmv
            ( Lower, Normal, NonUnit, 
              paddedA22.LockedLocalMatrix(), paddedz21_MC_Star.LocalMatrix() );
        }
        paddedz21_MR_Star.LocalMatrix() = paddeda21_MR_Star.LocalMatrix();
        if( g.MCRank() <= g.MRRank() ) 
        {
            // We're on or above the diagonal, zero our local diagonal
            u.resize( paddedA22.LocalWidth() );
            for( int j=0; j<paddedA22.LocalWidth(); ++j )
            {
                u[j] = paddedA22.GetLocalEntry(j,j);
                paddedA22.SetLocalEntry(j,j,0);
            }
            blas::Trmv
            ( Lower, Transpose, NonUnit,
              paddedA22.LockedLocalMatrix(), paddedz21_MR_Star.LocalMatrix() );
            for( int j=0; j<paddedA22.LocalWidth(); ++j )
                paddedA22.SetLocalEntry(j,j,u[j]);
        }
        else
        {
            blas::Trmv
            ( Lower, Transpose, NonUnit,
              paddedA22.LockedLocalMatrix(), paddedz21_MR_Star.LocalMatrix() );
        }

        a21_MC_Star.View
        ( paddeda21_MC_Star, 0, 0, a21_MC_Star.Height()-padding, 1 );
        z21_MC_Star.View
        ( paddedz21_MC_Star, 0, 0, z21_MC_Star.Height()-padding, 1 );
        z21_MR_Star.View
        ( paddedz21_MR_Star, 0, 0, z21_MR_Star.Height()-padding, 1 );

        blas::Gemv
        ( Transpose, 
          (R)1, W20.LockedLocalMatrix(),
                a21_MC_Star.LockedLocalMatrix(),
          (R)0, z10_Star_MR.LocalMatrix() );
        z10_Star_MR.SumOverCol();

        blas::Gemv
        ( Normal, 
          (R)-1, A20.LockedLocalMatrix(),
                 z10_Star_MR.LockedLocalMatrix(),
          (R)+1, z21_MC_Star.LocalMatrix() );

        blas::Gemv
        ( Transpose, 
          (R)1, A20.LockedLocalMatrix(),
                a21_MC_Star.LockedLocalMatrix(),
          (R)0, z10_Star_MR.LocalMatrix() );
        z10_Star_MR.SumOverCol();

        blas::Gemv
        ( Normal, 
          (R)-1, W20.LockedLocalMatrix(),
                 z10_Star_MR.LockedLocalMatrix(),
          (R)+1, z21_MC_Star.LocalMatrix() );

        w21.SumScatterFrom( z21_MC_Star ); // SumScatter col
        z21_MR_MC.SumScatterFrom( z21_MR_Star ); // SumScatter row
        z21 = z21_MR_MC; // pairwise exchange

        if( thisIsMyColumn )
        {
            blas::Axpy( (R)1, z21, w21 );
            blas::Scal( tau, w21 );

            R alpha;
            R myAlpha = -static_cast<R>(0.5)*tau*
                        blas::Dot( w21.LockedLocalMatrix(), 
                                   a21.LockedLocalMatrix() );            
            AllReduce( &myAlpha, &alpha, 1, MPI_SUM, g.MCComm() );
            blas::Axpy( alpha, a21, w21 );
        }
        //--------------------------------------------------------------------//
        paddeda21_MC_Star.FreeAlignments();
        paddeda21_MR_Star.FreeAlignments();
        z10_Star_MR.FreeAlignments();
        z21.FreeAlignments();
        paddedz21_MC_Star.FreeAlignments();
        paddedz21_MR_Star.FreeAlignments();
        z21_MR_MC.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL,       /**/ paddedATR,  A00,       a01,       /**/ paddedA02,
                     /**/             a10,       alpha11,   /**/ paddeda12,
         /*************************/ /************************************/
          paddedABL, /**/ paddedABR,  paddedA20, paddeda21, /**/ paddedA22 );

        SlidePartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );
        
        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::PanelTridiagLSquare
( DistMatrix<complex<R>,MC,MR  >& paddedA,
  DistMatrix<complex<R>,MC,MR  >& W,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<complex<R>,MD,Star>& t,
  int padding )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelTridiagLSquare");
#endif    
    const Grid& g = paddedA.Grid();
    typedef complex<R> C;

    // Separate out the padded and original A (the padded parts should be zero)
    DistMatrix<C,MC,MR>
        A(g),         paddedATR(g),
        paddedABL(g), paddedABR(g);
    PartitionUpLeftDiagonal
    ( paddedA, A,         paddedATR,
               paddedABL, paddedABR, padding );

#ifndef RELEASE
    if( A.Grid() != W.Grid() || 
        W.Grid() != e.Grid() || 
        e.Grid() != t.Grid() )
        throw logic_error
        ( "A, d, e, and t must be distributed over the same grid." );
    if( A.Height() != A.Width() )
        throw logic_error( "A must be square." );
    if( A.Height() != W.Height() )
        throw logic_error( "A and W must be the same height." );
    if( W.Height() < W.Width() )
        throw logic_error( "W must be a column panel." );
    if( W.ColAlignment() != A.ColAlignment() || 
        W.RowAlignment() != A.RowAlignment() )
        throw logic_error( "W and A must be aligned." );
    if( e.Height() != W.Width() || e.Width() != 1 )
        throw logic_error
        ( "e must be a column vector of the same length as W's width." );
    if( !e.AlignedWithDiag( A, -1 ) )
        throw logic_error( "e is not aligned with A's subdiagonal." );
    if( t.Height() != W.Width() || t.Width() != 1 )
        throw logic_error
        ( "t must be a column vector of the same length as W's width." );
    if( !t.AlignedWithDiag( A, -1 ) )
        throw logic_error( "t is not aligned with A's subdiagonal." );
    if( g.Height() != g.Width() || 
        A.ColAlignment() != A.RowAlignment() )
        throw logic_error("Square Tridiag is for square, diag aligned grids.");
    if( paddedA.Height() % g.Height() != 0 )
        throw logic_error("Square Tridiag requires a padded matrix.");
#endif

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ATL(g),  A00(g),       a01(g),       paddedA02(g), 
                 a10(g),       alpha11(g),   paddeda12(g),
                 paddedA20(g), paddeda21(g), paddedA22(g),
        ABL(g), A20(g), a21(g), A22(g), ACol(g), alpha21T(g), a21B(g); 
    DistMatrix<C,MC,MR> 
        WTL(g), WTR(g),  W00(g), w01(g),     W02(g),  WCol(g),
        WBL(g), WBR(g),  w10(g), omega11(g), w12(g),
                         W20(g), w21(g),     W22(g);
    DistMatrix<R,MD,Star> 
        eT(g),  e0(g),
        eB(g),  epsilon1(g),
                e2(g);
    DistMatrix<C,MD,Star>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    DistMatrix<C,MC,  MR  > a10Conj(g);
    DistMatrix<C,MC,  MR  > w10Conj(g);
    DistMatrix<C,MC,  Star> a21_MC_Star(g);
    DistMatrix<C,MC,  Star> paddeda21_MC_Star(g);
    DistMatrix<C,MR,  Star> paddeda21_MR_Star(g);
    DistMatrix<C,Star,MR  > z10_Star_MR(g);
    DistMatrix<C,MC,  MR  > z21(g);
    DistMatrix<C,MC,  Star> z21_MC_Star(g);
    DistMatrix<C,MC,  Star> paddedz21_MC_Star(g);
    DistMatrix<C,MR,  Star> z21_MR_Star(g);
    DistMatrix<C,MR,  Star> paddedz21_MR_Star(g);
    DistMatrix<C,MR,  MC  > z21_MR_MC(g);
    vector<C> u;

    // Push to the blocksize of 1, then pop at the end of the routine
    PushBlocksizeStack( 1 );

    PartitionDownLeftDiagonal
    ( paddedA, ATL,       paddedATR,
               paddedABL, paddedABR, 0 );
    PartitionDownLeftDiagonal
    ( W, WTL, WTR,
         WBL, WBR, 0 );
    PartitionDown
    ( e,  eT,
          eB, 0 );
    PartitionDown
    ( t,  tT,
          tB, 0 );
    while( WTL.Width() < W.Width() )
    {
        RepartitionDownDiagonal
        ( ATL,       /**/ paddedATR,  A00,       /**/ a01,       paddedA02,
         /*************************/ /************************************/
                     /**/             a10,       /**/ alpha11,   paddeda12, 
          paddedABL, /**/ paddedABR,  paddedA20, /**/ paddeda21, paddedA22 );
        
        RepartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, /**/ w01,     W02,
         /*************/ /**********************/
               /**/       w10, /**/ omega11, w12,
          WBL, /**/ WBR,  W20, /**/ w21,     W22 );

        RepartitionDown
        ( eT,  e0,
         /**/ /********/
               epsilon1,
          eB,  e2 );

        RepartitionDown
        ( tT,  t0,
         /**/ /****/
               tau1,
          tB,  t2 );

        ABL.View
        ( paddedABL, 0, 0, ABL.Height()-padding, ABL.Width() );
        A20.View
        ( paddedA20, 0, 0, A20.Height()-padding, A20.Width() );
        a21.View
        ( paddeda21, 0, 0, a21.Height()-padding, a21.Width() );
        A22.View
        ( paddedA22, 0, 0, A22.Height()-padding, A22.Width()-padding );

        ACol.View2x1
        ( alpha11,
          a21 );

        WCol.View2x1
        ( omega11,
          w21 );
        
        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        paddeda21_MC_Star.AlignWith( paddedA22 );
        paddeda21_MR_Star.AlignWith( paddedA22 );
        z10_Star_MR.AlignWith( W20 );
        z21.AlignWith( w21 );
        paddedz21_MC_Star.AlignWith( paddedA22 );
        paddedz21_MR_Star.AlignWith( paddedA22 );
        paddedz21_MC_Star.ResizeTo( paddeda21.Height(), 1 );
        paddedz21_MR_Star.ResizeTo( paddeda21.Height(), 1 );
        z21_MR_MC.AlignColsWith( A22 );
        z10_Star_MR.ResizeTo( 1, w10.Width() );
        //--------------------------------------------------------------------//
        alpha11.SetImag( 0, 0, (R)0 ); 
        blas::Conj( w10, w10Conj );
        blas::Gemv( Normal, (C)-1, ABL, w10Conj, (C)1, ACol );
        blas::Conj( a10, a10Conj );
        blas::Gemv( Normal, (C)-1, WBL, a10Conj, (C)1, ACol );
        alpha11.SetImag( 0, 0, (R)0 );

        C tau = 0; // Initializing avoids false compiler warnings
        const bool thisIsMyColumn = ( g.MRRank() == a21.RowAlignment() );
        if( thisIsMyColumn )
        {
            tau = lapack::internal::ColReflector( alpha21T, a21B );
            const bool thisIsMyRow = ( g.MCRank() == alpha21T.ColAlignment() );
            if( thisIsMyRow )
                tau1.SetLocalEntry(0,0,tau);
        }
            
        alpha21T.GetRealDiagonal( epsilon1 );
        alpha21T.Set( 0, 0, (C)1 );

        paddeda21_MC_Star = paddeda21; // Broadcast in row
        paddeda21_MR_Star = paddeda21; // pairwise exchange, Bcast in col

        // Perform the local portion of the Hemv
        paddedz21_MC_Star.LocalMatrix() = paddeda21_MC_Star.LocalMatrix();
        if( g.MCRank() < g.MRRank() )
        {
            // We're above the diagonal, so we must zero our local diagonal
            u.resize( paddedA22.LocalWidth() );
            for( int j=0; j<paddedA22.LocalWidth(); ++j )
            {
                u[j] = paddedA22.GetLocalEntry(j,j);
                paddedA22.SetLocalEntry(j,j,0);
            }
            blas::Trmv
            ( Lower, Normal, NonUnit,
              paddedA22.LockedLocalMatrix(), paddedz21_MC_Star.LocalMatrix() );
            for( int j=0; j<paddedA22.LocalWidth(); ++j )
                paddedA22.SetLocalEntry(j,j,u[j]);
        }
        else
        {
            blas::Trmv
            ( Lower, Normal, NonUnit,
              paddedA22.LockedLocalMatrix(), paddedz21_MC_Star.LocalMatrix() );
        }
        paddedz21_MR_Star.LocalMatrix() = paddeda21_MR_Star.LocalMatrix();
        if( g.MCRank() <= g.MRRank() )
        {
            // We're on or above the diagonal, zero our local diagonal
            u.resize( paddedA22.LocalWidth() );
            for( int j=0; j<paddedA22.LocalWidth(); ++j )
            {
                u[j] = paddedA22.GetLocalEntry(j,j);
                paddedA22.SetLocalEntry(j,j,0);
            }
            blas::Trmv
            ( Lower, ConjugateTranspose, NonUnit,
              paddedA22.LockedLocalMatrix(), paddedz21_MR_Star.LocalMatrix() );
            for( int j=0; j<paddedA22.LocalWidth(); ++j )
                paddedA22.SetLocalEntry(j,j,u[j]);
        }
        else
        {
            blas::Trmv
            ( Lower, ConjugateTranspose, NonUnit,
              paddedA22.LockedLocalMatrix(), paddedz21_MR_Star.LocalMatrix() );
        }

        a21_MC_Star.View
        ( paddeda21_MC_Star, 0, 0, a21_MC_Star.Height()-padding, 1 );
        z21_MC_Star.View
        ( paddedz21_MC_Star, 0, 0, z21_MC_Star.Height()-padding, 1 );
        z21_MR_Star.View
        ( paddedz21_MR_Star, 0, 0, z21_MR_Star.Height()-padding, 1 );

        blas::Gemv
        ( ConjugateTranspose, 
          (C)1, W20.LockedLocalMatrix(),
                a21_MC_Star.LockedLocalMatrix(),
          (C)0, z10_Star_MR.LocalMatrix() );
        z10_Star_MR.SumOverCol();

        blas::Gemv
        ( Normal, 
          (C)-1, A20.LockedLocalMatrix(),
                 z10_Star_MR.LockedLocalMatrix(),
          (C)+1, z21_MC_Star.LocalMatrix() );

        blas::Gemv
        ( ConjugateTranspose, 
          (C)1, A20.LockedLocalMatrix(),
                a21_MC_Star.LockedLocalMatrix(),
          (C)0, z10_Star_MR.LocalMatrix() );
        z10_Star_MR.SumOverCol();

        blas::Gemv
        ( Normal, 
          (C)-1, W20.LockedLocalMatrix(),
                 z10_Star_MR.LockedLocalMatrix(),
          (C)+1, z21_MC_Star.LocalMatrix() );

        w21.SumScatterFrom( z21_MC_Star ); // SumScatter col
        z21_MR_MC.SumScatterFrom( z21_MR_Star ); // SumScatter row
        z21 = z21_MR_MC; // pairwise exchange

        if( thisIsMyColumn )
        {
            blas::Axpy( (C)1, z21, w21 );
            blas::Scal( tau, w21 );

            C alpha;
            C myAlpha = -static_cast<C>(0.5)*tau*
                        blas::Dot( w21.LockedLocalMatrix(), 
                                   a21.LockedLocalMatrix() );            
            AllReduce( &myAlpha, &alpha, 1, MPI_SUM, g.MCComm() );
            blas::Axpy( alpha, a21, w21 );
        }
        //--------------------------------------------------------------------//
        paddeda21_MC_Star.FreeAlignments();
        paddeda21_MR_Star.FreeAlignments();
        z10_Star_MR.FreeAlignments();
        z21.FreeAlignments();
        paddedz21_MC_Star.FreeAlignments();
        paddedz21_MR_Star.FreeAlignments();
        z21_MR_MC.FreeAlignments();

        SlidePartitionDownDiagonal
        ( ATL,       /**/ paddedATR,  A00,       a01,       /**/ paddedA02,
                     /**/             a10,       alpha11,   /**/ paddeda12,
         /*************************/ /************************************/
          paddedABL, /**/ paddedABR,  paddedA20, paddeda21, /**/ paddedA22 );

        SlidePartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );
        
        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        SlidePartitionDown
        ( tT,  t0,
         /**/ /****/
               tau1,
          tB,  t2 );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::internal::PanelTridiagLSquare
( DistMatrix<float,MC,MR  >& paddedA,
  DistMatrix<float,MC,MR  >& W,
  DistMatrix<float,MD,Star>& e,
  int padding );

template void elemental::lapack::internal::PanelTridiagLSquare
( DistMatrix<double,MC,MR  >& paddedA,
  DistMatrix<double,MC,MR  >& W,
  DistMatrix<double,MD,Star>& e,
  int padding );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::PanelTridiagLSquare
( DistMatrix<scomplex,MC,MR  >& paddedA,
  DistMatrix<scomplex,MC,MR  >& W,
  DistMatrix<float,   MD,Star>& e,
  DistMatrix<scomplex,MD,Star>& t,
  int padding );

template void elemental::lapack::internal::PanelTridiagLSquare
( DistMatrix<dcomplex,MC,MR  >& paddedA,
  DistMatrix<dcomplex,MC,MR  >& W,
  DistMatrix<double,  MD,Star>& e,
  DistMatrix<dcomplex,MD,Star>& t,
  int padding );
#endif

