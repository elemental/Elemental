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
using namespace elemental::wrappers::mpi;

template<typename R>
void
elemental::lapack::internal::PanelTridiagL
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<R,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelTridiagL");
    if( A.GetGrid() != W.GetGrid() ||
        W.GetGrid() != e.GetGrid() ||
        e.GetGrid() != t.GetGrid() )
        throw "A, d, e, and t must be distributed over the same grid.";
    if( A.Height() != A.Width() )
        throw "A must be square.";
    if( A.Height() != W.Height() )
        throw "A and W must be the same height.";
    if( W.Height() <= W.Width() )
        throw "W must be a column panel.";
    if( W.ColAlignment() != A.ColAlignment() || 
        W.RowAlignment() != A.RowAlignment() )
        throw "W and A must be aligned.";
    if( e.Height() != W.Width() || e.Width() != 1 )
        throw "e must be a column vector of the same length as W's width.";
    if( t.Height() != W.Width() || t.Width() != 1 )
        throw "t must be a column vector of the same length as W's width.";
    if( e.ColAlignment() != ((A.ColAlignment()+1) % e.GetGrid().Height())
                            + A.RowAlignment() * e.GetGrid().Height() )
        throw "e is not aligned with A.";
    if( t.ColAlignment() != (A.ColAlignment()+
                             A.RowAlignment()*t.GetGrid().Height()) )
        throw "t is not aligned with A.";
#endif
    const Grid& grid = A.GetGrid();

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), a01(grid),     A02(grid),  ACol(grid),
        ABL(grid), ABR(grid),  a10(grid), alpha11(grid), a12(grid),
                               A20(grid), a21(grid),     A22(grid);
    DistMatrix<R,MC,MR> 
        WTL(grid), WTR(grid),  W00(grid), w01(grid),     W02(grid),  WCol(grid),
        WBL(grid), WBR(grid),  w10(grid), omega11(grid), w12(grid),
                               W20(grid), w21(grid),     W22(grid);
    DistMatrix<R,MD,Star> 
        eT(grid),  e0(grid),
        eB(grid),  epsilon1(grid),
                   e2(grid);
    DistMatrix<R,MD,Star> 
        tT(grid),  t0(grid),
        tB(grid),  tau1(grid),
                   t2(grid);

    // Temporary distributions
    DistMatrix<R,MC,  Star> a21_MC_Star(grid);
    DistMatrix<R,MR,  Star> a21_MR_Star(grid);
    DistMatrix<R,Star,MR  > z10_Star_MR(grid);
    DistMatrix<R,MC,  MR  > z21(grid);
    DistMatrix<R,MC,  Star> z21_MC_Star(grid);
    DistMatrix<R,MR,  Star> z21_MR_Star(grid);
    DistMatrix<R,MR,  MC  > z21_MR_MC(grid);

    // Push to the blocksize of 1, then pop at the end of the routine
    PushBlocksizeStack( 1 );

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR );
    PartitionDown
    ( e,  eT,
          eB );
    PartitionDown
    ( t,  tT,
          tB );
    while( WTL.Width() < W.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12, 
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );
        
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

        ACol.View2x1
        ( alpha11,
          a21 );

        WCol.View2x1
        ( omega11,
          w21 );

        a21_MC_Star.AlignWith( A22 );
        a21_MR_Star.AlignWith( A22 );
        z10_Star_MR.AlignWith( W20 );
        z21.AlignWith( w21 );
        z21_MC_Star.AlignWith( A22 );
        z21_MR_Star.AlignWith( A22 );
        z21_MR_MC.AlignColsWith( A22 );
        z21_MC_Star.ResizeTo( w21.Height(), 1 );
        z21_MR_Star.ResizeTo( w21.Height(), 1 );
        z10_Star_MR.ResizeTo( 1, w10.Width() );
        z21_MC_Star.SetToZero();
        z21_MR_Star.SetToZero();
        //--------------------------------------------------------------------//
        blas::Gemv( Normal, (R)-1, ABL, w10, (R)1, ACol );
        blas::Gemv( Normal, (R)-1, WBL, a10, (R)1, ACol );

        R tau = 0; // Initializing avoids false compiler warnings
        const bool thisIsMyColumn = ( grid.MRRank() == a21.RowAlignment() );
        if( thisIsMyColumn )
        {
            tau = lapack::internal::LocalColReflector( a21 );
            tau1.Set( 0, 0, tau );
        }
            
        a21.GetDiagonal( epsilon1 );
        a21.Set( 0, 0, (R)1 ); 

        a21_MR_Star = a21_MC_Star = a21;

        PopBlocksizeStack();
        blas::internal::SymvColAccumulate
        ( Lower, (R)1, A22, a21_MC_Star, a21_MR_Star, 
                            z21_MC_Star, z21_MR_Star );
        PushBlocksizeStack( 1 );

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

        w21.SumScatterFrom( z21_MC_Star );
        z21_MR_MC.SumScatterFrom( z21_MR_Star );
        z21 = z21_MR_MC;

        if( thisIsMyColumn )
        {
            blas::Axpy( (R)1, z21, w21 );
            blas::Scal( (R)tau, w21 );

            R alpha;
            R myAlpha = -static_cast<R>(0.5)*tau*
                        blas::Dot( w21.LockedLocalMatrix(), 
                                   a21.LockedLocalMatrix() );            
            AllReduce( &myAlpha, &alpha, 1, MPI_SUM, grid.MCComm() );
            blas::Axpy( alpha, a21, w21 );
        }
        //--------------------------------------------------------------------//
        a21_MC_Star.FreeAlignments();
        a21_MR_Star.FreeAlignments();
        z10_Star_MR.FreeAlignments();
        z21.FreeAlignments();
        z21_MC_Star.FreeAlignments();
        z21_MR_Star.FreeAlignments();
        z21_MR_MC.FreeAlignments();

        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R>
void
elemental::lapack::internal::PanelTridiagL
( DistMatrix<complex<R>,MC,MR  >& A,
  DistMatrix<complex<R>,MC,MR  >& W,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<complex<R>,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("lapack::internal::PanelTridiagL");
    if( A.GetGrid() != W.GetGrid() ||
        W.GetGrid() != e.GetGrid() ||
        e.GetGrid() != t.GetGrid() )
        throw "A, d, e, and t must be distributed over the same grid.";
    if( A.Height() != A.Width() )
        throw "A must be square.";
    if( A.Height() != W.Height() )
        throw "A and W must be the same height.";
    if( W.Height() <= W.Width() )
        throw "W must be a column panel.";
    if( W.ColAlignment() != A.ColAlignment() || 
        W.RowAlignment() != A.RowAlignment() )
        throw "W and A must be aligned.";
    if( e.Height() != W.Width() || e.Width() != 1 )
        throw "e must be a column vector of the same length as W's width.";
    if( t.Height() != W.Width() || t.Width() != 1 )
        throw "t must be a column vector of the same length as W's width.";
    if( e.ColAlignment() != ((A.ColAlignment()+1) % e.GetGrid().Height())
                            + A.RowAlignment() * e.GetGrid().Height() )
        throw "e is not aligned with A.";
    if( t.ColAlignment() != (A.ColAlignment()+
                             A.RowAlignment()*t.GetGrid().Height()) )
        throw "t is not aligned with A.";
#endif
    typedef complex<R> C;

    const Grid& grid = A.GetGrid();

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), a01(grid),     A02(grid),  ACol(grid),
        ABL(grid), ABR(grid),  a10(grid), alpha11(grid), a12(grid),
                               A20(grid), a21(grid),     A22(grid);
    DistMatrix<C,MC,MR> 
        WTL(grid), WTR(grid),  W00(grid), w01(grid),     W02(grid),  WCol(grid),
        WBL(grid), WBR(grid),  w10(grid), omega11(grid), w12(grid),
                               W20(grid), w21(grid),     W22(grid);
    DistMatrix<R,MD,Star> 
        eT(grid),  e0(grid),
        eB(grid),  epsilon1(grid),
                   e2(grid);
    DistMatrix<C,MD,Star> 
        tT(grid),  t0(grid),
        tB(grid),  tau1(grid),
                   t2(grid);

    // Temporary distributions
    DistMatrix<C,MC,  MR  > a10Conj(grid);
    DistMatrix<C,MC,  MR  > w10Conj(grid);
    DistMatrix<C,MC,  Star> a21_MC_Star(grid);
    DistMatrix<C,MR,  Star> a21_MR_Star(grid);
    DistMatrix<C,Star,MR  > z10_Star_MR(grid);
    DistMatrix<C,MC,  MR  > z21(grid);
    DistMatrix<C,MC,  Star> z21_MC_Star(grid);
    DistMatrix<C,MR,  Star> z21_MR_Star(grid);
    DistMatrix<C,MR,  MC  > z21_MR_MC(grid);

    // Push to the blocksize of 1, then pop at the end of the routine
    PushBlocksizeStack( 1 );

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR );
    PartitionDown
    ( e,  eT,
          eB );
    PartitionDown
    ( t,  tT,
          tB );
    while( WTL.Width() < W.Width() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12, 
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );
        
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

        ACol.View2x1
        ( alpha11,
          a21 );

        WCol.View2x1
        ( omega11,
          w21 );

        a21_MC_Star.AlignWith( A22 );
        a21_MR_Star.AlignWith( A22 );
        z10_Star_MR.AlignWith( W20 );
        z21.AlignWith( w21 );
        z21_MC_Star.AlignWith( A22 );
        z21_MR_Star.AlignWith( A22 );
        z21_MR_MC.AlignColsWith( A22 );
        z21_MC_Star.ResizeTo( w21.Height(), 1 );
        z21_MR_Star.ResizeTo( w21.Height(), 1 );
        z10_Star_MR.ResizeTo( 1, w10.Width() );
        z21_MC_Star.SetToZero();
        z21_MR_Star.SetToZero();
        //--------------------------------------------------------------------//
        alpha11.SetImag( 0, 0, (R)0 ); 
        blas::Conj( w10, w10Conj );
        blas::Gemv( Normal, (C)-1, ABL, w10Conj, (C)1, ACol );
        blas::Conj( a10, a10Conj );
        blas::Gemv( Normal, (C)-1, WBL, a10Conj, (C)1, ACol );
        alpha11.SetImag( 0, 0, (R)0 );

        C tau = 0; // Initializing avoids false compiler warnings
        const bool thisIsMyColumn = ( grid.MRRank() == a21.RowAlignment() );
        if( thisIsMyColumn )
        {
            tau = lapack::internal::LocalColReflector( a21 );
            tau1.Set( 0, 0, tau );
        }
            
        a21.GetRealDiagonal( epsilon1 );
        a21.Set( 0, 0, (C)1 ); 

        a21_MR_Star = a21_MC_Star = a21;

        PopBlocksizeStack();
        blas::internal::HemvColAccumulate
        ( Lower, (C)1, A22, a21_MC_Star, a21_MR_Star, 
                            z21_MC_Star, z21_MR_Star );
        PushBlocksizeStack( 1 );

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

        w21.SumScatterFrom( z21_MC_Star );
        z21_MR_MC.SumScatterFrom( z21_MR_Star );
        z21 = z21_MR_MC;

        if( thisIsMyColumn )
        {
            blas::Axpy( (C)1, z21, w21 );
            blas::Scal( tau, w21 );

            C alpha;
            C myAlpha = -static_cast<C>(0.5)*tau*
                        blas::Dot( w21.LockedLocalMatrix(), 
                                   a21.LockedLocalMatrix() );            
            AllReduce( &myAlpha, &alpha, 1, MPI_SUM, grid.MCComm() );
            blas::Axpy( alpha, a21, w21 );
        }
        //--------------------------------------------------------------------//
        a21_MC_Star.FreeAlignments();
        a21_MR_Star.FreeAlignments();
        z10_Star_MR.FreeAlignments();
        z21.FreeAlignments();
        z21_MC_Star.FreeAlignments();
        z21_MR_Star.FreeAlignments();
        z21_MR_MC.FreeAlignments();

        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

template void elemental::lapack::internal::PanelTridiagL
( DistMatrix<float,MC,MR  >& A,
  DistMatrix<float,MC,MR  >& W,
  DistMatrix<float,MD,Star>& e,
  DistMatrix<float,MD,Star>& t );

template void elemental::lapack::internal::PanelTridiagL
( DistMatrix<double,MC,MR  >& A,
  DistMatrix<double,MC,MR  >& W,
  DistMatrix<double,MD,Star>& e,
  DistMatrix<double,MD,Star>& t );

#ifndef WITHOUT_COMPLEX
template void elemental::lapack::internal::PanelTridiagL
( DistMatrix<scomplex,MC,MR  >& A,
  DistMatrix<scomplex,MC,MR  >& W,
  DistMatrix<float,   MD,Star>& e,
  DistMatrix<scomplex,MD,Star>& t );

template void elemental::lapack::internal::PanelTridiagL
( DistMatrix<dcomplex,MC,MR  >& A,
  DistMatrix<dcomplex,MC,MR  >& W,
  DistMatrix<double,  MD,Star>& e,
  DistMatrix<dcomplex,MD,Star>& t );
#endif

