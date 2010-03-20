/*
   Copyright 2009-2010 Jack Poulson

   This file is part of Elemental.

   Elemental is free software: you can redistribute it and/or modify it under
   the terms of the GNU Lesser General Public License as published by the
   Free Software Foundation; either version 3 of the License, or 
   (at your option) any later version.

   Elemental is distributed in the hope that it will be useful, but 
   WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public License
   along with Elemental. If not, see <http://www.gnu.org/licenses/>.
*/
#include "ElementalBLAS_Internal.h"
#include "ElementalLAPACK_Internal.h"
using namespace std;
using namespace Elemental;
using namespace Elemental::wrappers::MPI;

template<typename R>
void
Elemental::LAPACK::Internal::PanelTridiagL
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MD,Star>& e,
  DistMatrix<R,MD,Star>& t )
{
#ifndef RELEASE
    PushCallStack("LAPACK::Internal::PanelTridiagL");
#endif
    const Grid& grid = A.GetGrid();
#ifndef RELEASE
    if( A.GetGrid() != W.GetGrid() ||
        W.GetGrid() != e.GetGrid() ||
        e.GetGrid() != t.GetGrid()   )
    {
        if( grid.VCRank() == 0 )
            cerr << "A, d, e, and t must be distributed over the same grid."
                 << endl;
        DumpCallStack();
        throw exception();
    }
    if( A.Height() != A.Width() )
    {
        if( grid.VCRank() == 0 )
            cerr << "A must be square." << endl;
        DumpCallStack();
        throw exception();
    }
    if( A.Height() != W.Height() )
    {
        if( grid.VCRank() == 0 )
            cerr << "A and W must be the same height." << endl;
        DumpCallStack();
        throw exception();
    }
    if( W.Height() <= W.Width() )
    {
        if( grid.VCRank() == 0 )
            cerr << "W must be a column panel." << endl;
        DumpCallStack();
        throw exception();
    }
    if( W.ColAlignment() != A.ColAlignment() || 
        W.RowAlignment() != A.RowAlignment()   )
    {
        if( grid.VCRank() == 0 )
            cerr << "W and A must be aligned." << endl;
        DumpCallStack();
        throw exception();
    }
    if( e.Height() != W.Width() || e.Width() != 1 )
    {
        if( grid.VCRank() == 0 )
            cerr << "e must be a column vector of the same length as W's width."
                 << endl;
        DumpCallStack();
        throw exception();
    }
    if( t.Height() != W.Width() || t.Width() != 1 )
    {
        if( grid.VCRank() == 0 )
            cerr << "t must be a column vector of the same length as W's width."
                 << endl;
        DumpCallStack();
        throw exception();
    }
    if( e.ColAlignment() != ((A.ColAlignment()+1) % grid.Height())
                            + A.RowAlignment() * grid.Height()    )
    {
        if( grid.VCRank() == 0 )
            cerr << "e is not aligned with A." << endl;
        DumpCallStack();
        throw exception();
    }
    if( t.ColAlignment() != (A.ColAlignment()+A.RowAlignment()*grid.Height()) )
    {
        if( grid.VCRank() == 0 )
            cerr << "t is not aligned with A." << endl;
        DumpCallStack();
        throw exception();
    }
#endif
    const int myCol  = grid.MRRank();
    const int myRank = grid.VCRank();

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ATL(grid), ATR(grid),  A00(grid), a01(grid),     A02(grid),  ACol(grid),
        ABL(grid), ABR(grid),  a10(grid), alpha11(grid), a12(grid),
                               A20(grid), a21(grid),     A22(grid);
    DistMatrix<R,MC,MR> 
        WTL(grid), WTR(grid),  W00(grid), w01(grid),     W02(grid),  WCol(grid),
        WBL(grid), WBR(grid),  w10(grid), omega11(grid), w12(grid),
                               W20(grid), w21(grid),     W22(grid);
    DistMatrix<R,MD,Star> eT(grid),  e0(grid),
                          eB(grid),  epsilon1(grid),
                                                  e2(grid);
    DistMatrix<R,MD,Star> tT(grid),  t0(grid),
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

    PartitionDownDiagonal( A, ATL, ATR,
                              ABL, ABR );
    PartitionDownDiagonal( W, WTL, WTR,
                              WBL, WBR );
    PartitionDown( e,  eT,
                       eB );
    PartitionDown( t,  tT,
                       tB );
    while( WTL.Width() < W.Width() )
    {
        RepartitionDownDiagonal( ATL, /**/ ATR,  A00, /**/ a01,     A02,
                                /*************/ /**********************/
                                      /**/       a10, /**/ alpha11, a12, 
                                 ABL, /**/ ABR,  A20, /**/ a21,     A22 );
        
        RepartitionDownDiagonal( WTL, /**/ WTR,  W00, /**/ w01,     W02,
                                /*************/ /**********************/
                                      /**/       w10, /**/ omega11, w12,
                                 WBL, /**/ WBR,  W20, /**/ w21,     W22 );

        RepartitionDown( eT,  e0,
                        /**/ /********/
                              epsilon1,
                         eB,  e2       );

        RepartitionDown( tT,  t0,
                        /**/ /****/
                              tau1,
                         tB,  t2   );

        ACol.View2x1( alpha11,
                      a21     );

        WCol.View2x1( omega11,
                      w21     );

        a21_MC_Star.ConformWith( A22 );
        a21_MR_Star.ConformWith( A22 );
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
        BLAS::Gemv( Normal, (R)-1, ABL, w10, (R)1, ACol );
        BLAS::Gemv( Normal, (R)-1, WBL, a10, (R)1, ACol );

        R tau = 0; // Initializing avoids false compiler warnings
        const bool thisIsMyColumn = ( myCol == a21.RowAlignment() );
        if( thisIsMyColumn )
        {
            tau = LAPACK::Internal::LocalColReflector( a21 );
            if( myRank == tau1.ColAlignment() )
                tau1.LocalEntry(0,0) = tau;
        }
            
        a21.GetDiagonal( epsilon1 );
        a21.Set( 0, 0, (R)1 ); 

        // Set up for the W updates
        a21_MC_Star = a21;
        a21_MR_Star = a21_MC_Star;

        PopBlocksizeStack();
        BLAS::Internal::SymvColAccumulate
        ( Lower, (R)1, A22, a21_MC_Star, a21_MR_Star, 
                            z21_MC_Star, z21_MR_Star );
        PushBlocksizeStack( 1 );

        BLAS::Gemv
        ( Transpose, (R)1, W20.LockedLocalMatrix(),
                           a21_MC_Star.LockedLocalMatrix(),
                     (R)0, z10_Star_MR.LocalMatrix()       );
        z10_Star_MR.AllReduce();
        BLAS::Gemv
        ( Normal, (R)-1, A20.LockedLocalMatrix(),
                         z10_Star_MR.LockedLocalMatrix(),
                  (R)+1, z21_MC_Star.LocalMatrix()       );
        BLAS::Gemv
        ( Transpose, (R)1, A20.LockedLocalMatrix(),
                           a21_MC_Star.LockedLocalMatrix(),
                     (R)0, z10_Star_MR.LocalMatrix()       );
        z10_Star_MR.AllReduce();
        BLAS::Gemv
        ( Normal, (R)-1, W20.LockedLocalMatrix(),
                         z10_Star_MR.LockedLocalMatrix(),
                  (R)+1, z21_MC_Star.LocalMatrix()       );

        w21.ReduceScatterUpdate( (R)1, z21_MC_Star );
        z21_MR_MC.ReduceScatterFrom( z21_MR_Star );
        z21 = z21_MR_MC;

        if( thisIsMyColumn )
        {
            BLAS::Axpy( (R)1, z21, w21 );
            BLAS::Scal( tau, w21 );

            R alpha;
            R myAlpha = -static_cast<R>(0.5)*tau*
                        BLAS::Dot( w21.LockedLocalMatrix(), 
                                   a21.LockedLocalMatrix() );            
            AllReduce( &myAlpha, &alpha, 1, MPI_SUM, grid.MCComm() );
            BLAS::Axpy( alpha, a21, w21 );
        }
        //--------------------------------------------------------------------//
        a21_MC_Star.FreeConstraints();
        a21_MR_Star.FreeConstraints();
        z10_Star_MR.FreeConstraints();
        z21.FreeConstraints();
        z21_MC_Star.FreeConstraints();
        z21_MR_Star.FreeConstraints();
        z21_MR_MC.FreeConstraints();

        SlidePartitionDown( eT,  e0,
                                 epsilon1,
                           /**/ /********/
                            eB,  e2       );

        SlidePartitionDown( tT,  t0,
                                 tau1,
                           /**/ /****/
                            tB,  t2   );

        SlidePartitionDownDiagonal( ATL, /**/ ATR,  A00, a01,     /**/ A02,
                                         /**/       a10, alpha11, /**/ a12,
                                   /*************/ /**********************/
                                    ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionDownDiagonal( WTL, /**/ WTR,  W00, w01,     /**/ W02,
                                         /**/       w10, omega11, /**/ w12,
                                   /*************/ /**********************/
                                    WBL, /**/ WBR,  W20, w21,     /**/ W22 );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template void Elemental::LAPACK::Internal::PanelTridiagL
( DistMatrix<float,MC,MR  >& A,
  DistMatrix<float,MC,MR  >& W,
  DistMatrix<float,MD,Star>& e,
  DistMatrix<float,MD,Star>& t );

template void Elemental::LAPACK::Internal::PanelTridiagL
( DistMatrix<double,MC,MR  >& A,
  DistMatrix<double,MC,MR  >& W,
  DistMatrix<double,MD,Star>& e,
  DistMatrix<double,MD,Star>& t );

