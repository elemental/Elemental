/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef LAPACK_HERMITIANTRIDIAG_PANELLSQUARE_HPP
#define LAPACK_HERMITIANTRIDIAG_PANELLSQUARE_HPP

#include "elemental/blas-like/level1/Zero.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/lapack-like/Reflector/Col.hpp"

namespace elem {
namespace hermitian_tridiag {

template<typename R> 
void PanelLSquare
( DistMatrix<R>& A,
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR, 
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR )
{
    const int panelSize = W.Width();
    const int bottomSize = W.Height()-panelSize;
    
    const Grid& g = A.Grid();

#ifndef RELEASE
    PushCallStack("hermitian_tridiag::PanelLSquare");
    if( A.Grid() != W.Grid() )
        throw std::logic_error
        ("A and W must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != W.Height() )
        throw std::logic_error("A and W must be the same height");
    if( W.Height() < panelSize )
        throw std::logic_error("W must be a column panel");
    if( W.ColAlignment() != A.ColAlignment() || 
        W.RowAlignment() != A.RowAlignment() )
        throw std::logic_error("W and A must be aligned");
#endif
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Find the process holding our transposed data
    int transposeRank;
    const int r = g.Height();
    {
        const int colAlignment = A.ColAlignment();
        const int rowAlignment = A.RowAlignment();
        const int colShift = A.ColShift();
        const int rowShift = A.RowShift();

        const int transposeRow = (colAlignment+rowShift) % r;
        const int transposeCol = (rowAlignment+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    // Create a distributed matrix for storing the subdiagonal
    DistMatrix<R,MD,STAR> e(g);
    e.AlignWithDiagonal( A, -1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<R> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g),  alpha21T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),            a21B(g),
                         A20(g), a21(g),     A22(g),  A20B(g);
    DistMatrix<R> 
        WTL(g), WTR(g),  W00(g), w01(g),     W02(g),  WCol(g),
        WBL(g), WBR(g),  w10(g), omega11(g), w12(g),
                         W20(g), w21(g),     W22(g),  W20B(g), w21Last(g);
    DistMatrix<R,MD,STAR> eT(g),  e0(g),
                          eB(g),  epsilon1(g),
                                  e2(g);

    // Temporary distributions
    std::vector<R> w21LastBuffer(A.Height()/r+1);
    DistMatrix<R,MC,STAR> a21_MC_STAR(g), a21B_MC_STAR(g), a21Last_MC_STAR(g);
    DistMatrix<R,MR,STAR> a21_MR_STAR(g), a21Last_MR_STAR(g);
    DistMatrix<R,MC,STAR> p21_MC_STAR(g), p21B_MC_STAR(g);
    DistMatrix<R,MR,STAR> p21_MR_STAR(g);
    DistMatrix<R,MR,STAR> q21_MR_STAR(g);
    DistMatrix<R,MR,STAR> x01_MR_STAR(g);
    DistMatrix<R,MR,STAR> y01_MR_STAR(g);
    DistMatrix<R,MC,STAR> w21Last_MC_STAR(g);
    DistMatrix<R,MR,STAR> w21Last_MR_STAR(g);

    // Push to the blocksize of 1, then pop at the end of the routine
    PushBlocksizeStack( 1 );

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDownLeftDiagonal
    ( W, WTL, WTR,
         WBL, WBR, 0 );
    PartitionDown
    ( e, eT,
         eB, 0 );
    bool firstIteration = true;
    R tau = 0;
    R w21LastFirstEntry = 0;
    while( WTL.Width() < panelSize )
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

        View2x1
        ( ACol, alpha11,
                a21 );

        View2x1
        ( WCol, omega11,
                w21 );

        // View the portions of A20 and W20 outside of this panel's square
        View( A20B, A, panelSize, 0, bottomSize, A20.Width() );
        View( W20B, W, panelSize, 0, bottomSize, W20.Width() );

        if( !firstIteration )
        {
            View
            ( a21Last_MC_STAR,
              APan_MC_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            View
            ( a21Last_MR_STAR,
              APan_MR_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            View( w21Last, W, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
        }
            
        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        a21_MC_STAR.AlignWith( A22 );
        a21_MR_STAR.AlignWith( A22 );
        p21_MC_STAR.AlignWith( A22 );
        p21_MR_STAR.AlignWith( A22 );
        q21_MR_STAR.AlignWith( A22 );
        x01_MR_STAR.AlignWith( W20B );
        y01_MR_STAR.AlignWith( W20B );

        a21_MC_STAR.ResizeTo( a21.Height(), 1 );
        a21_MR_STAR.ResizeTo( a21.Height(), 1 );
        p21_MC_STAR.ResizeTo( a21.Height(), 1 );
        p21_MR_STAR.ResizeTo( a21.Height(), 1 );
        q21_MR_STAR.ResizeTo( a21.Height(), 1 );
        x01_MR_STAR.ResizeTo( W20B.Width(), 1 );
        y01_MR_STAR.ResizeTo( W20B.Width(), 1 );

        // View the portions of a21[MC,* ] and p21[MC,* ] below the current
        // panel's square
        View
        ( a21B_MC_STAR,
          a21_MC_STAR, a21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        View
        ( p21B_MC_STAR,
          p21_MC_STAR, p21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const int AColLocalHeight = ACol.LocalHeight();
                R* AColBuffer = ACol.Buffer();
                const R* a21Last_MC_STAR_Buffer = a21Last_MC_STAR.Buffer();
                for( int i=0; i<AColLocalHeight; ++i )
                    AColBuffer[i] -=
                        w21LastBuffer[i] + 
                        a21Last_MC_STAR_Buffer[i]*w21LastFirstEntry;
            }
            // Compute the Householder reflector
            tau = reflector::Col( alpha21T, a21B );
        }

        // Store the subdiagonal value and turn a21 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry
        alpha21T.GetDiagonal( epsilon1 );
        alpha21T.Set( 0, 0, R(1) );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a21 within its process row. 
        // Otherwise, also add w21 into the broadcast.
        if( firstIteration )
        {
            const int a21LocalHeight = a21.LocalHeight();
            std::vector<R> rowBroadcastBuffer(a21LocalHeight+1);
            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a21 and tau
                MemCopy( &rowBroadcastBuffer[0], a21.Buffer(), a21LocalHeight );
                rowBroadcastBuffer[a21LocalHeight] = tau;
            }
            // Broadcast a21 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+1, a21.RowAlignment(), g.RowComm() );
            // Store a21[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy
            ( a21_MC_STAR.Buffer(), &rowBroadcastBuffer[0], a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,0), 
              &rowBroadcastBuffer[0],
              APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset );
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight];
            
            // Take advantage of the square process grid in order to form 
            // a21[MR,* ] from a21[MC,* ]
            if( onDiagonal )
            {
                MemCopy
                ( a21_MR_STAR.Buffer(), a21_MC_STAR.Buffer(), a21LocalHeight );
            }
            else
            {
                // Pairwise exchange
                const int sendSize = A22.LocalHeight();
                const int recvSize = A22.LocalWidth();
                mpi::SendRecv
                ( a21_MC_STAR.Buffer(), sendSize, transposeRank, 0,
                  a21_MR_STAR.Buffer(), recvSize, transposeRank, 0, 
                  g.VCComm() );
            }
            // Store a21[MR,* ]
            const int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            MemCopy
            ( APan_MR_STAR.Buffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.Buffer(),
              APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset );
        }
        else
        {
            const int a21LocalHeight = a21.LocalHeight();
            const int w21LastLocalHeight = ACol.LocalHeight();
            std::vector<R> 
                rowBroadcastBuffer(a21LocalHeight+w21LastLocalHeight+1);
            if( thisIsMyCol ) 
            {
                // Pack the broadcast buffer with a21, w21Last, and tau
                MemCopy( &rowBroadcastBuffer[0], a21.Buffer(), a21LocalHeight );
                MemCopy
                ( &rowBroadcastBuffer[a21LocalHeight], 
                  &w21LastBuffer[0], w21LastLocalHeight );
                rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight] = tau;
            }
            // Broadcast a21, w21Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+w21LastLocalHeight+1, 
              a21.RowAlignment(), g.RowComm() );
            // Store a21[MC,* ] into its DistMatrix class 
            MemCopy
            ( a21_MC_STAR.Buffer(), &rowBroadcastBuffer[0], a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,A00.Width()), 
              &rowBroadcastBuffer[0],
              APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset );
            // Store w21Last[MC,* ] into its DistMatrix class
            w21Last_MC_STAR.AlignWith( alpha11 );
            w21Last_MC_STAR.ResizeTo( a21.Height()+1, 1 );
            MemCopy
            ( w21Last_MC_STAR.Buffer(), 
              &rowBroadcastBuffer[a21LocalHeight], w21LastLocalHeight );
            // Store the bottom part of w21Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w21.
            const int W_MC_STAR_Offset = 
                W_MC_STAR.LocalHeight()-w21LastLocalHeight;
            MemCopy
            ( W_MC_STAR.Buffer(W_MC_STAR_Offset,A00.Width()-1),
              &rowBroadcastBuffer[a21LocalHeight],
              W_MC_STAR.LocalHeight()-W_MC_STAR_Offset );
            if( g.Col() == w21Last.RowAlignment() )
            {
                MemCopy
                ( w21Last.Buffer(),
                  &rowBroadcastBuffer[a21LocalHeight], w21LastLocalHeight );
            }
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight];

            // Take advantage of the square process grid in order to quickly
            // form a21[MR,* ] and w21Last[MR,* ] from their [MC,* ] 
            // counterparts
            w21Last_MR_STAR.AlignWith( ABR );
            w21Last_MR_STAR.ResizeTo( w21Last.Height(), 1 );
            if( onDiagonal )
            {
                MemCopy
                ( a21_MR_STAR.Buffer(), a21_MC_STAR.Buffer(), a21LocalHeight );
                MemCopy
                ( w21Last_MR_STAR.Buffer(),
                  w21Last_MC_STAR.Buffer(), w21LastLocalHeight );
            }
            else
            {
                const int sendSize = A22.LocalHeight()+ABR.LocalHeight();
                const int recvSize = A22.LocalWidth()+ABR.LocalWidth();
                std::vector<R> sendBuffer(sendSize), recvBuffer(recvSize);

                // Pack the send buffer
                MemCopy
                ( &sendBuffer[0], a21_MC_STAR.Buffer(), A22.LocalHeight() );
                MemCopy
                ( &sendBuffer[A22.LocalHeight()],
                  w21Last_MC_STAR.Buffer(), ABR.LocalHeight() );

                // Pairwise exchange
                mpi::SendRecv
                ( &sendBuffer[0], sendSize, transposeRank, 0,
                  &recvBuffer[0], recvSize, transposeRank, 0,
                  g.VCComm() );

                // Unpack the recv buffer
                MemCopy
                ( a21_MR_STAR.Buffer(), &recvBuffer[0], A22.LocalWidth() );
                MemCopy
                ( w21Last_MR_STAR.Buffer(),
                  &recvBuffer[A22.LocalWidth()], ABR.LocalWidth() );
            }

            // Store w21Last[MR,* ]
            const int W_MR_STAR_Offset = 
                W_MR_STAR.LocalHeight()-w21Last_MR_STAR.LocalHeight();
            MemCopy
            ( W_MR_STAR.Buffer(W_MR_STAR_Offset,A00.Width()-1),
              w21Last_MR_STAR.Buffer(),
              W_MR_STAR.LocalHeight()-W_MR_STAR_Offset );
            // Store a21[MR,* ]
            const int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            MemCopy
            ( APan_MR_STAR.Buffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.Buffer(),
              APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset );

            // Update the portion of A22 that is in our current panel with 
            // w21Last and a21Last using two gers. We do not need their top 
            // entries. We trash the upper triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            DistMatrix<R,MC,STAR> a21Last_MC_STAR_Bottom(g),
                                  w21Last_MC_STAR_Bottom(g);
            DistMatrix<R,MR,STAR> a21Last_MR_STAR_Bottom(g),
                                  w21Last_MR_STAR_Bottom(g);
            View
            ( a21Last_MC_STAR_Bottom,
              a21Last_MC_STAR, 1, 0, a21Last_MC_STAR.Height()-1, 1 );
            View
            ( w21Last_MC_STAR_Bottom,
              w21Last_MC_STAR, 1, 0, w21Last_MC_STAR.Height()-1, 1 );
            View
            ( a21Last_MR_STAR_Bottom,
              a21Last_MR_STAR, 1, 0, a21Last_MR_STAR.Height()-1, 1 );
            View
            ( w21Last_MR_STAR_Bottom,
              w21Last_MR_STAR, 1, 0, w21Last_MR_STAR.Height()-1, 1 );
            const R* a21_MC_STAR_Buffer = a21Last_MC_STAR_Bottom.Buffer();
            const R* w21_MC_STAR_Buffer = w21Last_MC_STAR_Bottom.Buffer();
            const R* a21_MR_STAR_Buffer = a21Last_MR_STAR_Bottom.Buffer();
            const R* w21_MR_STAR_Buffer = w21Last_MR_STAR_Bottom.Buffer();
            R* A22Buffer = A22.Buffer();
            const int localHeight = W22.LocalHeight();
            const int localWidth = W22.LocalWidth();
            const int lDim = A22.LDim();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    A22Buffer[iLocal+jLocal*lDim] -=
                        w21_MC_STAR_Buffer[iLocal]*a21_MR_STAR_Buffer[jLocal] +
                        a21_MC_STAR_Buffer[iLocal]*w21_MR_STAR_Buffer[jLocal];

            // We are through with the last iteration's w21
            w21Last_MC_STAR.FreeAlignments();
            w21Last_MR_STAR.FreeAlignments();
        }

        // Form the local portions of (A22 a21) into p21[MC,* ] and q21[MR,* ]:
        //   p21[MC,* ] := tril(A22)[MC,MR] a21[MR,* ]
        //   q21[MR,* ] := tril(A22,-1)'[MR,MC] a21[MC,* ]
        if( A22.ColShift() > A22.RowShift() )
        {
            // We are below the diagonal, so we can multiply without an 
            // offset for tril(A22)[MC,MR] and tril(A22,-1)'[MR,MC]
            if( A22.LocalHeight() != 0 )
            {
                MemCopy
                ( p21_MC_STAR.Buffer(),
                  a21_MR_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(), 
                  p21_MC_STAR.Buffer(), 1 );
            }
            if( A22.LocalWidth() != 0 )
            {
                // Our local portion of q21[MR,* ] might be one entry longer 
                // than A22.LocalHeight(), so go ahead and set the last entry 
                // to 0.
                R* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer,
                  a21_MC_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'T', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }
        else if( A22.ColShift() < A22.RowShift() )
        {
            // We are above the diagonal, so we need to use an offset of +1
            // for both tril(A22)[MC,MR] and tril(A22,-1)'[MR,MC]
            const R* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
            const R* a21_MR_STAR_Buffer = a21_MR_STAR.Buffer();
            const R* A22Buffer = A22.Buffer();
            if( A22.LocalHeight() != 0 )
            {
                // The first entry of p21[MC,* ] will always be zero due to 
                // the forced offset.
                R* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
                p21_MC_STAR_Buffer[0] = 0;
                MemCopy
                ( &p21_MC_STAR_Buffer[1], a21_MR_STAR_Buffer, 
                  A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  &p21_MC_STAR_Buffer[1], 1 );
             }
             if( A22.LocalWidth() != 0 )
             {
                // The last entry of q21[MR,* ] must be zero due to the forced
                // offset.
                R* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer,
                  &a21_MC_STAR_Buffer[1], A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'T', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }
        else
        {
            // We are on the diagonal, so we only need an offset of +1 for
            // tril(A22,-1)'[MR,MC]
            if( A22.LocalHeight() != 0 )
            {
                MemCopy
                ( p21_MC_STAR.Buffer(),
                  a21_MR_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(), 
                  p21_MC_STAR.Buffer(), 1 );
 
                // The last entry of q21[MR,* ] must be zero due to the forced
                // offset.
                const R* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
                const R* A22Buffer = A22.Buffer();
                R* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer,
                  &a21_MC_STAR_Buffer[1], A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'T', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }

        LocalGemv( TRANSPOSE, R(1), W20B, a21B_MC_STAR, R(0), x01_MR_STAR );
        LocalGemv( TRANSPOSE, R(1), A20B, a21B_MC_STAR, R(0), y01_MR_STAR );

        // Combine the AllReduce column summations of x01[MR,* ] and y01[MR,* ]
        {
            const int x01LocalHeight = x01_MR_STAR.LocalHeight();
            std::vector<R> colSumSendBuffer(2*x01LocalHeight),
                           colSumRecvBuffer(2*x01LocalHeight);
            MemCopy
            ( &colSumSendBuffer[0], x01_MR_STAR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x01LocalHeight],
              y01_MR_STAR.Buffer(), x01LocalHeight );
            mpi::AllReduce
            ( &colSumSendBuffer[0], 
              &colSumRecvBuffer[0],
              2*x01LocalHeight, mpi::SUM, g.ColComm() );
            MemCopy
            ( x01_MR_STAR.Buffer(), &colSumRecvBuffer[0], x01LocalHeight );
            MemCopy
            ( y01_MR_STAR.Buffer(), 
              &colSumRecvBuffer[x01LocalHeight], x01LocalHeight );
        }

        LocalGemv( NORMAL, R(-1), A20B, x01_MR_STAR, R(1), p21B_MC_STAR );
        LocalGemv( NORMAL, R(-1), W20B, y01_MR_STAR, R(1), p21B_MC_STAR );

        // Fast transpose the unsummed q21[MR,* ] -> q21[MC,* ], so that
        // it needs to be summed over process rows instead of process 
        // columns. We immediately add it onto p21[MC,* ], which also needs
        // to be summed within process rows.
        if( onDiagonal )
        {
            const int a21LocalHeight = a21.LocalHeight();
            R* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
            const R* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
            for( int i=0; i<a21LocalHeight; ++i )
                p21_MC_STAR_Buffer[i] += q21_MR_STAR_Buffer[i];
        }
        else
        {
            // Pairwise exchange with the transpose process
            const int sendSize = A22.LocalWidth();
            const int recvSize = A22.LocalHeight();
            std::vector<R> recvBuffer(recvSize);
            mpi::SendRecv
            ( q21_MR_STAR.Buffer(), sendSize, transposeRank, 0,
              &recvBuffer[0],       recvSize, transposeRank, 0, g.VCComm() );

            // Unpack the recv buffer directly onto p21[MC,* ]
            R* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
            for( int i=0; i<recvSize; ++i )
                p21_MC_STAR_Buffer[i] += recvBuffer[i];
        }

        if( W22.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p21[MC,* ] to the next process column.
            const int a21LocalHeight = a21.LocalHeight();

            const int nextProcessRow = (alpha11.ColAlignment()+1) % r;
            const int nextProcessCol = (alpha11.RowAlignment()+1) % r;

            std::vector<R> reduceToOneRecvBuffer(a21LocalHeight);
            mpi::Reduce
            ( p21_MC_STAR.Buffer(),
              &reduceToOneRecvBuffer[0],
              a21LocalHeight, mpi::SUM, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Finish computing w21. During its computation, ensure that 
                // every process has a copy of the first element of the w21.
                // We know a priori that the first element of a21 is one.
                R myDotProduct = blas::Dot
                    ( a21LocalHeight, &reduceToOneRecvBuffer[0], 1, 
                                      a21_MC_STAR.Buffer(), 1 );
                R sendBuffer[2], recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.Row()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[0] : 0 );
                mpi::AllReduce
                ( sendBuffer, recvBuffer, 2, mpi::SUM, g.ColComm() );
                R dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w21LastBuffer
                // - w21LastFirstEntry
                R scale = dotProduct*tau / R(2);
                const R* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
                for( int i=0; i<a21LocalHeight; ++i )
                    w21LastBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a21_MC_STAR_Buffer[i] );
                w21LastFirstEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, so we just need to form w21[MC,* ]
            // and w21[MR,* ]. 
            const int a21LocalHeight = a21.LocalHeight();

            // AllReduce sum p21[MC,* ] over process rows
            std::vector<R> allReduceRecvBuffer(a21LocalHeight);
            mpi::AllReduce
            ( p21_MC_STAR.Buffer(), &allReduceRecvBuffer[0],
              a21LocalHeight, mpi::SUM, g.RowComm() );

            // Finish computing w21.
            R myDotProduct = blas::Dot
                ( a21LocalHeight, &allReduceRecvBuffer[0], 1, 
                                  a21_MC_STAR.Buffer(),    1 );
            R dotProduct;
            mpi::AllReduce
            ( &myDotProduct, &dotProduct, 1, mpi::SUM, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<R,MC,STAR> w21_MC_STAR(g);
            DistMatrix<R,MR,STAR> w21_MR_STAR(g);
            View
            ( w21_MC_STAR,
              W_MC_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );
            View
            ( w21_MR_STAR, 
              W_MR_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );

            // Store w21[MC,* ]
            R scale = dotProduct*tau / R(2);
            R* w21_MC_STAR_Buffer = w21_MC_STAR.Buffer();
            const R* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
            for( int i=0; i<a21LocalHeight; ++i )
                w21_MC_STAR_Buffer[i] = 
                    tau*( allReduceRecvBuffer[i]-scale*a21_MC_STAR_Buffer[i] );

            // Fast transpose w21[MC,* ] -> w21[MR,* ]
            if( onDiagonal )
            {
                MemCopy
                ( w21_MR_STAR.Buffer(), w21_MC_STAR.Buffer(), a21LocalHeight );
            }
            else
            {
                // Pairwise exchange with the transpose process
                const int sendSize = A22.LocalHeight();
                const int recvSize = A22.LocalWidth();
                mpi::SendRecv
                ( w21_MC_STAR.Buffer(), sendSize, transposeRank, 0,
                  w21_MR_STAR.Buffer(), recvSize, transposeRank, 0,
                  g.VCComm() );
            }
        }
        //--------------------------------------------------------------------//
        a21_MC_STAR.FreeAlignments();
        a21_MR_STAR.FreeAlignments();
        p21_MC_STAR.FreeAlignments();
        p21_MR_STAR.FreeAlignments();
        q21_MR_STAR.FreeAlignments();
        x01_MR_STAR.FreeAlignments();
        y01_MR_STAR.FreeAlignments();

        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        SlidePartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        firstIteration = false;
    }
    PopBlocksizeStack();

    // View the portion of A that e is the subdiagonal of, then place e into it
    DistMatrix<R> expandedATL(g);
    View( expandedATL, A, 0, 0, panelSize+1, panelSize+1 );
    expandedATL.SetDiagonal( e, -1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void PanelLSquare
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR, 
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR )
{
    typedef Complex<R> C;

    const int panelSize = W.Width();
    const int bottomSize = W.Height()-panelSize;

    const Grid& g = A.Grid();

#ifndef RELEASE
    PushCallStack("hermitian_tridiag::PanelLSquare");
    if( A.Grid() != W.Grid() || W.Grid() != t.Grid() )
        throw std::logic_error
        ("A, W, and t must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != W.Height() )
        throw std::logic_error("A and W must be the same height");
    if( W.Height() < panelSize )
        throw std::logic_error("W must be a column panel");
    if( W.ColAlignment() != A.ColAlignment() || 
        W.RowAlignment() != A.RowAlignment() )
        throw std::logic_error("W and A must be aligned");
    if( t.Height() != W.Width() || t.Width() != 1 )
        throw std::logic_error
        ("t must be a column vector of the same length as W's width");
    if( !t.AlignedWithDiagonal(A,-1) )
        throw std::logic_error("t is not aligned with A's subdiagonal");
#endif
    if( !g.InGrid() )
    {
#ifndef RELEASE
        PopCallStack();
#endif
        return;
    }

    // Find the process holding our transposed data
    const int r = g.Height();
    int transposeRank;
    {
        const int colAlignment = A.ColAlignment();
        const int rowAlignment = A.RowAlignment();
        const int colShift = A.ColShift();
        const int rowShift = A.RowShift();

        const int transposeRow = (colAlignment+rowShift) % r;
        const int transposeCol = (rowAlignment+colShift) % r;
        transposeRank = transposeRow + r*transposeCol;
    }
    const bool onDiagonal = ( transposeRank == g.VCRank() );

    // Create a distributed matrix for storing the subdiagonal
    DistMatrix<R,MD,STAR> e(g);
    e.AlignWithDiagonal( A.DistData(), -1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<C> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g),  alpha21T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),            a21B(g),
                         A20(g), a21(g),     A22(g),  A20B(g);
    DistMatrix<C> 
        WTL(g), WTR(g),  W00(g), w01(g),     W02(g),  WCol(g),
        WBL(g), WBR(g),  w10(g), omega11(g), w12(g),
                         W20(g), w21(g),     W22(g),  W20B(g), w21Last(g);
    DistMatrix<R,MD,STAR> eT(g),  e0(g),
                          eB(g),  epsilon1(g),
                                  e2(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    std::vector<C> w21LastBuffer(A.Height()/r+1);
    DistMatrix<C,MC,STAR> a21_MC_STAR(g), a21B_MC_STAR(g), a21Last_MC_STAR(g);
    DistMatrix<C,MR,STAR> a21_MR_STAR(g), a21Last_MR_STAR(g);
    DistMatrix<C,MC,STAR> p21_MC_STAR(g), p21B_MC_STAR(g);
    DistMatrix<C,MR,STAR> p21_MR_STAR(g);
    DistMatrix<C,MR,STAR> q21_MR_STAR(g);
    DistMatrix<C,MR,STAR> x01_MR_STAR(g);
    DistMatrix<C,MR,STAR> y01_MR_STAR(g);
    DistMatrix<C,MC,STAR> w21Last_MC_STAR(g);
    DistMatrix<C,MR,STAR> w21Last_MR_STAR(g);

    // Push to the blocksize of 1, then pop at the end of the routine
    PushBlocksizeStack( 1 );

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionDownLeftDiagonal
    ( W, WTL, WTR,
         WBL, WBR, 0 );
    PartitionDown
    ( e, eT,
         eB, 0 );
    PartitionDown
    ( t, tT,
         tB, 0 );
    bool firstIteration = true;
    C tau = 0;
    C w21LastFirstEntry = 0;
    while( WTL.Width() < panelSize )
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

        View2x1
        ( ACol, alpha11,
                a21 );

        View2x1
        ( WCol, omega11,
                w21 );

        // View the portions of A20 and W20 outside of this panel's square
        View( A20B, A, panelSize, 0, bottomSize, A20.Width() );
        View( W20B, W, panelSize, 0, bottomSize, W20.Width() );

        if( !firstIteration )
        {
            View
            ( a21Last_MC_STAR,
              APan_MC_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            View
            ( a21Last_MR_STAR,
              APan_MR_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            View( w21Last, W, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
        }
            
        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        a21_MC_STAR.AlignWith( A22 );
        a21_MR_STAR.AlignWith( A22 );
        p21_MC_STAR.AlignWith( A22 );
        p21_MR_STAR.AlignWith( A22 );
        q21_MR_STAR.AlignWith( A22 );
        x01_MR_STAR.AlignWith( W20B );
        y01_MR_STAR.AlignWith( W20B );

        a21_MC_STAR.ResizeTo( a21.Height(), 1 );
        a21_MR_STAR.ResizeTo( a21.Height(), 1 );
        p21_MC_STAR.ResizeTo( a21.Height(), 1 );
        p21_MR_STAR.ResizeTo( a21.Height(), 1 );
        q21_MR_STAR.ResizeTo( a21.Height(), 1 );
        x01_MR_STAR.ResizeTo( W20B.Width(), 1 );
        y01_MR_STAR.ResizeTo( W20B.Width(), 1 );

        // View the portions of a21[MC,* ] and p21[MC,* ] below the current
        // panel's square
        View
        ( a21B_MC_STAR,
          a21_MC_STAR, a21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        View
        ( p21B_MC_STAR,
          p21_MC_STAR, p21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const int AColLocalHeight = ACol.LocalHeight();
                C* AColBuffer = ACol.Buffer();
                const C* a21Last_MC_STAR_Buffer = a21Last_MC_STAR.Buffer();
                for( int i=0; i<AColLocalHeight; ++i )
                    AColBuffer[i] -=
                        w21LastBuffer[i] + 
                        a21Last_MC_STAR_Buffer[i]*Conj(w21LastFirstEntry);
            }
            // Compute the Householder reflector
            tau = reflector::Col( alpha21T, a21B );
            if( g.Row() == alpha21T.ColAlignment() )
                tau1.SetLocal(0,0,tau);
        }

        // Store the subdiagonal value and turn a21 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        alpha21T.GetRealPartOfDiagonal( epsilon1 );
        alpha21T.Set( 0, 0, C(1) );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a21 within its process row. 
        // Otherwise, also add w21 into the broadcast.
        if( firstIteration )
        {
            const int a21LocalHeight = a21.LocalHeight();
            std::vector<C> rowBroadcastBuffer(a21LocalHeight+1);
            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a21 and tau
                MemCopy
                ( &rowBroadcastBuffer[0], a21.Buffer(), a21LocalHeight );
                rowBroadcastBuffer[a21LocalHeight] = tau;
            }
            // Broadcast a21 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+1, a21.RowAlignment(), g.RowComm() );
            // Store a21[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy
            ( a21_MC_STAR.Buffer(), &rowBroadcastBuffer[0], a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,0), 
              &rowBroadcastBuffer[0],
              APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset );
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight];
            
            // Take advantage of the square process grid in order to form 
            // a21[MR,* ] from a21[MC,* ]
            if( onDiagonal )
            {
                MemCopy
                ( a21_MR_STAR.Buffer(), a21_MC_STAR.Buffer(), a21LocalHeight );
            }
            else
            {
                // Pairwise exchange
                const int sendSize = A22.LocalHeight();
                const int recvSize = A22.LocalWidth();
                mpi::SendRecv
                ( a21_MC_STAR.Buffer(), sendSize, transposeRank, 0,
                  a21_MR_STAR.Buffer(), recvSize, transposeRank, 0,
                  g.VCComm() );
            }
            // Store a21[MR,* ]
            const int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            MemCopy
            ( APan_MR_STAR.Buffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.Buffer(),
              APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset );
        }
        else
        {
            const int a21LocalHeight = a21.LocalHeight();
            const int w21LastLocalHeight = ACol.LocalHeight();
            std::vector<C> 
                rowBroadcastBuffer(a21LocalHeight+w21LastLocalHeight+1);
            if( thisIsMyCol ) 
            {
                // Pack the broadcast buffer with a21, w21Last, and tau
                MemCopy
                ( &rowBroadcastBuffer[0], a21.Buffer(), a21LocalHeight );
                MemCopy
                ( &rowBroadcastBuffer[a21LocalHeight], 
                  &w21LastBuffer[0], w21LastLocalHeight );
                rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight] = tau;
            }
            // Broadcast a21, w21Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+w21LastLocalHeight+1, 
              a21.RowAlignment(), g.RowComm() );
            // Store a21[MC,* ] into its DistMatrix class 
            MemCopy
            ( a21_MC_STAR.Buffer(), &rowBroadcastBuffer[0], a21LocalHeight );
            // Store a21[MC,* ] into APan[MC,* ]
            const int APan_MC_STAR_Offset = 
                APan_MC_STAR.LocalHeight()-a21LocalHeight;
            MemCopy
            ( APan_MC_STAR.Buffer(APan_MC_STAR_Offset,A00.Width()), 
              &rowBroadcastBuffer[0],
              APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset );
            // Store w21Last[MC,* ] into its DistMatrix class
            w21Last_MC_STAR.AlignWith( alpha11 );
            w21Last_MC_STAR.ResizeTo( a21.Height()+1, 1 );
            MemCopy
            ( w21Last_MC_STAR.Buffer(), 
              &rowBroadcastBuffer[a21LocalHeight], w21LastLocalHeight );
            // Store the bottom part of w21Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w21.
            const int W_MC_STAR_Offset = 
                W_MC_STAR.LocalHeight()-w21LastLocalHeight;
            MemCopy
            ( W_MC_STAR.Buffer(W_MC_STAR_Offset,A00.Width()-1),
              &rowBroadcastBuffer[a21LocalHeight],
              W_MC_STAR.LocalHeight()-W_MC_STAR_Offset );
            if( g.Col() == w21Last.RowAlignment() )
            {
                MemCopy
                ( w21Last.Buffer(),
                  &rowBroadcastBuffer[a21LocalHeight], w21LastLocalHeight );
            }
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight];

            // Take advantage of the square process grid in order to quickly
            // form a21[MR,* ] and w21Last[MR,* ] from their [MC,* ] 
            // counterparts
            w21Last_MR_STAR.AlignWith( ABR );
            w21Last_MR_STAR.ResizeTo( w21Last.Height(), 1 );
            if( onDiagonal )
            {
                MemCopy
                ( a21_MR_STAR.Buffer(), a21_MC_STAR.Buffer(), a21LocalHeight );
                MemCopy
                ( w21Last_MR_STAR.Buffer(),
                  w21Last_MC_STAR.Buffer(), w21LastLocalHeight );
            }
            else
            {
                const int sendSize = A22.LocalHeight()+ABR.LocalHeight();
                const int recvSize = A22.LocalWidth()+ABR.LocalWidth();
                std::vector<C> sendBuffer(sendSize), recvBuffer(recvSize);

                // Pack the send buffer
                MemCopy
                ( &sendBuffer[0], a21_MC_STAR.Buffer(), A22.LocalHeight() );
                MemCopy
                ( &sendBuffer[A22.LocalHeight()],
                  w21Last_MC_STAR.Buffer(), ABR.LocalHeight() );

                // Pairwise exchange
                mpi::SendRecv
                ( &sendBuffer[0], sendSize, transposeRank, 0,
                  &recvBuffer[0], recvSize, transposeRank, 0,
                  g.VCComm() );

                // Unpack the recv buffer
                MemCopy
                ( a21_MR_STAR.Buffer(), &recvBuffer[0], A22.LocalWidth() );
                MemCopy
                ( w21Last_MR_STAR.Buffer(),
                  &recvBuffer[A22.LocalWidth()], ABR.LocalWidth() );
            }

            // Store w21Last[MR,* ]
            const int W_MR_STAR_Offset = 
                W_MR_STAR.LocalHeight()-w21Last_MR_STAR.LocalHeight();
            MemCopy
            ( W_MR_STAR.Buffer(W_MR_STAR_Offset,A00.Width()-1),
              w21Last_MR_STAR.Buffer(),
              W_MR_STAR.LocalHeight()-W_MR_STAR_Offset );
            // Store a21[MR,* ]
            const int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            MemCopy
            ( APan_MR_STAR.Buffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.Buffer(),
              APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset );

            // Update the portion of A22 that is in our current panel with 
            // w21Last and a21Last using two gers. We do not need their top 
            // entries. We trash the upper triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            DistMatrix<C,MC,STAR> a21Last_MC_STAR_Bottom(g),
                                  w21Last_MC_STAR_Bottom(g);
            DistMatrix<C,MR,STAR> a21Last_MR_STAR_Bottom(g),
                                  w21Last_MR_STAR_Bottom(g);
            View
            ( a21Last_MC_STAR_Bottom,
              a21Last_MC_STAR, 1, 0, a21Last_MC_STAR.Height()-1, 1 );
            View
            ( w21Last_MC_STAR_Bottom,
              w21Last_MC_STAR, 1, 0, w21Last_MC_STAR.Height()-1, 1 );
            View
            ( a21Last_MR_STAR_Bottom,
              a21Last_MR_STAR, 1, 0, a21Last_MR_STAR.Height()-1, 1 );
            View
            ( w21Last_MR_STAR_Bottom,
              w21Last_MR_STAR, 1, 0, w21Last_MR_STAR.Height()-1, 1 );
            const C* a21_MC_STAR_Buffer = a21Last_MC_STAR_Bottom.Buffer();
            const C* w21_MC_STAR_Buffer = w21Last_MC_STAR_Bottom.Buffer();
            const C* a21_MR_STAR_Buffer = a21Last_MR_STAR_Bottom.Buffer();
            const C* w21_MR_STAR_Buffer = w21Last_MR_STAR_Bottom.Buffer();
            C* A22Buffer = A22.Buffer();
            const int localHeight = W22.LocalHeight();
            const int localWidth = W22.LocalWidth();
            const int lDim = A22.LDim();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    A22Buffer[iLocal+jLocal*lDim] -=
                        w21_MC_STAR_Buffer[iLocal]*
                        Conj(a21_MR_STAR_Buffer[jLocal]) +
                        a21_MC_STAR_Buffer[iLocal]*
                        Conj(w21_MR_STAR_Buffer[jLocal]);

            // We are through with the last iteration's w21
            w21Last_MC_STAR.FreeAlignments();
            w21Last_MR_STAR.FreeAlignments();
        }

        // Form the local portions of (A22 a21) into p21[MC,* ] and q21[MR,* ]:
        //   p21[MC,* ] := tril(A22)[MC,MR] a21[MR,* ]
        //   q21[MR,* ] := tril(A22,-1)'[MR,MC] a21[MC,* ]
        if( A22.ColShift() > A22.RowShift() )
        {
            // We are below the diagonal, so we can multiply without an 
            // offset for tril(A22)[MC,MR] and tril(A22,-1)'[MR,MC]
            if( A22.LocalHeight() != 0 )
            {
                MemCopy
                ( p21_MC_STAR.Buffer(),
                  a21_MR_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(), 
                  p21_MC_STAR.Buffer(), 1 );
            }
            if( A22.LocalWidth() != 0 )
            {
                // Our local portion of q21[MR,* ] might be one entry longer 
                // than A22.LocalHeight(), so go ahead and set the last entry 
                // to 0.
                C* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer, a21_MC_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'C', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }
        else if( A22.ColShift() < A22.RowShift() )
        {
            // We are above the diagonal, so we need to use an offset of +1
            // for both tril(A22)[MC,MR] and tril(A22,-1)'[MR,MC]
            const C* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
            const C* a21_MR_STAR_Buffer = a21_MR_STAR.Buffer();
            const C* A22Buffer = A22.Buffer();
            if( A22.LocalHeight() != 0 )
            {
                // The first entry of p21[MC,* ] will always be zero due to 
                // the forced offset.
                C* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
                p21_MC_STAR_Buffer[0] = 0;
                MemCopy
                ( &p21_MC_STAR_Buffer[1],
                  a21_MR_STAR_Buffer, A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  &p21_MC_STAR_Buffer[1], 1 );
            }
            if( A22.LocalWidth() != 0 )
            {
                // The last entry of q21[MR,* ] will be zero if the local
                // height and width are equal.
                C* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer,
                  &a21_MC_STAR_Buffer[1], A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'C', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }
        else
        {
            // We are on the diagonal, so we only need an offset of +1 for
            // tril(A22,-1)'[MR,MC]
            if( A22.LocalHeight() != 0 )
            {
                MemCopy
                ( p21_MC_STAR.Buffer(),
                  a21_MR_STAR.Buffer(), A22.LocalHeight() );
                blas::Trmv
                ( 'L', 'N', 'N', A22.LocalHeight(), A22.Buffer(), A22.LDim(), 
                  p21_MC_STAR.Buffer(), 1 );

                // The last entry of q21[MR,* ] will be zero if the local 
                // height and width are equal.
                const C* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
                const C* A22Buffer = A22.Buffer();
                C* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
                q21_MR_STAR_Buffer[A22.LocalWidth()-1] = 0;
                MemCopy
                ( q21_MR_STAR_Buffer,
                  &a21_MC_STAR_Buffer[1], A22.LocalHeight()-1 );
                blas::Trmv
                ( 'L', 'C', 'N', A22.LocalHeight()-1, &A22Buffer[1], A22.LDim(),
                  q21_MR_STAR_Buffer, 1 );
            }
        }

        LocalGemv( ADJOINT, C(1), W20B, a21B_MC_STAR, C(0), x01_MR_STAR );
        LocalGemv( ADJOINT, C(1), A20B, a21B_MC_STAR, C(0), y01_MR_STAR );

        // Combine the AllReduce column summations of x01[MR,* ] and 
        // y01[MR,* ]
        {
            const int x01LocalHeight = x01_MR_STAR.LocalHeight();
            std::vector<C> colSumSendBuffer(2*x01LocalHeight),
                           colSumRecvBuffer(2*x01LocalHeight);
            MemCopy
            ( &colSumSendBuffer[0], 
              x01_MR_STAR.Buffer(), x01LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x01LocalHeight],
              y01_MR_STAR.Buffer(), x01LocalHeight );
            mpi::AllReduce
            ( &colSumSendBuffer[0], 
              &colSumRecvBuffer[0],
              2*x01LocalHeight, mpi::SUM, g.ColComm() );
            MemCopy
            ( x01_MR_STAR.Buffer(), 
              &colSumRecvBuffer[0], x01LocalHeight );
            MemCopy
            ( y01_MR_STAR.Buffer(), 
              &colSumRecvBuffer[x01LocalHeight], x01LocalHeight );
        }

        LocalGemv( NORMAL, C(-1), A20B, x01_MR_STAR, C(1), p21B_MC_STAR );
        LocalGemv( NORMAL, C(-1), W20B, y01_MR_STAR, C(1), p21B_MC_STAR );

        // Fast transpose the unsummed q21[MR,* ] -> q21[MC,* ], so that
        // it needs to be summed over process rows instead of process 
        // columns. We immediately add it onto p21[MC,* ], which also needs
        // to be summed within process rows.
        if( onDiagonal )
        {
            const int a21LocalHeight = a21.LocalHeight();
            C* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
            const C* q21_MR_STAR_Buffer = q21_MR_STAR.Buffer();
            for( int i=0; i<a21LocalHeight; ++i )
                p21_MC_STAR_Buffer[i] += q21_MR_STAR_Buffer[i];
        }
        else
        {
            // Pairwise exchange with the transpose process
            const int sendSize = A22.LocalWidth();
            const int recvSize = A22.LocalHeight();
            std::vector<C> recvBuffer(recvSize);
            mpi::SendRecv
            ( q21_MR_STAR.Buffer(), sendSize, transposeRank, 0,
              &recvBuffer[0],       recvSize, transposeRank, 0,
              g.VCComm() );

            // Unpack the recv buffer directly onto p21[MC,* ]
            C* p21_MC_STAR_Buffer = p21_MC_STAR.Buffer();
            for( int i=0; i<recvSize; ++i )
                p21_MC_STAR_Buffer[i] += recvBuffer[i];
        }

        if( W22.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p21[MC,* ] to the next process column.
            const int a21LocalHeight = a21.LocalHeight();

            const int nextProcessRow = (alpha11.ColAlignment()+1) % r;
            const int nextProcessCol = (alpha11.RowAlignment()+1) % r;

            std::vector<C> reduceToOneRecvBuffer(a21LocalHeight);
            mpi::Reduce
            ( p21_MC_STAR.Buffer(),
              &reduceToOneRecvBuffer[0],
              a21LocalHeight, mpi::SUM, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Finish computing w21. During its computation, ensure that 
                // every process has a copy of the first element of the w21.
                // We know a priori that the first element of a21 is one.
                const C* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
                C myDotProduct = blas::Dot
                    ( a21LocalHeight, &reduceToOneRecvBuffer[0], 1,
                                      &a21_MC_STAR_Buffer[0],    1 );
                C sendBuffer[2], recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.Row()==nextProcessRow ?
                                  reduceToOneRecvBuffer[0] : 0 );
                mpi::AllReduce
                ( sendBuffer, recvBuffer, 2, mpi::SUM, g.ColComm() );
                C dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w21LastBuffer
                // - w21LastFirstEntry
                C scale = dotProduct*Conj(tau) / C(2);
                for( int i=0; i<a21LocalHeight; ++i )
                    w21LastBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a21_MC_STAR_Buffer[i] );
                w21LastFirstEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, so we just need to form w21[MC,* ]
            // and w21[MR,* ]. 
            const int a21LocalHeight = a21.LocalHeight();

            // AllReduce sum p21[MC,* ] over process rows
            std::vector<C> allReduceRecvBuffer(a21LocalHeight);
            mpi::AllReduce
            ( p21_MC_STAR.Buffer(),
              &allReduceRecvBuffer[0],
              a21LocalHeight, mpi::SUM, g.RowComm() );

            // Finish computing w21.
            C myDotProduct = blas::Dot
                ( a21LocalHeight, &allReduceRecvBuffer[0], 1,
                                  a21_MC_STAR.Buffer(),    1 );
            C dotProduct;
            mpi::AllReduce
            ( &myDotProduct, &dotProduct, 1, mpi::SUM, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<C,MC,STAR> w21_MC_STAR(g);
            DistMatrix<C,MR,STAR> w21_MR_STAR(g);
            View
            ( w21_MC_STAR,
              W_MC_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );
            View
            ( w21_MR_STAR,
              W_MR_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );

            // Store w21[MC,* ]
            C scale = dotProduct*Conj(tau) / C(2);
            C* w21_MC_STAR_Buffer = w21_MC_STAR.Buffer();
            const C* a21_MC_STAR_Buffer = a21_MC_STAR.Buffer();
            for( int i=0; i<a21LocalHeight; ++i )
                w21_MC_STAR_Buffer[i] = 
                    tau*( allReduceRecvBuffer[i]-scale*a21_MC_STAR_Buffer[i] );

            // Fast transpose w21[MC,* ] -> w21[MR,* ]
            if( onDiagonal )
            {
                MemCopy
                ( w21_MR_STAR.Buffer(),
                  w21_MC_STAR.Buffer(), a21LocalHeight );
            }
            else
            {
                // Pairwise exchange with the transpose process
                const int sendSize = A22.LocalHeight();
                const int recvSize = A22.LocalWidth();
                mpi::SendRecv
                ( w21_MC_STAR.Buffer(), sendSize, transposeRank, 0,
                  w21_MR_STAR.Buffer(), recvSize, transposeRank, 0,
                  g.VCComm() );
            }
        }
        //--------------------------------------------------------------------//
        a21_MC_STAR.FreeAlignments();
        a21_MR_STAR.FreeAlignments();
        p21_MC_STAR.FreeAlignments();
        p21_MR_STAR.FreeAlignments();
        q21_MR_STAR.FreeAlignments();
        x01_MR_STAR.FreeAlignments();
        y01_MR_STAR.FreeAlignments();

        SlidePartitionDown
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );

        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        SlidePartitionDownDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        firstIteration = false;
    }
    PopBlocksizeStack();

    // View the portion of A that e is the subdiagonal of, then place e into it
    DistMatrix<C> expandedATL(g);
    View( expandedATL, A, 0, 0, panelSize+1, panelSize+1 );
    expandedATL.SetRealPartOfDiagonal( e, -1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace hermitian_tridiag
} // namespace elem

#endif // ifndef LAPACK_HERMITIANTRIDIAG_PANELLSQUARE_HPP
