/*
   Copyright (c) 2009-2012, Jack Poulson
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

namespace elem {
namespace internal {

template<typename R> 
inline void
HermitianPanelTridiagUSquare
( DistMatrix<R>& A,
  DistMatrix<R>& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR, 
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR )
{
    const int panelSize = W.Width();
    const int topSize = W.Height()-panelSize;
#ifndef RELEASE
    PushCallStack("internal::HermitianPanelTridiagUSquare");
    if( A.Grid() != W.Grid() )
        throw std::logic_error
        ("A and W must be distributed over the same grid");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square");
    if( A.Height() != W.Height() )
        throw std::logic_error("A and W must be the same height");
    if( W.Height() < panelSize )
        throw std::logic_error("W must be a column panel");
#endif
    const Grid& g = A.Grid();
    const int r = g.Height();

    // Find the process holding our transposed data
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

    // Create a distributed matrix for storing the superdiagonal
    DistMatrix<R,MD,STAR> e(g);
    DistMatrix<R> expandedABR(g);
    expandedABR.View( A, topSize-1, topSize-1, panelSize+1, panelSize+1 );
    e.AlignWithDiagonal( expandedABR, 1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<R> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g), a01T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),           alpha01B(g),
                         A20(g), a21(g),     A22(g),  A02T(g), A00Pan(g);
    DistMatrix<R> 
        WTL(g), WTR(g),  W00(g), w01(g),     W02(g),  WCol(g),
        WBL(g), WBR(g),  w10(g), omega11(g), w12(g),
                         W20(g), w21(g),     W22(g),  W02T(g), w01Last(g);
    DistMatrix<R,MD,STAR> eT(g),  e0(g),
                          eB(g),  epsilon1(g),
                                  e2(g);

    // Temporary distributions
    std::vector<R> w01LastLocalBuffer(A.Height()/r+1);
    DistMatrix<R,MC,STAR> a01_MC_STAR(g), a01T_MC_STAR(g), a01Last_MC_STAR(g);
    DistMatrix<R,MR,STAR> a01_MR_STAR(g), a01Last_MR_STAR(g);
    DistMatrix<R,MC,STAR> p01_MC_STAR(g), p01T_MC_STAR(g);
    DistMatrix<R,MR,STAR> p01_MR_STAR(g);
    DistMatrix<R,MR,STAR> q01_MR_STAR(g);
    DistMatrix<R,MR,STAR> x21_MR_STAR(g);
    DistMatrix<R,MR,STAR> y21_MR_STAR(g);
    DistMatrix<R,MC,STAR> w01Last_MC_STAR(g);
    DistMatrix<R,MR,STAR> w01Last_MR_STAR(g);

    // Push to the blocksize of 1, then pop at the end of the routine
    PushBlocksizeStack( 1 );

    PartitionUpRightDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUpRightDiagonal
    ( W, WTL, WTR,
         WBL, WBR, 0 );
    PartitionUp
    ( e, eT,
         eB, 0 );
    bool firstIteration = true;
    R tau = 0;
    R w01LastBottomEntry = 0;
    while( WBR.Width() < panelSize )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12, 
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
        
        RepartitionUpDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );

        RepartitionUp
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        ACol.View2x1
        ( a01,
          alpha11 );

        WCol.View2x1
        ( w01,
          omega11 );

        // View the portions of A02 and W02 outside of this panel's square
        A02T.View( A02, 0, 0, topSize, A02.Width() );
        W02T.View( W02, 0, 0, topSize, W02.Width() );

        // View the portion of A00 that is inside this panel
        A00Pan.View( A00, 0, topSize, A00.Height(), A00.Width()-topSize );

        if( !firstIteration )
        {
            a01Last_MC_STAR.View
            ( APan_MC_STAR, 0, WTL.Width(), ACol.Height(), 1 );
            a01Last_MR_STAR.View
            ( APan_MR_STAR, 0, WTL.Width(), ACol.Height(), 1 );
            w01Last.View
            ( W, 0, WTL.Width(), ACol.Height(), 1 );
        }
            
        PartitionUp
        ( a01, a01T,
               alpha01B, 1 );

        a01_MC_STAR.AlignWith( A00 );
        a01_MR_STAR.AlignWith( A00 );
        p01_MC_STAR.AlignWith( A00 );
        p01_MR_STAR.AlignWith( A00 );
        q01_MR_STAR.AlignWith( A00 );
        x21_MR_STAR.AlignWith( A02T );
        y21_MR_STAR.AlignWith( A02T );

        a01_MC_STAR.ResizeTo( a01.Height(), 1 );
        a01_MR_STAR.ResizeTo( a01.Height(), 1 );
        p01_MC_STAR.ResizeTo( a01.Height(), 1 );
        p01_MR_STAR.ResizeTo( a01.Height(), 1 );
        q01_MR_STAR.ResizeTo( a01.Height(), 1 );
        x21_MR_STAR.ResizeTo( A02.Width(), 1 );
        y21_MR_STAR.ResizeTo( A02.Width(), 1 );

        // View the portions of a01[MC,* ] and p01[MC,* ] above the current
        // panel's square
        a01T_MC_STAR.View( a01_MC_STAR, 0, 0, topSize, 1 );
        p01T_MC_STAR.View( p01_MC_STAR, 0, 0, topSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const int AColLocalHeight = ACol.LocalHeight();
                R* AColLocalBuffer = ACol.LocalBuffer();
                const R* a01Last_MC_STAR_LocalBuffer = 
                    a01Last_MC_STAR.LocalBuffer();
                for( int i=0; i<AColLocalHeight; ++i )
                    AColLocalBuffer[i] -=
                        w01LastLocalBuffer[i] + 
                        a01Last_MC_STAR_LocalBuffer[i]*w01LastBottomEntry;
            }
        }
        if( thisIsMyCol )
        {
            // Compute the Householder reflector
            tau = ColReflector( alpha01B, a01T );
        }
        // Store the subdiagonal value and turn a01 into a proper scaled 
        // reflector by explicitly placing the implicit one in its bottom entry
        alpha01B.GetDiagonal( epsilon1 );
        alpha01B.Set( 0, 0, (R)1 );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a01 within its process row. 
        // Otherwise, also add w01 into the broadcast.
        if( firstIteration )
        {
            const int a01LocalHeight = a01.LocalHeight();
            std::vector<R> rowBroadcastBuffer(a01LocalHeight+1);
            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a01 and tau
                MemCopy
                ( &rowBroadcastBuffer[0], a01.LocalBuffer(), a01LocalHeight );
                rowBroadcastBuffer[a01LocalHeight] = tau;
            }
            // Broadcast a01 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a01LocalHeight+1, a01.RowAlignment(), g.RowComm() );
            // Store a01[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy
            ( a01_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0], a01LocalHeight );
            // Store a01[MC,* ] into APan[MC,* ]
            MemCopy
            ( APan_MC_STAR.LocalBuffer(0,W00.Width()), 
              &rowBroadcastBuffer[0], a01LocalHeight );
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight];
            
            // Take advantage of the square process grid in order to form
            // a01[MR,* ] from a01[MC,* ]
            if( onDiagonal )
            {
                MemCopy
                ( a01_MR_STAR.LocalBuffer(),
                  a01_MC_STAR.LocalBuffer(), a01LocalHeight );
            }
            else
            {
                // Pairwise exchange
                const int sendSize = A00.LocalHeight();
                const int recvSize = A00.LocalWidth();
                mpi::SendRecv
                ( a01_MC_STAR.LocalBuffer(), sendSize, transposeRank, 0,
                  a01_MR_STAR.LocalBuffer(), recvSize, transposeRank, 0,
                  g.VCComm() );
            }
            // Store a01[MR,* ]
            MemCopy
            ( APan_MR_STAR.LocalBuffer(0,W00.Width()),
              a01_MR_STAR.LocalBuffer(),
              a01_MR_STAR.LocalHeight() );
        }
        else
        {
            const int a01LocalHeight = a01.LocalHeight();
            const int w01LastLocalHeight = ACol.LocalHeight();
            std::vector<R> 
                rowBroadcastBuffer(a01LocalHeight+w01LastLocalHeight+1);
            if( thisIsMyCol ) 
            {
                // Pack the broadcast buffer with a01, w01Last, and tau
                MemCopy
                ( &rowBroadcastBuffer[0], 
                  a01.LocalBuffer(), a01LocalHeight );
                MemCopy
                ( &rowBroadcastBuffer[a01LocalHeight], 
                  &w01LastLocalBuffer[0], w01LastLocalHeight );
                rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight] = tau;
            }
            // Broadcast a01, w01Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a01LocalHeight+w01LastLocalHeight+1, 
              a01.RowAlignment(), g.RowComm() );
            // Store a01[MC,* ] into its DistMatrix class 
            MemCopy
            ( a01_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0], a01LocalHeight );
            // Store a01[MC,* ] into APan[MC,* ]
            MemCopy
            ( APan_MC_STAR.LocalBuffer(0,W00.Width()),
              &rowBroadcastBuffer[0], a01LocalHeight );
            // Store w01Last[MC,* ] into its DistMatrix class
            w01Last_MC_STAR.AlignWith( A00 );
            w01Last_MC_STAR.ResizeTo( a01.Height()+1, 1 );
            MemCopy
            ( w01Last_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            // Store the bottom part of w01Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w01.
            MemCopy
            ( W_MC_STAR.LocalBuffer(0,W00.Width()+1),
              &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            if( g.Col() == w01Last.RowAlignment() )
            {
                MemCopy
                ( w01Last.LocalBuffer(),
                  &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            }
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight];

            // Take advantage of the square process grid in order to quickly
            // form a01[MR,* ] and w01Last[MR,* ] from their [MC,* ] 
            // counterparts
            w01Last_MR_STAR.AlignWith( A00 );
            w01Last_MR_STAR.ResizeTo( w01Last.Height(), 1 );
            if( onDiagonal )
            {
                MemCopy
                ( a01_MR_STAR.LocalBuffer(),
                  a01_MC_STAR.LocalBuffer(), a01LocalHeight );
                MemCopy
                ( w01Last_MR_STAR.LocalBuffer(),
                  w01Last_MC_STAR.LocalBuffer(), w01LastLocalHeight );
            }
            else
            {
                const int sendSize = A00.LocalHeight()+ATL.LocalHeight();
                const int recvSize = A00.LocalWidth()+ATL.LocalWidth();
                std::vector<R> sendBuffer(sendSize);
                std::vector<R> recvBuffer(recvSize);

                // Pack the send buffer
                MemCopy
                ( &sendBuffer[0],
                  a01_MC_STAR.LocalBuffer(), A00.LocalHeight() );
                MemCopy
                ( &sendBuffer[A00.LocalHeight()],
                  w01Last_MC_STAR.LocalBuffer(), ATL.LocalHeight() );

                // Pairwise exchange
                mpi::SendRecv
                ( &sendBuffer[0], sendSize, transposeRank, 0,
                  &recvBuffer[0], recvSize, transposeRank, 0,
                  g.VCComm() );

                // Unpack the recv buffer
                MemCopy
                ( a01_MR_STAR.LocalBuffer(),
                  &recvBuffer[0], A00.LocalWidth() );
                MemCopy
                ( w01Last_MR_STAR.LocalBuffer(),
                  &recvBuffer[A00.LocalWidth()], ATL.LocalWidth() );
            }

            // Store w01Last[MR,* ]
            MemCopy
            ( W_MR_STAR.LocalBuffer(0,W00.Width()+1),
              w01Last_MR_STAR.LocalBuffer(),
              w01Last_MR_STAR.LocalHeight() );
            // Store a01[MR,* ]
            MemCopy
            ( APan_MR_STAR.LocalBuffer(0,W00.Width()),
              a01_MR_STAR.LocalBuffer(),
              a01_MR_STAR.LocalHeight() );

            // Update the portion of A00 that is in our current panel with 
            // w01Last and a01Last using two gers. We do not need their bottom
            // entries. We trash the lower triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            DistMatrix<R,MC,STAR> a01Last_MC_STAR_Top(g),
                                  w01Last_MC_STAR_Top(g);
            DistMatrix<R,MR,STAR> a01Last_MR_STAR_TopPan(g),
                                  w01Last_MR_STAR_TopPan(g);
            a01Last_MC_STAR_Top.View( a01Last_MC_STAR, 0, 0, a01.Height(), 1 );
            w01Last_MC_STAR_Top.View( w01Last_MC_STAR, 0, 0, a01.Height(), 1 );
            a01Last_MR_STAR_TopPan.View
            ( a01Last_MR_STAR, topSize, 0, a01.Height()-topSize, 1 );
            w01Last_MR_STAR_TopPan.View
            ( w01Last_MR_STAR, topSize, 0, a01.Height()-topSize, 1 );
            const R* a01_MC_STAR_Buffer = a01Last_MC_STAR_Top.LocalBuffer();
            const R* w01_MC_STAR_Buffer = w01Last_MC_STAR_Top.LocalBuffer();
            const R* a01_MR_STAR_Buffer = a01Last_MR_STAR_TopPan.LocalBuffer();
            const R* w01_MR_STAR_Buffer = w01Last_MR_STAR_TopPan.LocalBuffer();
            R* A00PanBuffer = A00Pan.LocalBuffer();
            const int localHeight = A00Pan.LocalHeight();
            const int localWidth = A00Pan.LocalWidth();
            const int lDim = A00Pan.LocalLDim();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    A00PanBuffer[iLocal+jLocal*lDim] -=
                        w01_MC_STAR_Buffer[iLocal]*a01_MR_STAR_Buffer[jLocal] +
                        a01_MC_STAR_Buffer[iLocal]*w01_MR_STAR_Buffer[jLocal];

            // We are through with the last iteration's w01
            w01Last_MC_STAR.FreeAlignments();
            w01Last_MR_STAR.FreeAlignments();
        }

        // Form the local portions of (A00 a01) into p01[MC,* ] and q01[MR,* ]:
        //   p01[MC,* ] := triu(A00)[MC,MR] a01[MR,* ]
        //   q01[MR,* ] := triu(A00,+1)'[MR,MC] a01[MC,* ]
        if( A00.ColShift() < A00.RowShift() )
        {
            // We are above the diagonal, so we can multiply without an 
            // offset for triu(A00)[MC,MR] and triu(A00,1)'[MR,MC]
            if( A00.LocalHeight() != 0 )
            {
                // Our local portion of p01[MC,* ] might be one entry longer
                // than A00.LocalWidth(), so go ahead and set the last entry
                // to 0.
                R* p01_MC_STAR_LocalBuffer = p01_MC_STAR.LocalBuffer();
                p01_MC_STAR_LocalBuffer[A00.LocalHeight()-1] = 0;
                MemCopy
                ( p01_MC_STAR.LocalBuffer(),
                  a01_MR_STAR.LocalBuffer(), A00.LocalWidth() );
                blas::Trmv
                ( 'U', 'N', 'N', A00.LocalWidth(),
                  A00.LocalBuffer(), A00.LocalLDim(),
                  p01_MC_STAR.LocalBuffer(), 1 );
            }
            if( A00.LocalWidth() != 0 )
            {
                MemCopy
                ( q01_MR_STAR.LocalBuffer(),
                  a01_MC_STAR.LocalBuffer(), A00.LocalWidth() );
                blas::Trmv
                ( 'U', 'T', 'N', A00.LocalWidth(),
                  A00.LocalBuffer(), A00.LocalLDim(),
                  q01_MR_STAR.LocalBuffer(), 1 );
            }
        }
        else if( A00.ColShift() > A00.RowShift() )
        {
            // We are below the diagonal, so we need to use an offset 
            // for both triu(A00)[MC,MR] and triu(A00,+1)'[MR,MC]
            const R* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
            const R* a01_MR_STAR_LocalBuffer = a01_MR_STAR.LocalBuffer();
            const R* A00LocalBuffer = A00.LocalBuffer();
            if( A00.LocalHeight() != 0 )
            {
                // The last entry of p01[MC,* ] will be zero due to the forced
                // offset
                R* p01_MC_STAR_LocalBuffer = p01_MC_STAR.LocalBuffer();
                p01_MC_STAR_LocalBuffer[A00.LocalHeight()-1] = 0;
                MemCopy
                ( &p01_MC_STAR_LocalBuffer[0],
                  &a01_MR_STAR_LocalBuffer[1], A00.LocalWidth()-1 );
                blas::Trmv
                ( 'U', 'N', 'N', A00.LocalWidth()-1,
                  &A00LocalBuffer[A00.LocalLDim()], A00.LocalLDim(),
                  &p01_MC_STAR_LocalBuffer[0], 1 );
            }
            if( A00.LocalWidth() != 0 )
            {
                // The first entry of q01[MR,* ] will be zero due to the forced
                // offset
                R* q01_MR_STAR_LocalBuffer = q01_MR_STAR.LocalBuffer();
                q01_MR_STAR_LocalBuffer[0] = 0;
                MemCopy
                ( &q01_MR_STAR_LocalBuffer[1],
                  &a01_MC_STAR_LocalBuffer[0], A00.LocalWidth()-1 );
                blas::Trmv
                ( 'U', 'T', 'N', A00.LocalWidth()-1,
                  &A00LocalBuffer[A00.LocalLDim()], A00.LocalLDim(),
                  &q01_MR_STAR_LocalBuffer[1], 1 );
            }
        }
        else
        {
            // We are on the diagonal, so we only need an offset for
            // triu(A00,+1)'[MR,MC]
            if( A00.LocalWidth() != 0 )
            {
                MemCopy
                ( p01_MC_STAR.LocalBuffer(),
                  a01_MR_STAR.LocalBuffer(), A00.LocalHeight() );
                blas::Trmv
                ( 'U', 'N', 'N', A00.LocalHeight(),
                  A00.LocalBuffer(), A00.LocalLDim(),
                  p01_MC_STAR.LocalBuffer(), 1 );

                // The first entry of q01[MR,* ] must be zero due to the forced
                // offset
                const R* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
                const R* A00LocalBuffer = A00.LocalBuffer();
                R* q01_MR_STAR_LocalBuffer = q01_MR_STAR.LocalBuffer();
                q01_MR_STAR_LocalBuffer[0] = 0;
                MemCopy
                ( &q01_MR_STAR_LocalBuffer[1],
                  &a01_MC_STAR_LocalBuffer[0], A00.LocalWidth()-1 );
                blas::Trmv
                ( 'U', 'T', 'N', A00.LocalWidth()-1,
                  &A00LocalBuffer[A00.LocalLDim()], A00.LocalLDim(),
                  &q01_MR_STAR_LocalBuffer[1], 1 );
            }
        }

        LocalGemv( TRANSPOSE, (R)1, W02T, a01T_MC_STAR, (R)0, x21_MR_STAR );
        LocalGemv( TRANSPOSE, (R)1, A02T, a01T_MC_STAR, (R)0, y21_MR_STAR );

        // Combine the AllReduce column summations of x21[MR,* ] and y21[MR,* ]
        {
            const int x21LocalHeight = x21_MR_STAR.LocalHeight();
            const int y21LocalHeight = y21_MR_STAR.LocalHeight();
            const int reduceSize = x21LocalHeight+y21LocalHeight;
            std::vector<R> colSumSendBuffer(reduceSize);
            std::vector<R> colSumRecvBuffer(reduceSize);
            MemCopy
            ( &colSumSendBuffer[0], 
              x21_MR_STAR.LocalBuffer(), x21LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x21LocalHeight],
              y21_MR_STAR.LocalBuffer(), y21LocalHeight );
            mpi::AllReduce
            ( &colSumSendBuffer[0], 
              &colSumRecvBuffer[0],
              reduceSize, mpi::SUM, g.ColComm() );
            MemCopy
            ( x21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[0], x21LocalHeight );
            MemCopy
            ( y21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[x21LocalHeight], y21LocalHeight );
        }

        LocalGemv( NORMAL, (R)-1, A02T, x21_MR_STAR, (R)1, p01T_MC_STAR );
        LocalGemv( NORMAL, (R)-1, W02T, y21_MR_STAR, (R)1, p01T_MC_STAR );

        // Fast transpose the unsummed q01[MR,* ] -> q01[MC,* ], so that
        // it needs to be summed over process rows instead of process
        // columns. We immediately add it onto p01[MC,* ], which also needs to
        // be summed within process rows.
        if( onDiagonal )
        {
            const int a01LocalHeight = a01.LocalHeight();
            R* p01_MC_STAR_LocalBuffer = p01_MC_STAR.LocalBuffer();
            const R* q01_MR_STAR_LocalBuffer = q01_MR_STAR.LocalBuffer();
            for( int i=0; i<a01LocalHeight; ++i )
                p01_MC_STAR_LocalBuffer[i] += q01_MR_STAR_LocalBuffer[i];
        }
        else
        {
            // Pairwise exchange with the transpose process
            const int sendSize = A00.LocalWidth();
            const int recvSize = A00.LocalHeight();
            std::vector<R> recvBuffer(recvSize);
            mpi::SendRecv
            ( q01_MR_STAR.LocalBuffer(), sendSize, transposeRank, 0,
              &recvBuffer[0],            recvSize, transposeRank, 0,
              g.VCComm() );

            // Unpack the recv buffer directly onto p01[MC,* ]
            R* p01_MC_STAR_LocalBuffer = p01_MC_STAR.LocalBuffer();
            for( int i=0; i<recvSize; ++i )
                p01_MC_STAR_LocalBuffer[i] += recvBuffer[i];
        }

        if( W00.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p01[MC,* ] to the next process column.
            const int a01LocalHeight = a01.LocalHeight();

            const int nextProcessRow = (alpha11.ColAlignment()+r-1) % r;
            const int nextProcessCol = (alpha11.RowAlignment()+r-1) % r;

            std::vector<R> reduceToOneRecvBuffer(a01LocalHeight);
            mpi::Reduce
            ( p01_MC_STAR.LocalBuffer(),
              &reduceToOneRecvBuffer[0],
              a01LocalHeight, mpi::SUM, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Finish computing w01. During its computation, ensure that 
                // every process has a copy of the bottom element of the w01.
                // We know a priori that the bottom element of a01 is one.
                const R* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
                R myDotProduct = blas::Dot
                    ( a01LocalHeight, &reduceToOneRecvBuffer[0],   1, 
                                      &a01_MC_STAR_LocalBuffer[0], 1 );
                R sendBuffer[2];
                R recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.Row()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[a01LocalHeight-1] : 0 );
                mpi::AllReduce
                ( sendBuffer, recvBuffer, 2, mpi::SUM, g.ColComm() );
                R dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w01LastLocalBuffer
                // - w01LastBottomEntry
                R scale = 0.5*dotProduct*tau;
                for( int i=0; i<a01LocalHeight; ++i )
                    w01LastLocalBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a01_MC_STAR_LocalBuffer[i] );
                w01LastBottomEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, so we just need to form w01[MC,* ]
            // and w01[MR,* ].
            const int a01LocalHeight = a01.LocalHeight();

            // AllReduce sum p01[MC,* ] over process rows
            std::vector<R> allReduceRecvBuffer(a01LocalHeight);
            mpi::AllReduce
            ( p01_MC_STAR.LocalBuffer(),
              &allReduceRecvBuffer[0],
              a01LocalHeight, mpi::SUM, g.RowComm() );

            // Finish computing w01. 
            const R* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
            R myDotProduct = blas::Dot
                ( a01LocalHeight, &allReduceRecvBuffer[0],     1, 
                                  &a01_MC_STAR_LocalBuffer[0], 1 );
            R dotProduct;
            mpi::AllReduce
            ( &myDotProduct, &dotProduct, 1, mpi::SUM, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<R,MC,STAR> w01_MC_STAR(g);
            DistMatrix<R,MR,STAR> w01_MR_STAR(g);
            w01_MC_STAR.View( W_MC_STAR, 0, W00.Width(), w01.Height(), 1 );
            w01_MR_STAR.View( W_MR_STAR, 0, W00.Width(), w01.Height(), 1 );

            // Store w01[MC,* ]
            R scale = 0.5*dotProduct*tau;
            R* w01_MC_STAR_LocalBuffer = w01_MC_STAR.LocalBuffer();
            for( int i=0; i<a01LocalHeight; ++i )
                w01_MC_STAR_LocalBuffer[i] = tau*
                    ( allReduceRecvBuffer[i]-
                      scale*a01_MC_STAR_LocalBuffer[i] );

            // Fast transpose w01[MC,* ] -> w01[MR,* ]
            if( onDiagonal )
            {
                MemCopy
                ( w01_MR_STAR.LocalBuffer(),
                  w01_MC_STAR.LocalBuffer(), a01LocalHeight );
            }
            else
            {
                // Pairwise exchange with the transpose process
                const int sendSize = A00.LocalHeight();
                const int recvSize = A00.LocalWidth();
                mpi::SendRecv
                ( w01_MC_STAR.LocalBuffer(), sendSize, transposeRank, 0,
                  w01_MR_STAR.LocalBuffer(), recvSize, transposeRank, 0,
                  g.VCComm() );
            }
        }
        //--------------------------------------------------------------------//
        a01_MC_STAR.FreeAlignments();
        a01_MR_STAR.FreeAlignments();
        p01_MC_STAR.FreeAlignments();
        p01_MR_STAR.FreeAlignments();
        q01_MR_STAR.FreeAlignments();
        x21_MR_STAR.FreeAlignments();
        y21_MR_STAR.FreeAlignments();
        
        SlidePartitionUp
        ( eT,  e0,
         /**/ /********/
               epsilon1,
          eB,  e2 );
 
        SlidePartitionUpDiagonal
        ( WTL, /**/ WTR,  W00, /**/ w01,     W02,
         /*************/ /**********************/
               /**/       w10, /**/ omega11, w12,
          WBL, /**/ WBR,  W20, /**/ w21,     W22 );

        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

       firstIteration = false;
    }
    PopBlocksizeStack();

    expandedABR.SetDiagonal( e, 1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
HermitianPanelTridiagUSquare
( DistMatrix<Complex<R> >& A,
  DistMatrix<Complex<R> >& W,
  DistMatrix<Complex<R>,MD,STAR>& t,
  DistMatrix<Complex<R>,MC,STAR>& APan_MC_STAR, 
  DistMatrix<Complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<Complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<Complex<R>,MR,STAR>& W_MR_STAR )
{
    const int panelSize = W.Width();
    const int topSize = W.Height()-panelSize;
#ifndef RELEASE
    PushCallStack("internal::HermitianPanelTridiagUSquare");
    if( A.Grid() != W.Grid() || W.Grid() != t.Grid() )
        throw std::logic_error
        ("A, W, and t must be distributed over the same grid.");
    if( A.Height() != A.Width() )
        throw std::logic_error("A must be square.");
    if( A.Height() != W.Height() )
        throw std::logic_error( "A and W must be the same height.");
    if( W.Height() < panelSize )
        throw std::logic_error("W must be a column panel.");
    if( t.Height() != W.Width() || t.Width() != 1 )
        throw std::logic_error
        ("t must be a column vector of the same length as W's width.");
#endif
    typedef Complex<R> C;

    const Grid& g = A.Grid();
    const int r = g.Height();

    // Find the process holding our transposed data
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

    // Create a distributed matrix for storing the superdiagonal
    DistMatrix<R,MD,STAR> e(g);
    DistMatrix<C> expandedABR(g);
    expandedABR.View( A, topSize-1, topSize-1, panelSize+1, panelSize+1 );
    e.AlignWithDiagonal( expandedABR, 1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<C> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g), a01T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),           alpha01B(g),
                         A20(g), a21(g),     A22(g),  A02T(g), A00Pan(g);
    DistMatrix<C> 
        WTL(g), WTR(g),  W00(g), w01(g),     W02(g),  WCol(g),
        WBL(g), WBR(g),  w10(g), omega11(g), w12(g),
                         W20(g), w21(g),     W22(g),  W02T(g), w01Last(g);
    DistMatrix<R,MD,STAR> eT(g),  e0(g),
                          eB(g),  epsilon1(g),
                                  e2(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  tau1(g),
                t2(g);

    // Temporary distributions
    std::vector<C> w01LastLocalBuffer(A.Height()/r+1);
    DistMatrix<C,MC,STAR> a01_MC_STAR(g), a01T_MC_STAR(g), a01Last_MC_STAR(g);
    DistMatrix<C,MR,STAR> a01_MR_STAR(g), a01Last_MR_STAR(g);
    DistMatrix<C,MC,STAR> p01_MC_STAR(g), p01T_MC_STAR(g);
    DistMatrix<C,MR,STAR> p01_MR_STAR(g);
    DistMatrix<C,MR,STAR> q01_MR_STAR(g);
    DistMatrix<C,MR,STAR> x21_MR_STAR(g);
    DistMatrix<C,MR,STAR> y21_MR_STAR(g);
    DistMatrix<C,MC,STAR> w01Last_MC_STAR(g);
    DistMatrix<C,MR,STAR> w01Last_MR_STAR(g);

    // Push to the blocksize of 1, then pop at the end of the routine
    PushBlocksizeStack( 1 );

    PartitionUpRightDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionUpRightDiagonal
    ( W, WTL, WTR,
         WBL, WBR, 0 );
    PartitionUp
    ( e, eT,
         eB, 0 );
    PartitionUp
    ( t, tT,
         tB, 0 );
    bool firstIteration = true;
    C tau = 0;
    C w01LastBottomEntry = 0;
    while( WBR.Width() < panelSize )
    {
        RepartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12, 
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
        
        RepartitionUpDiagonal
        ( WTL, /**/ WTR,  W00, w01,     /**/ W02,
               /**/       w10, omega11, /**/ w12,
         /*************/ /**********************/
          WBL, /**/ WBR,  W20, w21,     /**/ W22 );

        RepartitionUp
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );

        RepartitionUp
        ( tT,  t0,
               tau1,
         /**/ /****/
          tB,  t2 );

        ACol.View2x1
        ( a01,
          alpha11 );

        WCol.View2x1
        ( w01,
          omega11 );

        // View the portions of A02 and W0T outside of this panel's square
        A02T.View( A02, 0, 0, topSize, A02.Width() );
        W02T.View( W02, 0, 0, topSize, W02.Width() );

        // View the portion of A00 inside the current panel
        A00Pan.View( A00, 0, topSize, A00.Height(), A00.Width()-topSize );

        if( !firstIteration )
        {
            a01Last_MC_STAR.View
            ( APan_MC_STAR, 0, WTL.Width(), ACol.Height(), 1 );
            a01Last_MR_STAR.View
            ( APan_MR_STAR, 0, WTL.Width(), ACol.Height(), 1 );
            w01Last.View
            ( W, 0, WTL.Width(), ACol.Height(), 1 );
        }
            
        PartitionUp
        ( a01, a01T,
               alpha01B, 1 );

        a01_MC_STAR.AlignWith( A00 );
        a01_MR_STAR.AlignWith( A00 );
        p01_MC_STAR.AlignWith( A00 );
        p01_MR_STAR.AlignWith( A00 );
        q01_MR_STAR.AlignWith( A00 );
        x21_MR_STAR.AlignWith( A02T );
        y21_MR_STAR.AlignWith( A02T );
        
        a01_MC_STAR.ResizeTo( a01.Height(), 1 );
        a01_MR_STAR.ResizeTo( a01.Height(), 1 );
        p01_MC_STAR.ResizeTo( a01.Height(), 1 );
        p01_MR_STAR.ResizeTo( a01.Height(), 1 );
        q01_MR_STAR.ResizeTo( a01.Height(), 1 );
        x21_MR_STAR.ResizeTo( A02.Width(), 1 );
        y21_MR_STAR.ResizeTo( A02.Width(), 1 );

        // View the portions of a01[MC,* ] and p01[MC,* ] above the current
        // panel's square
        a01T_MC_STAR.View( a01_MC_STAR, 0, 0, topSize, 1 );
        p01T_MC_STAR.View( p01_MC_STAR, 0, 0, topSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.Col() == alpha11.RowAlignment() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                const int AColLocalHeight = ACol.LocalHeight();
                C* AColLocalBuffer = ACol.LocalBuffer();
                const C* a01Last_MC_STAR_LocalBuffer = 
                    a01Last_MC_STAR.LocalBuffer();
                for( int i=0; i<AColLocalHeight; ++i )
                    AColLocalBuffer[i] -=
                        w01LastLocalBuffer[i] + 
                        a01Last_MC_STAR_LocalBuffer[i]*Conj(w01LastBottomEntry);
            }
        }
        if( thisIsMyCol )
        {
            // Compute the Householder reflector
            tau = ColReflector( alpha01B, a01T );
            if( g.Row() == alpha01B.ColAlignment() )
                tau1.SetLocal(0,0,tau);
        }
        // Store the subdiagonal value and turn a01 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        alpha01B.GetRealPartOfDiagonal( epsilon1 );
        alpha01B.Set( 0, 0, (C)1 );

        // If this is the first iteration, have each member of the owning 
        // process column broadcast tau and a01 within its process row. 
        // Otherwise, also add w01 into the broadcast.
        if( firstIteration )
        {
            const int a01LocalHeight = a01.LocalHeight();
            std::vector<C> rowBroadcastBuffer(a01LocalHeight+1);
            if( thisIsMyCol )
            {
                // Pack the broadcast buffer with a01 and tau
                MemCopy
                ( &rowBroadcastBuffer[0], 
                  a01.LocalBuffer(), a01LocalHeight );
                rowBroadcastBuffer[a01LocalHeight] = tau;
            }
            // Broadcast a01 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a01LocalHeight+1, a01.RowAlignment(), g.RowComm() );
            // Store a01[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            MemCopy
            ( a01_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0], a01LocalHeight );
            // Store a01[MC,* ] into APan[MC,* ]
            MemCopy
            ( APan_MC_STAR.LocalBuffer(0,W00.Width()), 
              &rowBroadcastBuffer[0], a01LocalHeight );
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight];
            
            // Take advantage of the square process grid in order to form
            // a01[MR,* ] from a01[MC,* ]
            if( onDiagonal )
            {
                MemCopy
                ( a01_MR_STAR.LocalBuffer(),
                  a01_MC_STAR.LocalBuffer(), a01LocalHeight );
            }
            else
            {
                // Pairwise exchange
                const int sendSize = A00.LocalHeight();
                const int recvSize = A00.LocalWidth();
                mpi::SendRecv
                ( a01_MC_STAR.LocalBuffer(), sendSize, transposeRank, 0,
                  a01_MR_STAR.LocalBuffer(), recvSize, transposeRank, 0,
                  g.VCComm() );
            }
            // Store a01[MR,* ]
            MemCopy
            ( APan_MR_STAR.LocalBuffer(0,W00.Width()),
              a01_MR_STAR.LocalBuffer(), a01_MR_STAR.LocalHeight() );
        }
        else
        {
            const int a01LocalHeight = a01.LocalHeight();
            const int w01LastLocalHeight = ACol.LocalHeight();
            std::vector<C> 
                rowBroadcastBuffer(a01LocalHeight+w01LastLocalHeight+1);
            if( thisIsMyCol ) 
            {
                // Pack the broadcast buffer with a01, w01Last, and tau
                MemCopy
                ( &rowBroadcastBuffer[0], 
                  a01.LocalBuffer(), a01LocalHeight );
                MemCopy
                ( &rowBroadcastBuffer[a01LocalHeight], 
                  &w01LastLocalBuffer[0], w01LastLocalHeight );
                rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight] = tau;
            }
            // Broadcast a01, w01Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a01LocalHeight+w01LastLocalHeight+1, 
              a01.RowAlignment(), g.RowComm() );
            // Store a01[MC,* ] into its DistMatrix class 
            MemCopy
            ( a01_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0], a01LocalHeight );
            // Store a01[MC,* ] into APan[MC,* ]
            MemCopy
            ( APan_MC_STAR.LocalBuffer(0,W00.Width()), 
              &rowBroadcastBuffer[0], a01LocalHeight );
            // Store w01Last[MC,* ] into its DistMatrix class
            w01Last_MC_STAR.AlignWith( A00 );
            w01Last_MC_STAR.ResizeTo( a01.Height()+1, 1 );
            MemCopy
            ( w01Last_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            // Store the bottom part of w01Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w01.
            MemCopy
            ( W_MC_STAR.LocalBuffer(0,W00.Width()+1),
              &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            if( g.Col() == w01Last.RowAlignment() )
            {
                MemCopy
                ( w01Last.LocalBuffer(),
                  &rowBroadcastBuffer[a01LocalHeight], w01LastLocalHeight );
            }
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight];

            // Take advantage of the square process grid in order to quickly
            // form a01[MR,* ] and w01Last[MR,* ] from their [MC,* ]
            // counterparts
            w01Last_MR_STAR.AlignWith( A00 );
            w01Last_MR_STAR.ResizeTo( w01Last.Height(), 1 );
            if( onDiagonal )
            {
                MemCopy
                ( a01_MR_STAR.LocalBuffer(),
                  a01_MC_STAR.LocalBuffer(), a01LocalHeight );
                MemCopy
                ( w01Last_MR_STAR.LocalBuffer(),
                  w01Last_MC_STAR.LocalBuffer(), w01LastLocalHeight );
            }
            else
            {
                const int sendSize = A00.LocalHeight()+ATL.LocalHeight();
                const int recvSize = A00.LocalWidth()+ATL.LocalWidth();
                std::vector<C> sendBuffer(sendSize);
                std::vector<C> recvBuffer(recvSize);

                // Pack the send buffer
                MemCopy
                ( &sendBuffer[0],
                  a01_MC_STAR.LocalBuffer(), A00.LocalHeight() );
                MemCopy
                ( &sendBuffer[A00.LocalHeight()],
                  w01Last_MC_STAR.LocalBuffer(), ATL.LocalHeight() );

                // Pairwise exchange
                mpi::SendRecv
                ( &sendBuffer[0], sendSize, transposeRank, 0,
                  &recvBuffer[0], recvSize, transposeRank, 0,
                  g.VCComm() );

                // Unpack the recv buffer
                MemCopy
                ( a01_MR_STAR.LocalBuffer(),
                  &recvBuffer[0], A00.LocalWidth() );
                MemCopy
                ( w01Last_MR_STAR.LocalBuffer(),
                  &recvBuffer[A00.LocalWidth()], ATL.LocalWidth() );
            }

            // Store w01Last[MR,* ]
            MemCopy
            ( W_MR_STAR.LocalBuffer(0,W00.Width()+1),
              w01Last_MR_STAR.LocalBuffer(),
              w01Last_MR_STAR.LocalHeight() );
            // Store a01[MR,* ]
            MemCopy
            ( APan_MR_STAR.LocalBuffer(0,W00.Width()),
              a01_MR_STAR.LocalBuffer(), a01_MR_STAR.LocalHeight() );

            // Update the portion of A00 that is in our current panel with 
            // w01Last and a01Last using two gers. We do not need their bottom
            // entries. We trash the lower triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            DistMatrix<C,MC,STAR> a01Last_MC_STAR_Top(g),
                                  w01Last_MC_STAR_Top(g);
            DistMatrix<C,MR,STAR> a01Last_MR_STAR_TopPan(g),
                                  w01Last_MR_STAR_TopPan(g);
            a01Last_MC_STAR_Top.View( a01Last_MC_STAR, 0, 0, a01.Height(), 1 ); 
            w01Last_MC_STAR_Top.View( w01Last_MC_STAR, 0, 0, a01.Height(), 1 );
            a01Last_MR_STAR_TopPan.View
            ( a01Last_MR_STAR, topSize, 0, a01.Height()-topSize, 1 );
            w01Last_MR_STAR_TopPan.View
            ( w01Last_MR_STAR, topSize, 0, a01.Height()-topSize, 1 );
            const C* a01_MC_STAR_Buffer = a01Last_MC_STAR_Top.LocalBuffer();
            const C* w01_MC_STAR_Buffer = w01Last_MC_STAR_Top.LocalBuffer();
            const C* a01_MR_STAR_Buffer = a01Last_MR_STAR_TopPan.LocalBuffer();
            const C* w01_MR_STAR_Buffer = w01Last_MR_STAR_TopPan.LocalBuffer();
            C* A00PanBuffer = A00Pan.LocalBuffer();
            const int localHeight = A00Pan.LocalHeight();
            const int localWidth = A00Pan.LocalWidth();
            const int lDim = A00Pan.LocalLDim();
            for( int jLocal=0; jLocal<localWidth; ++jLocal )
                for( int iLocal=0; iLocal<localHeight; ++iLocal )
                    A00PanBuffer[iLocal+jLocal*lDim] -=
                        w01_MC_STAR_Buffer[iLocal]*
                        Conj(a01_MR_STAR_Buffer[jLocal]) +
                        a01_MC_STAR_Buffer[iLocal]*
                        Conj(w01_MR_STAR_Buffer[jLocal]);

            // We are through with the last iteration's w01
            w01Last_MC_STAR.FreeAlignments();
            w01Last_MR_STAR.FreeAlignments();
        }

        // Form the local portions of (A00 a01) into p01[MC,* ] and q01[MR,* ]:
        //   p01[MC,* ] := triu(A00)[MC,MR] a01[MR,* ]
        //   q01[MR,* ] := triu(A00,+1)'[MR,MC] a01[MC,* ]
        if( A00.ColShift() < A00.RowShift() )
        {
            // We are above the diagonal, so we can multiply without an 
            // offset for triu(A00)[MC,MR] and triu(A00,1)'[MR,MC]
            if( A00.LocalHeight() != 0 )
            {
                // Our local portion of p01[MC,* ] might be one entry longer
                // than A00.LocalWidth(), so go ahead and set the last entry
                // to 0.
                C* p01_MC_STAR_LocalBuffer = p01_MC_STAR.LocalBuffer();
                p01_MC_STAR_LocalBuffer[A00.LocalHeight()-1] = 0;
                MemCopy
                ( p01_MC_STAR.LocalBuffer(),
                  a01_MR_STAR.LocalBuffer(),
                  A00.LocalWidth() );
                blas::Trmv
                ( 'U', 'N', 'N', A00.LocalWidth(),
                  A00.LocalBuffer(), A00.LocalLDim(),
                  p01_MC_STAR.LocalBuffer(), 1 );
            }
            if( A00.LocalWidth() != 0 )
            {
                MemCopy
                ( q01_MR_STAR.LocalBuffer(),
                  a01_MC_STAR.LocalBuffer(), A00.LocalWidth() );
                blas::Trmv
                ( 'U', 'C', 'N', A00.LocalWidth(),
                  A00.LocalBuffer(), A00.LocalLDim(),
                  q01_MR_STAR.LocalBuffer(), 1 );
            }
        }
        else if( A00.ColShift() > A00.RowShift() )
        {
            // We are below the diagonal, so we need to use an offset 
            // for both triu(A00)[MC,MR] and triu(A00,+1)'[MR,MC]
            const C* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
            const C* a01_MR_STAR_LocalBuffer = a01_MR_STAR.LocalBuffer();
            const C* A00LocalBuffer = A00.LocalBuffer();
            if( A00.LocalHeight() != 0 )
            {
                // The last entry of p01[MC,* ] will be zero due to the forced
                // offset
                C* p01_MC_STAR_LocalBuffer = p01_MC_STAR.LocalBuffer();
                p01_MC_STAR_LocalBuffer[A00.LocalHeight()-1] = 0;
                MemCopy
                ( &p01_MC_STAR_LocalBuffer[0],
                  &a01_MR_STAR_LocalBuffer[1], A00.LocalWidth()-1 );
                blas::Trmv
                ( 'U', 'N', 'N', A00.LocalWidth()-1,
                  &A00LocalBuffer[A00.LocalLDim()], A00.LocalLDim(),
                  &p01_MC_STAR_LocalBuffer[0], 1 );
            }
            if( A00.LocalWidth() != 0 )
            {
                // The first entry of q01[MR,* ] will be zero due to the forced
                // offset
                C* q01_MR_STAR_LocalBuffer = q01_MR_STAR.LocalBuffer();
                q01_MR_STAR_LocalBuffer[0] = 0;
                MemCopy
                ( &q01_MR_STAR_LocalBuffer[1],
                  &a01_MC_STAR_LocalBuffer[0], A00.LocalWidth()-1 );
                blas::Trmv
                ( 'U', 'C', 'N', A00.LocalWidth()-1,
                  &A00LocalBuffer[A00.LocalLDim()], A00.LocalLDim(),
                  &q01_MR_STAR_LocalBuffer[1], 1 );
            }
        }
        else
        {
            // We are on the diagonal, so we only need an offset for
            // triu(A00,+1)'[MR,MC]
            if( A00.LocalWidth() != 0 )
            {
                MemCopy
                ( p01_MC_STAR.LocalBuffer(),
                  a01_MR_STAR.LocalBuffer(), A00.LocalHeight() );
                blas::Trmv
                ( 'U', 'N', 'N', A00.LocalHeight(),
                  A00.LocalBuffer(), A00.LocalLDim(),
                  p01_MC_STAR.LocalBuffer(), 1 );

                // The first entry of q01[MR,* ] must be zero due to the forced
                // offset
                const C* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
                const C* A00LocalBuffer = A00.LocalBuffer();
                C* q01_MR_STAR_LocalBuffer = q01_MR_STAR.LocalBuffer();
                q01_MR_STAR_LocalBuffer[0] = 0;
                MemCopy
                ( &q01_MR_STAR_LocalBuffer[1],
                  &a01_MC_STAR_LocalBuffer[0], A00.LocalWidth()-1 );
                blas::Trmv
                ( 'U', 'C', 'N', A00.LocalWidth()-1,
                  &A00LocalBuffer[A00.LocalLDim()], A00.LocalLDim(),
                  &q01_MR_STAR_LocalBuffer[1], 1 );
            }
        }

        LocalGemv( ADJOINT, (C)1, W02T, a01T_MC_STAR, (C)0, x21_MR_STAR );
        LocalGemv( ADJOINT, (C)1, A02T, a01T_MC_STAR, (C)0, y21_MR_STAR );

        // Combine the AllReduce column summations of x21[MR,* ] and y21[MR,* ]
        {
            const int x21LocalHeight = x21_MR_STAR.LocalHeight();
            const int y21LocalHeight = y21_MR_STAR.LocalHeight();
            const int reduceSize = x21LocalHeight+y21LocalHeight;
            std::vector<C> colSumSendBuffer(reduceSize);
            std::vector<C> colSumRecvBuffer(reduceSize);
            MemCopy
            ( &colSumSendBuffer[0], 
              x21_MR_STAR.LocalBuffer(), x21LocalHeight );
            MemCopy
            ( &colSumSendBuffer[x21LocalHeight],
              y21_MR_STAR.LocalBuffer(), y21LocalHeight );
            mpi::AllReduce
            ( &colSumSendBuffer[0], 
              &colSumRecvBuffer[0],
              reduceSize, mpi::SUM, g.ColComm() );
            MemCopy
            ( x21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[0], x21LocalHeight );
            MemCopy
            ( y21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[x21LocalHeight], y21LocalHeight );
        }

        LocalGemv( NORMAL, (C)-1, A02T, x21_MR_STAR, (C)1, p01T_MC_STAR );
        LocalGemv( NORMAL, (C)-1, W02T, y21_MR_STAR, (C)1, p01T_MC_STAR );

        // Fast transpose the unsummed q01[MR,* ] -> q01[MC,* ], so that
        // it needs to be summed over process rows instead of process 
        // columns. We immediately add it onto p01[MC,* ], which also needs
        // to be summed within process rows.
        if( onDiagonal )
        {
            const int a01LocalHeight = a01.LocalHeight();
            C* p01_MC_STAR_LocalBuffer = p01_MC_STAR.LocalBuffer();
            const C* q01_MR_STAR_LocalBuffer = q01_MR_STAR.LocalBuffer();
            for( int i=0; i<a01LocalHeight; ++i )
                p01_MC_STAR_LocalBuffer[i] += q01_MR_STAR_LocalBuffer[i];
        }
        else
        {
            // Pairwise exchange with the transpose process
            const int sendSize = A00.LocalWidth();
            const int recvSize = A00.LocalHeight();
            std::vector<C> recvBuffer(recvSize);
            mpi::SendRecv
            ( q01_MR_STAR.LocalBuffer(), sendSize, transposeRank, 0,
              &recvBuffer[0],            recvSize, transposeRank, 0,
              g.VCComm() );

            // Unpack the recv buffer directly onto p01[MC,* ]
            C* p01_MC_STAR_LocalBuffer = p01_MC_STAR.LocalBuffer();
            for( int i=0; i<recvSize; ++i )
                p01_MC_STAR_LocalBuffer[i] += recvBuffer[i];
        }

        if( W00.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // Reduce to one p01[MC,* ] to the next process column.
            const int a01LocalHeight = a01.LocalHeight();

            const int nextProcessRow = (alpha11.ColAlignment()+r-1) % r;
            const int nextProcessCol = (alpha11.RowAlignment()+r-1) % r;

            std::vector<C> reduceToOneRecvBuffer(a01LocalHeight);
            mpi::Reduce
            ( p01_MC_STAR.LocalBuffer(),
              &reduceToOneRecvBuffer[0],
              a01LocalHeight, mpi::SUM, nextProcessCol, g.RowComm() );
            if( g.Col() == nextProcessCol )
            {
                // Finish computing w01. During its computation, ensure that 
                // every process has a copy of the last element of the w01.
                // We know a priori that the last element of a01 is one.
                const C* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
                C myDotProduct = blas::Dot
                    ( a01LocalHeight, &reduceToOneRecvBuffer[0],   1, 
                                      &a01_MC_STAR_LocalBuffer[0], 1 );
                C sendBuffer[2];
                C recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.Row()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[a01LocalHeight-1] : 0 );
                mpi::AllReduce
                ( sendBuffer, recvBuffer, 2, mpi::SUM, g.ColComm() );
                C dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w01LastLocalBuffer
                // - w01LastBottomEntry
                C scale = static_cast<C>(0.5)*dotProduct*Conj(tau);
                for( int i=0; i<a01LocalHeight; ++i )
                    w01LastLocalBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a01_MC_STAR_LocalBuffer[i] );
                w01LastBottomEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w01[MC,* ] and w01[MR,* ]
            const int a01LocalHeight = a01.LocalHeight();

            // AllReduce sum p01[MC,* ] over process rows
            std::vector<C> allReduceRecvBuffer(a01LocalHeight);
            mpi::AllReduce
            ( p01_MC_STAR.LocalBuffer(),
              &allReduceRecvBuffer[0],
              a01LocalHeight, mpi::SUM, g.RowComm() );

            // Finish computing w01. During its computation, ensure that 
            // every process has a copy of the last element of the w01.
            // We know a priori that the last element of a01 is one.
            const C* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
            C myDotProduct = blas::Dot
                ( a01LocalHeight, &allReduceRecvBuffer[0],     1, 
                                  &a01_MC_STAR_LocalBuffer[0], 1 );
            C dotProduct;
            mpi::AllReduce
            ( &myDotProduct, &dotProduct, 1, mpi::SUM, g.ColComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<C,MC,STAR> w01_MC_STAR(g);
            DistMatrix<C,MR,STAR> w01_MR_STAR(g);
            w01_MC_STAR.View( W_MC_STAR, 0, W00.Width(), w01.Height(), 1 );
            w01_MR_STAR.View( W_MR_STAR, 0, W00.Width(), w01.Height(), 1 );

            // Store w01[MC,* ]
            C scale = static_cast<C>(0.5)*dotProduct*Conj(tau);
            C* w01_MC_STAR_LocalBuffer = w01_MC_STAR.LocalBuffer();
            for( int i=0; i<a01LocalHeight; ++i )
                w01_MC_STAR_LocalBuffer[i] = tau*
                    ( allReduceRecvBuffer[i]-
                      scale*a01_MC_STAR_LocalBuffer[i] );

            // Fast transpose w01[MC,* ] -> w01[MR,* ]
            if( onDiagonal )
            {
                MemCopy
                ( w01_MR_STAR.LocalBuffer(),
                  w01_MC_STAR.LocalBuffer(), a01LocalHeight );
            }
            else
            {
                // Pairwise exchange with the transpose process
                const int sendSize = A00.LocalHeight();
                const int recvSize = A00.LocalWidth();
                mpi::SendRecv
                ( w01_MC_STAR.LocalBuffer(), sendSize, transposeRank, 0,
                  w01_MR_STAR.LocalBuffer(), recvSize, transposeRank, 0,
                  g.VCComm() );
            }
        }
        //--------------------------------------------------------------------//
        a01_MC_STAR.FreeAlignments();
        a01_MR_STAR.FreeAlignments();
        p01_MC_STAR.FreeAlignments();
        p01_MR_STAR.FreeAlignments();
        q01_MR_STAR.FreeAlignments();
        x21_MR_STAR.FreeAlignments();
        y21_MR_STAR.FreeAlignments();

        SlidePartitionUp
        ( tT,  t0,
         /**/ /****/
               tau1,
          tB,  t2 );

        SlidePartitionUp
        ( eT,  e0,
         /**/ /********/
               epsilon1,
          eB,  e2 );

        SlidePartitionUpDiagonal
        ( WTL, /**/ WTR,  W00, /**/ w01,     W02,
         /*************/ /**********************/
               /**/       w10, /**/ omega11, w12,
          WBL, /**/ WBR,  W20, /**/ w21,     W22 );
 
        SlidePartitionUpDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        firstIteration = false;
    }
    PopBlocksizeStack();

    expandedABR.SetRealPartOfDiagonal( e, 1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
