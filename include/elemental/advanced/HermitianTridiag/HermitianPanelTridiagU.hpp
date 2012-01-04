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

namespace elemental {

template<typename R>
inline void
internal::HermitianPanelTridiagU
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR, 
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR )
{
    const int panelSize = W.Width();
    const int topSize = W.Height()-panelSize;
#ifndef RELEASE
    PushCallStack("internal::HermitianPanelTridiagU");
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
    const int c = g.Width();
    const int p = g.Size();

    // Create a distributed matrix for storing the superdiagonal
    DistMatrix<R,MD,STAR> e(g);
    DistMatrix<R,MC,MR> expandedABR(g);
    expandedABR.View( A, topSize-1, topSize-1, panelSize+1, panelSize+1 );
    e.AlignWithDiagonal( expandedABR, 1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g), a01T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),           alpha01B(g),
                         A20(g), a21(g),     A22(g),  A02T(g), A00Pan(g);
    DistMatrix<R,MC,MR> 
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
    DistMatrix<R,MC,STAR> q01_MC_STAR(g);
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
        q01_MC_STAR.AlignWith( A00 );
        q01_MR_STAR.AlignWith( A00 );
        x21_MR_STAR.AlignWith( A02T );
        y21_MR_STAR.AlignWith( A02T );

        a01_MC_STAR.ResizeTo( a01.Height(), 1 );
        a01_MR_STAR.ResizeTo( a01.Height(), 1 );
        p01_MC_STAR.ResizeTo( a01.Height(), 1 );
        p01_MR_STAR.ResizeTo( a01.Height(), 1 );
        q01_MC_STAR.ResizeTo( a01.Height(), 1 );
        q01_MR_STAR.ResizeTo( a01.Height(), 1 );
        x21_MR_STAR.ResizeTo( A02.Width(), 1 );
        y21_MR_STAR.ResizeTo( A02.Width(), 1 );

        // View the portions of a01[MC,* ] and p01[MC,* ] above the current
        // panel's square
        a01T_MC_STAR.View( a01_MC_STAR, 0, 0, topSize, 1 );
        p01T_MC_STAR.View( p01_MC_STAR, 0, 0, topSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.MRRank() == alpha11.RowAlignment() );
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
            tau = internal::ColReflector( alpha01B, a01T );
        }
        // Store the subdiagonal value and turn a01 into a proper reflector
        // by explicitly placing the implicit one in its bottom entry
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
                std::memcpy
                ( &rowBroadcastBuffer[0], 
                  a01.LocalBuffer(), 
                  a01LocalHeight*sizeof(R) );
                rowBroadcastBuffer[a01LocalHeight] = tau;
            }
            // Broadcast a01 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a01LocalHeight+1, a01.RowAlignment(), g.MRComm() );
            // Store a01[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            std::memcpy
            ( a01_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0],
              a01LocalHeight*sizeof(R) );
            // Store a01[MC,* ] into APan[MC,* ]
            std::memcpy
            ( APan_MC_STAR.LocalBuffer(0,W00.Width()), 
              &rowBroadcastBuffer[0],
              a01LocalHeight*sizeof(R) );
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight];
            
            a01_MR_STAR = a01_MC_STAR;
            // Store a01[MR,* ]
            std::memcpy
            ( APan_MR_STAR.LocalBuffer(0,W00.Width()),
              a01_MR_STAR.LocalBuffer(),
              a01_MR_STAR.LocalHeight()*sizeof(R) );
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
                std::memcpy
                ( &rowBroadcastBuffer[0], 
                  a01.LocalBuffer(),
                  a01LocalHeight*sizeof(R) );
                std::memcpy
                ( &rowBroadcastBuffer[a01LocalHeight], 
                  &w01LastLocalBuffer[0],
                  w01LastLocalHeight*sizeof(R) );
                rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight] = tau;
            }
            // Broadcast a01, w01Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a01LocalHeight+w01LastLocalHeight+1, 
              a01.RowAlignment(), g.MRComm() );
            // Store a01[MC,* ] into its DistMatrix class 
            std::memcpy
            ( a01_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0],
              a01LocalHeight*sizeof(R) );
            // Store a01[MC,* ] into APan[MC,* ]
            std::memcpy
            ( APan_MC_STAR.LocalBuffer(0,W00.Width()),
              &rowBroadcastBuffer[0],
              a01LocalHeight*sizeof(R) );
            // Store w01Last[MC,* ] into its DistMatrix class
            w01Last_MC_STAR.AlignWith( A00 );
            w01Last_MC_STAR.ResizeTo( a01.Height()+1, 1 );
            std::memcpy
            ( w01Last_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[a01LocalHeight], 
              w01LastLocalHeight*sizeof(R) );
            // Store the bottom part of w01Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w01.
            std::memcpy
            ( W_MC_STAR.LocalBuffer(0,W00.Width()+1),
              &rowBroadcastBuffer[a01LocalHeight],
              w01LastLocalHeight*sizeof(R) );
            if( g.MRRank() == w01Last.RowAlignment() )
            {
                std::memcpy
                ( w01Last.LocalBuffer(),
                  &rowBroadcastBuffer[a01LocalHeight],
                  w01LastLocalHeight*sizeof(R) );
            }
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight];

            // Form a01[MR,* ] and w01Last[MR,* ] by combining the 
            // communications needed for taking a vector from 
            // [MC,* ] -> [MR,* ]: 
            //   local copy to [VC,* ], 
            //   Send/Recv to [VR,* ], 
            //   AllGather to [MR,* ]
            // We can combine the two by treating a01 as [a01; 0].

            const int colAlignSource = A00.ColAlignment();
            const int colAlignDest = A00.RowAlignment();
            const int colShiftSource = A00.ColShift();
            const int colShiftDest = A00.RowShift();

            const int height = a01.Height()+1;
            const int portionSize = 
                std::max(2*MaxLocalLength(height,p),mpi::MIN_COLL_MSG);

            const int colShiftVRDest = Shift(g.VRRank(),colAlignDest,p);
            const int colShiftVCSource = Shift(g.VCRank(),colAlignSource,p);
            const int sendRankRM = 
                (g.VRRank()+(p+colShiftVCSource-colShiftVRDest))%p;
            const int recvRankCM = 
                (g.VCRank()+(p+colShiftVRDest-colShiftVCSource))%p;
            const int recvRankRM = 
                (recvRankCM/r)+c*(recvRankCM%r);

            std::vector<R> transposeBuffer( (r+1)*portionSize );
            R* sendBuf = &transposeBuffer[0];
            R* recvBuf = &transposeBuffer[r*portionSize];

            // (w01Last[VC,* ] <- w01Last[MC,* ]) and
            // ([a01; 0][VC,* ] <- [a01; 0][MC,* ])
            {
                // Pack the necessary portion of w01Last[MC,* ]
                const int shift = Shift(g.VCRank(),colAlignSource,p);
                const int offset = (shift-colShiftSource)/r;
                const int w01LocalHeight = LocalLength(height,shift,p);
                const R* w01LastBuffer = w01Last_MC_STAR.LocalBuffer(offset,0);
                for( int i=0; i<w01LocalHeight; ++i )
                    sendBuf[i] = w01LastBuffer[i*c];
                
                // Pack the necessary portion of a01[MC,* ]
                const int a01LocalHeight = LocalLength(height-1,shift,p);
                const R* a01Buffer = a01_MC_STAR.LocalBuffer(offset,0);
                for( int i=0; i<a01LocalHeight; ++i )
                    sendBuf[w01LocalHeight+i] = a01Buffer[i*c];
            }

            // [VR,* ] <- [VC,* ]
            mpi::SendRecv
            ( sendBuf, portionSize, sendRankRM, 0,
              recvBuf, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

            // [MR,* ] <- [VR,* ]
            mpi::AllGather
            ( recvBuf, portionSize,
              sendBuf, portionSize, g.MCComm() );

            // Unpack
            w01Last_MR_STAR.AlignWith( A00 );
            w01Last_MR_STAR.ResizeTo( a01.Height()+1, 1 );
            for( int k=0; k<r; ++k )
            {
                // Unpack into w01Last[MR,* ]
                const R* w01Data = &sendBuf[k*portionSize];
                const int shift = Shift(g.MRRank()+c*k,colAlignDest,p);
                const int offset = (shift-colShiftDest) / c;
                const int w01LocalHeight = LocalLength(height,shift,p);
                R* w01LastBuffer = w01Last_MR_STAR.LocalBuffer(offset,0);
                for( int i=0; i<w01LocalHeight; ++i )
                    w01LastBuffer[i*r] = w01Data[i];

                // Unpack into a01[MR,* ]
                const R* a01Data = &sendBuf[k*portionSize+w01LocalHeight];
                const int a01LocalHeight = LocalLength(height-1,shift,p);
                R* a01Buffer = a01_MR_STAR.LocalBuffer(offset,0);
                for( int i=0; i<a01LocalHeight; ++i )
                    a01Buffer[i*r] = a01Data[i];
            }
            // Store w01Last[MR,* ]
            std::memcpy
            ( W_MR_STAR.LocalBuffer(0,W00.Width()+1),
              w01Last_MR_STAR.LocalBuffer(),
              w01Last_MR_STAR.LocalHeight()*sizeof(R) );
            // Store a01[MR,* ]
            std::memcpy
            ( APan_MR_STAR.LocalBuffer(0,W00.Width()),
              a01_MR_STAR.LocalBuffer(),
              a01_MR_STAR.LocalHeight()*sizeof(R) );

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
        PopBlocksizeStack();
        p01_MC_STAR.SetToZero();
        q01_MR_STAR.SetToZero();
        internal::LocalSymvColAccumulateU
        ( (R)1, A00, a01_MC_STAR, a01_MR_STAR, p01_MC_STAR, q01_MR_STAR );
        PushBlocksizeStack( 1 );

        Gemv
        ( TRANSPOSE, 
          (R)1, W02T.LockedLocalMatrix(),
                a01T_MC_STAR.LockedLocalMatrix(),
          (R)0, x21_MR_STAR.LocalMatrix() );
        Gemv
        ( TRANSPOSE, 
          (R)1, A02T.LockedLocalMatrix(),
                a01T_MC_STAR.LockedLocalMatrix(),
          (R)0, y21_MR_STAR.LocalMatrix() );
        // Combine the AllReduce column summations of x21[MR,* ], y21[MR,* ],
        // and q01[MR,* ]
        {
            const int x21LocalHeight = x21_MR_STAR.LocalHeight();
            const int y21LocalHeight = y21_MR_STAR.LocalHeight();
            const int q01LocalHeight = q01_MR_STAR.LocalHeight();
            const int reduceSize = x21LocalHeight+y21LocalHeight+q01LocalHeight;
            std::vector<R> colSumSendBuffer(reduceSize);
            std::vector<R> colSumRecvBuffer(reduceSize);
            std::memcpy
            ( &colSumSendBuffer[0], 
              x21_MR_STAR.LocalBuffer(), 
              x21LocalHeight*sizeof(R) );
            std::memcpy
            ( &colSumSendBuffer[x21LocalHeight],
              y21_MR_STAR.LocalBuffer(), 
              y21LocalHeight*sizeof(R) );
            std::memcpy
            ( &colSumSendBuffer[x21LocalHeight+y21LocalHeight],
              q01_MR_STAR.LocalBuffer(), 
              q01LocalHeight*sizeof(R) );
            mpi::AllReduce
            ( &colSumSendBuffer[0], 
              &colSumRecvBuffer[0],
              reduceSize, mpi::SUM, g.MCComm() );
            std::memcpy
            ( x21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[0], 
              x21LocalHeight*sizeof(R) );
            std::memcpy
            ( y21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[x21LocalHeight], 
              y21LocalHeight*sizeof(R) );
            std::memcpy
            ( q01_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[x21LocalHeight+y21LocalHeight], 
              q01LocalHeight*sizeof(R) );
        }

        Gemv
        ( NORMAL, 
          (R)-1, A02T.LockedLocalMatrix(),
                 x21_MR_STAR.LockedLocalMatrix(),
          (R)+1, p01T_MC_STAR.LocalMatrix() );
        Gemv
        ( NORMAL, 
          (R)-1, W02T.LockedLocalMatrix(),
                 y21_MR_STAR.LockedLocalMatrix(),
          (R)+1, p01T_MC_STAR.LocalMatrix() );

        if( W00.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // combine the Reduce to one of p01[MC,* ] with the redistribution 
            // of q01[MR,* ] -> q01[MC,MR] to the next process column.
            const int localHeight = p01_MC_STAR.LocalHeight();
            std::vector<R> reduceToOneSendBuffer(2*localHeight);
            std::vector<R> reduceToOneRecvBuffer(2*localHeight);

            // Pack p01[MC,* ]
            std::memcpy
            ( &reduceToOneSendBuffer[0], 
              p01_MC_STAR.LocalBuffer(),
              localHeight*sizeof(R) );

            // Fill in contributions to q01[MC,MR] from q01[MR,* ]
            const bool contributing = 
                ( q01_MR_STAR.ColShift() % g.GCD() ==
                  p01_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    std::memcpy
                    ( &reduceToOneSendBuffer[localHeight],
                      q01_MR_STAR.LocalBuffer(), 
                      localHeight*sizeof(R) );
                }
                else
                {
                    // Zero the entire buffer first
                    std::memset
                    ( &reduceToOneSendBuffer[localHeight], 0,
                      localHeight*sizeof(R) );
                    // Fill in the entries that we contribute to.
                    // We seek to find the minimum s in N such that
                    //   s*c = a0-b0 (mod r)
                    // where a0 is the column shift of MC, b0 is the row shift
                    // of MR, and s is our first local entry of MR that will 
                    // contribute to MC. I cannot think of an O(1) method, so
                    // I will instead use the worst-case O(lcm(c,r)/c) method.
                    const int sourcePeriod = g.LCM() / c;
                    const int targetPeriod = g.LCM() / r;
                    const int a0 = p01_MC_STAR.ColShift();
                    const int b0 = q01_MR_STAR.ColShift();

                    int sourceStart = 0;
                    const int f = (r+a0-b0) % r;
                    for( int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    const int globalShift = b0+sourceStart*c;
                    const int targetStart = (globalShift-a0)/r;
                    const int localLength =
                        LocalLength(localHeight,targetStart,targetPeriod);
                    const R* q01_MR_STAR_LocalBuffer = 
                        q01_MR_STAR.LocalBuffer();
                    const int offset = localHeight + targetStart;
                    for( int i=0; i<localLength; ++i )                        
                        reduceToOneSendBuffer[offset+i*targetPeriod] = 
                            q01_MR_STAR_LocalBuffer[sourceStart+i*sourcePeriod];
                }
            }
            else
            {
                std::memset
                ( &reduceToOneSendBuffer[localHeight], 0, 
                  localHeight*sizeof(R) );
            }

            const int nextProcessRow = (alpha11.ColAlignment()+r-1) % r;
            const int nextProcessCol = (alpha11.RowAlignment()+c-1) % c;
            mpi::Reduce
            ( &reduceToOneSendBuffer[0], 
              &reduceToOneRecvBuffer[0],
              2*localHeight, mpi::SUM, nextProcessCol, g.MRComm() );
            if( g.MRRank() == nextProcessCol )
            {
                // Combine the second half into the first half        
                for( int i=0; i<localHeight; ++i )
                    reduceToOneRecvBuffer[i] +=
                        reduceToOneRecvBuffer[i+localHeight];

                // Finish computing w01. During its computation, ensure that 
                // every process has a copy of the bottom element of the w01.
                // We know a priori that the bottom element of a01 is one.
                const R* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
                R myDotProduct = blas::Dot
                    ( localHeight, &reduceToOneRecvBuffer[0], 1, 
                                   a01_MC_STAR_LocalBuffer,   1 );
                R sendBuffer[2];
                R recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.MCRank()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[localHeight-1] : 0 );
                mpi::AllReduce
                ( sendBuffer, recvBuffer, 2, mpi::SUM, g.MCComm() );
                R dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w01LastLocalBuffer
                // - w01LastBottomEntry
                R scale = 0.5*dotProduct*tau;
                for( int i=0; i<localHeight; ++i )
                    w01LastLocalBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a01_MC_STAR_LocalBuffer[i] );
                w01LastBottomEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w01[MC,* ] and w01[MR,* ] so that we may place them into W[MC,* ]
            // and W[MR,* ]
            const int localHeight = p01_MC_STAR.LocalHeight();
            std::vector<R> allReduceSendBuffer(2*localHeight),
                           allReduceRecvBuffer(2*localHeight);

            // Pack p01[MC,* ]
            std::memcpy
            ( &allReduceSendBuffer[0], 
              p01_MC_STAR.LocalBuffer(),
              localHeight*sizeof(R) );

            // Fill in contributions to q01[MC,* ] from q01[MR,* ]
            const bool contributing = 
                ( q01_MR_STAR.ColShift() % g.GCD() ==
                  p01_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    std::memcpy
                    ( &allReduceSendBuffer[localHeight],
                      q01_MR_STAR.LocalBuffer(), 
                      localHeight*sizeof(R) );
                }
                else
                {
                    // Zero the entire buffer first
                    std::memset
                    ( &allReduceSendBuffer[localHeight], 0, 
                      localHeight*sizeof(R) );
                    // Fill in the entries that we contribute to.
                    // We seek to find the minimum s in N such that
                    //   s*c = a0-b0 (mod r)
                    // where a0 is the column shift of MC, b0 is the row shift
                    // of MR, and s is our first local entry of MR that will 
                    // contribute to MC. I cannot think of an O(1) method, so
                    // I will instead use the worst-case O(lcm(c,r)/c) method.
                    const int sourcePeriod = g.LCM() / c;
                    const int targetPeriod = g.LCM() / r;
                    const int a0 = p01_MC_STAR.ColShift();
                    const int b0 = q01_MR_STAR.ColShift();

                    int sourceStart = 0;
                    const int f = (r+a0-b0) % r;
                    for( int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    const int globalShift = b0+sourceStart*c;
                    const int targetStart = (globalShift-a0)/r;
                    const int localLength = 
                        LocalLength(localHeight,targetStart,targetPeriod);
                    const R* q01_MR_STAR_LocalBuffer = 
                        q01_MR_STAR.LocalBuffer();
                    const int offset = localHeight + targetStart;
                    for( int i=0; i<localLength; ++i )
                        allReduceSendBuffer[offset+i*targetPeriod] = 
                            q01_MR_STAR_LocalBuffer[sourceStart+i*sourcePeriod];
                }
            }
            else
            {
                std::memset
                ( &allReduceSendBuffer[localHeight], 0, 
                  localHeight*sizeof(R) );
            }

            mpi::AllReduce
            ( &allReduceSendBuffer[0], 
              &allReduceRecvBuffer[0],
              2*localHeight, mpi::SUM, g.MRComm() );

            // Combine the second half into the first half        
            for( int i=0; i<localHeight; ++i )
                allReduceRecvBuffer[i] += allReduceRecvBuffer[i+localHeight];
 
            // Finish computing w01. 
            const R* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
            R myDotProduct = blas::Dot
                ( localHeight, &allReduceRecvBuffer[0], 1, 
                               a01_MC_STAR_LocalBuffer, 1 );
            R dotProduct;
            mpi::AllReduce
            ( &myDotProduct, &dotProduct, 1, mpi::SUM, g.MCComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<R,MC,STAR> w01_MC_STAR(g);
            DistMatrix<R,MR,STAR> w01_MR_STAR(g);
            w01_MC_STAR.View( W_MC_STAR, 0, W00.Width(), w01.Height(), 1 );
            w01_MR_STAR.View( W_MR_STAR, 0, W00.Width(), w01.Height(), 1 );

            // Store w01[MC,* ]
            R scale = 0.5*dotProduct*tau;
            R* w01_MC_STAR_LocalBuffer = w01_MC_STAR.LocalBuffer();
            for( int i=0; i<localHeight; ++i )
                w01_MC_STAR_LocalBuffer[i] = tau*
                    ( allReduceRecvBuffer[i]-
                      scale*a01_MC_STAR_LocalBuffer[i] );

            // Form w01[MR,* ]
            w01_MR_STAR = w01_MC_STAR;
        }
        //--------------------------------------------------------------------//
        a01_MC_STAR.FreeAlignments();
        a01_MR_STAR.FreeAlignments();
        p01_MC_STAR.FreeAlignments();
        p01_MR_STAR.FreeAlignments();
        q01_MC_STAR.FreeAlignments();
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
internal::HermitianPanelTridiagU
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MC,MR  >& W,
  DistMatrix<std::complex<R>,MD,STAR>& t,
  DistMatrix<std::complex<R>,MC,STAR>& APan_MC_STAR, 
  DistMatrix<std::complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<std::complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<std::complex<R>,MR,STAR>& W_MR_STAR )
{
    const int panelSize = W.Width();
    const int topSize = W.Height()-panelSize;
#ifndef RELEASE
    PushCallStack("internal::HermitianPanelTridiagU");
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
    typedef std::complex<R> C;

    const Grid& g = A.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();

    // Create a distributed matrix for storing the superdiagonal
    DistMatrix<R,MD,STAR> e(g);
    DistMatrix<C,MC,MR> expandedABR(g);
    expandedABR.View( A, topSize-1, topSize-1, panelSize+1, panelSize+1 );
    e.AlignWithDiagonal( expandedABR, 1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g), a01T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),           alpha01B(g),
                         A20(g), a21(g),     A22(g),  A02T(g), A00Pan(g);
    DistMatrix<C,MC,MR> 
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
    DistMatrix<C,MC,STAR> a01_MC_STAR(g);
    DistMatrix<C,MC,STAR> a01T_MC_STAR(g);
    DistMatrix<C,MR,STAR> a01_MR_STAR(g);
    DistMatrix<C,MC,STAR> p01_MC_STAR(g);
    DistMatrix<C,MC,STAR> p01T_MC_STAR(g);
    DistMatrix<C,MR,STAR> p01_MR_STAR(g);
    DistMatrix<C,MC,STAR> q01_MC_STAR(g);
    DistMatrix<C,MR,STAR> q01_MR_STAR(g);
    DistMatrix<C,MR,STAR> x21_MR_STAR(g);
    DistMatrix<C,MR,STAR> y21_MR_STAR(g);
    DistMatrix<C,MC,STAR> a01Last_MC_STAR(g);
    DistMatrix<C,MR,STAR> a01Last_MR_STAR(g);
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
        q01_MC_STAR.AlignWith( A00 );
        q01_MR_STAR.AlignWith( A00 );
        x21_MR_STAR.AlignWith( A02T );
        y21_MR_STAR.AlignWith( A02T );
        
        a01_MC_STAR.ResizeTo( a01.Height(), 1 );
        a01_MR_STAR.ResizeTo( a01.Height(), 1 );
        p01_MC_STAR.ResizeTo( a01.Height(), 1 );
        p01_MR_STAR.ResizeTo( a01.Height(), 1 );
        q01_MC_STAR.ResizeTo( a01.Height(), 1 );
        q01_MR_STAR.ResizeTo( a01.Height(), 1 );
        x21_MR_STAR.ResizeTo( A02.Width(), 1 );
        y21_MR_STAR.ResizeTo( A02.Width(), 1 );

        // View the portions of a01[MC,* ] and p01[MC,* ] above the current
        // panel's square
        a01T_MC_STAR.View( a01_MC_STAR, 0, 0, topSize, 1 );
        p01T_MC_STAR.View( p01_MC_STAR, 0, 0, topSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.MRRank() == alpha11.RowAlignment() );
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
            tau = internal::ColReflector( alpha01B, a01T );
            if( g.MCRank() == alpha01B.ColAlignment() )
                tau1.SetLocalEntry(0,0,tau);
        }
        // Store the subdiagonal value and turn a01 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        alpha01B.GetRealDiagonal( epsilon1 );
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
                std::memcpy
                ( &rowBroadcastBuffer[0], 
                  a01.LocalBuffer(), 
                  a01LocalHeight*sizeof(C) );
                rowBroadcastBuffer[a01LocalHeight] = tau;
            }
            // Broadcast a01 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a01LocalHeight+1, a01.RowAlignment(), g.MRComm() );
            // Store a01[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            std::memcpy
            ( a01_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0],
              a01LocalHeight*sizeof(C) );
            // Store a01[MC,* ] into APan[MC,* ]
            std::memcpy
            ( APan_MC_STAR.LocalBuffer(0,W00.Width()), 
              &rowBroadcastBuffer[0],
              a01LocalHeight*sizeof(C) );
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight];
            
            a01_MR_STAR = a01_MC_STAR;
            // Store a01[MR,* ]
            std::memcpy
            ( APan_MR_STAR.LocalBuffer(0,W00.Width()),
              a01_MR_STAR.LocalBuffer(),
              a01_MR_STAR.LocalHeight()*sizeof(C) );
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
                std::memcpy
                ( &rowBroadcastBuffer[0], 
                  a01.LocalBuffer(),
                  a01LocalHeight*sizeof(C) );
                std::memcpy
                ( &rowBroadcastBuffer[a01LocalHeight], 
                  &w01LastLocalBuffer[0],
                  w01LastLocalHeight*sizeof(C) );
                rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight] = tau;
            }
            // Broadcast a01, w01Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a01LocalHeight+w01LastLocalHeight+1, 
              a01.RowAlignment(), g.MRComm() );
            // Store a01[MC,* ] into its DistMatrix class 
            std::memcpy
            ( a01_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0],
              a01LocalHeight*sizeof(C) );
            // Store a01[MC,* ] into APan[MC,* ]
            std::memcpy
            ( APan_MC_STAR.LocalBuffer(0,W00.Width()), 
              &rowBroadcastBuffer[0],
              a01LocalHeight*sizeof(C) );
            // Store w01Last[MC,* ] into its DistMatrix class
            w01Last_MC_STAR.AlignWith( A00 );
            w01Last_MC_STAR.ResizeTo( a01.Height()+1, 1 );
            std::memcpy
            ( w01Last_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[a01LocalHeight], 
              w01LastLocalHeight*sizeof(C) );
            // Store the bottom part of w01Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w01.
            std::memcpy
            ( W_MC_STAR.LocalBuffer(0,W00.Width()+1),
              &rowBroadcastBuffer[a01LocalHeight],
              w01LastLocalHeight*sizeof(C) );
            if( g.MRRank() == w01Last.RowAlignment() )
            {
                std::memcpy
                ( w01Last.LocalBuffer(),
                  &rowBroadcastBuffer[a01LocalHeight],
                  w01LastLocalHeight*sizeof(C) );
            }
            // Store tau
            tau = rowBroadcastBuffer[a01LocalHeight+w01LastLocalHeight];

            // Form a01[MR,* ] and w01Last[MR,* ] by combining the 
            // communications needed for taking a vector from 
            // [MC,* ] -> [MR,* ]: 
            //   local copy to [VC,* ], 
            //   Send/Recv to [VR,* ], 
            //   AllGather to [MR,* ]
            // We can combine the two by treating a01 as [a01; 0] 

            const int colAlignSource = A00.ColAlignment();
            const int colAlignDest = A00.RowAlignment();
            const int colShiftSource = A00.ColShift();
            const int colShiftDest = A00.RowShift();

            const int height = a01.Height()+1;
            const int portionSize = 
                std::max(2*MaxLocalLength(height,p),mpi::MIN_COLL_MSG);

            const int colShiftVRDest = Shift(g.VRRank(),colAlignDest,p);
            const int colShiftVCSource = Shift(g.VCRank(),colAlignSource,p);
            const int sendRankRM = 
                (g.VRRank()+(p+colShiftVCSource-colShiftVRDest))%p;
            const int recvRankCM = 
                (g.VCRank()+(p+colShiftVRDest-colShiftVCSource))%p;
            const int recvRankRM = 
                (recvRankCM/r)+c*(recvRankCM%r);

            std::vector<C> transposeBuffer( (r+1)*portionSize );
            C* sendBuf = &transposeBuffer[0];
            C* recvBuf = &transposeBuffer[r*portionSize];

            // (w01Last[VC,* ] <- w01Last[MC,* ]) and
            // ([a01; 0][VC,* ] <- [a01; 0][MC,* ])
            {
                // Pack the necessary portion of w01Last[MC,* ]
                const int shift = Shift(g.VCRank(),colAlignSource,p);
                const int offset = (shift-colShiftSource)/r;
                const int w01LocalHeight = LocalLength(height,shift,p);
                const C* w01LastBuffer = w01Last_MC_STAR.LocalBuffer(offset,0);
                for( int i=0; i<w01LocalHeight; ++i )
                    sendBuf[i] = w01LastBuffer[i*c];
                
                // Pack the necessary portion of a01[MC,* ]
                const int a01LocalHeight = LocalLength(height-1,shift,p);
                const C* a01Buffer = a01_MC_STAR.LocalBuffer(offset,0);
                for( int i=0; i<a01LocalHeight; ++i )
                    sendBuf[w01LocalHeight+i] = a01Buffer[i*c];
            }

            // [VR,* ] <- [VC,* ]
            mpi::SendRecv
            ( sendBuf, portionSize, sendRankRM, 0,
              recvBuf, portionSize, recvRankRM, mpi::ANY_TAG, g.VRComm() );

            // [MR,* ] <- [VR,* ]
            mpi::AllGather
            ( recvBuf, portionSize,
              sendBuf, portionSize, g.MCComm() );

            // Unpack
            w01Last_MR_STAR.AlignWith( A00 );
            w01Last_MR_STAR.ResizeTo( a01.Height()+1, 1 );
            for( int k=0; k<r; ++k )
            {
                // Unpack into w01Last[MR,* ]
                const C* w01Data = &sendBuf[k*portionSize];
                const int shift = Shift(g.MRRank()+c*k,colAlignDest,p);
                const int offset = (shift-colShiftDest) / c;
                const int w01LocalHeight = LocalLength(height,shift,p);
                C* w01LastBuffer = w01Last_MR_STAR.LocalBuffer(offset,0);
                for( int i=0; i<w01LocalHeight; ++i )
                    w01LastBuffer[i*r] = w01Data[i];

                // Unpack into a01[MR,* ]
                const C* a01Data = &sendBuf[k*portionSize+w01LocalHeight];
                const int a01LocalHeight = LocalLength(height-1,shift,p);
                C* a01Buffer = a01_MR_STAR.LocalBuffer(offset,0);
                for( int i=0; i<a01LocalHeight; ++i )
                    a01Buffer[i*r] = a01Data[i];
            }
            // Store w01Last[MR,* ]
            std::memcpy
            ( W_MR_STAR.LocalBuffer(0,W00.Width()+1),
              w01Last_MR_STAR.LocalBuffer(),
              w01Last_MR_STAR.LocalHeight()*sizeof(C) );
            // Store a01[MR,* ]
            std::memcpy
            ( APan_MR_STAR.LocalBuffer(0,W00.Width()),
              a01_MR_STAR.LocalBuffer(),
              a01_MR_STAR.LocalHeight()*sizeof(C) );

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
        PopBlocksizeStack();
        p01_MC_STAR.SetToZero();
        q01_MR_STAR.SetToZero();
        internal::LocalHemvColAccumulateU
        ( (C)1, A00, a01_MC_STAR, a01_MR_STAR, p01_MC_STAR, q01_MR_STAR );
        PushBlocksizeStack( 1 );

        Gemv
        ( ADJOINT, 
          (C)1, W02T.LockedLocalMatrix(),
                a01T_MC_STAR.LockedLocalMatrix(),
          (C)0, x21_MR_STAR.LocalMatrix() );
        Gemv
        ( ADJOINT, 
          (C)1, A02T.LockedLocalMatrix(),
                a01T_MC_STAR.LockedLocalMatrix(),
          (C)0, y21_MR_STAR.LocalMatrix() );
        // Combine the AllReduce column summations of x21[MR,* ], y21[MR,* ],
        // and q01[MR,* ]
        {
            const int x21LocalHeight = x21_MR_STAR.LocalHeight();
            const int y21LocalHeight = y21_MR_STAR.LocalHeight();
            const int q01LocalHeight = q01_MR_STAR.LocalHeight();
            const int reduceSize = x21LocalHeight+y21LocalHeight+q01LocalHeight;
            std::vector<C> colSumSendBuffer(reduceSize);
            std::vector<C> colSumRecvBuffer(reduceSize);
            std::memcpy
            ( &colSumSendBuffer[0], 
              x21_MR_STAR.LocalBuffer(), 
              x21LocalHeight*sizeof(C) );
            std::memcpy
            ( &colSumSendBuffer[x21LocalHeight],
              y21_MR_STAR.LocalBuffer(), 
              y21LocalHeight*sizeof(C) );
            std::memcpy
            ( &colSumSendBuffer[x21LocalHeight+y21LocalHeight],
              q01_MR_STAR.LocalBuffer(), 
              q01LocalHeight*sizeof(C) );
            mpi::AllReduce
            ( &colSumSendBuffer[0], 
              &colSumRecvBuffer[0],
              reduceSize, mpi::SUM, g.MCComm() );
            std::memcpy
            ( x21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[0], 
              x21LocalHeight*sizeof(C) );
            std::memcpy
            ( y21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[x21LocalHeight], 
              y21LocalHeight*sizeof(C) );
            std::memcpy
            ( q01_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[x21LocalHeight+y21LocalHeight], 
              q01LocalHeight*sizeof(C) );
        }

        Gemv
        ( NORMAL, 
          (C)-1, A02T.LockedLocalMatrix(),
                 x21_MR_STAR.LockedLocalMatrix(),
          (C)+1, p01T_MC_STAR.LocalMatrix() );
        Gemv
        ( NORMAL, 
          (C)-1, W02T.LockedLocalMatrix(),
                 y21_MR_STAR.LockedLocalMatrix(),
          (C)+1, p01T_MC_STAR.LocalMatrix() );

        if( W00.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // combine the Reduce to one of p01[MC,* ] with the redistribution 
            // of q01[MR,* ] -> q01[MC,MR] to the next process column.
            const int localHeight = p01_MC_STAR.LocalHeight();
            std::vector<C> reduceToOneSendBuffer(2*localHeight);
            std::vector<C> reduceToOneRecvBuffer(2*localHeight);

            // Pack p01[MC,* ]
            std::memcpy
            ( &reduceToOneSendBuffer[0], 
              p01_MC_STAR.LocalBuffer(),
              localHeight*sizeof(C) );

            // Fill in contributions to q01[MC,MR] from q01[MR,* ]
            const bool contributing = 
                ( q01_MR_STAR.ColShift() % g.GCD() ==
                  p01_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    std::memcpy
                    ( &reduceToOneSendBuffer[localHeight],
                      q01_MR_STAR.LocalBuffer(), 
                      localHeight*sizeof(C) );
                }
                else
                {
                    // Zero the entire buffer first
                    std::memset
                    ( &reduceToOneSendBuffer[localHeight], 0,
                      localHeight*sizeof(C) );
                    // Fill in the entries that we contribute to.
                    // We seek to find the minimum s in N such that
                    //   s*c = a0-b0 (mod r)
                    // where a0 is the column shift of MC, b0 is the row shift
                    // of MR, and s is our first local entry of MR that will 
                    // contribute to MC. I cannot think of an O(1) method, so
                    // I will instead use the worst-case O(lcm(c,r)/c) method.
                    const int sourcePeriod = g.LCM() / c;
                    const int targetPeriod = g.LCM() / r;
                    const int a0 = p01_MC_STAR.ColShift();
                    const int b0 = q01_MR_STAR.ColShift();

                    int sourceStart = 0;
                    const int f = (r+a0-b0) % r;
                    for( int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    const int globalShift = b0+sourceStart*c;
                    const int targetStart = (globalShift-a0)/r;
                    const int localLength =
                        LocalLength(localHeight,targetStart,targetPeriod);
                    const C* q01_MR_STAR_LocalBuffer = 
                        q01_MR_STAR.LocalBuffer();
                    const int offset = localHeight + targetStart;
                    for( int i=0; i<localLength; ++i )                        
                        reduceToOneSendBuffer[offset+i*targetPeriod] = 
                            q01_MR_STAR_LocalBuffer[sourceStart+i*sourcePeriod];
                }
            }
            else
            {
                std::memset
                ( &reduceToOneSendBuffer[localHeight], 0, 
                  localHeight*sizeof(C) );
            }

            const int nextProcessRow = (alpha11.ColAlignment()+r-1) % r;
            const int nextProcessCol = (alpha11.RowAlignment()+c-1) % c;
            mpi::Reduce
            ( &reduceToOneSendBuffer[0], 
              &reduceToOneRecvBuffer[0],
              2*localHeight, mpi::SUM, nextProcessCol, g.MRComm() );
            if( g.MRRank() == nextProcessCol )
            {
                // Combine the second half into the first half        
                for( int i=0; i<localHeight; ++i )
                    reduceToOneRecvBuffer[i] +=
                        reduceToOneRecvBuffer[i+localHeight];

                // Finish computing w01. During its computation, ensure that 
                // every process has a copy of the last element of the w01.
                // We know a priori that the last element of a01 is one.
                const C* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
                C myDotProduct = blas::Dot
                    ( localHeight, &reduceToOneRecvBuffer[0], 1, 
                                   &a01_MC_STAR_LocalBuffer[0], 1 );
                C sendBuffer[2];
                C recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.MCRank()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[localHeight-1] : 0 );
                mpi::AllReduce
                ( sendBuffer, recvBuffer, 2, mpi::SUM, g.MCComm() );
                C dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w01LastLocalBuffer
                // - w01LastBottomEntry
                C scale = static_cast<C>(0.5)*dotProduct*Conj(tau);
                for( int i=0; i<localHeight; ++i )
                    w01LastLocalBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a01_MC_STAR_LocalBuffer[i] );
                w01LastBottomEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w01[MC,* ] and w01[MR,* ] so that we may place them into W[MC,* ]
            // and W[MR,* ]
            const int localHeight = p01_MC_STAR.LocalHeight();
            std::vector<C> allReduceSendBuffer(2*localHeight);
            std::vector<C> allReduceRecvBuffer(2*localHeight);

            // Pack p01[MC,* ]
            std::memcpy
            ( &allReduceSendBuffer[0], 
              p01_MC_STAR.LocalBuffer(),
              localHeight*sizeof(C) );

            // Fill in contributions to q01[MC,* ] from q01[MR,* ]
            const bool contributing = 
                ( q01_MR_STAR.ColShift() % g.GCD() ==
                  p01_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    std::memcpy
                    ( &allReduceSendBuffer[localHeight],
                      q01_MR_STAR.LocalBuffer(), 
                      localHeight*sizeof(C) );
                }
                else
                {
                    // Zero the entire buffer first
                    std::memset
                    ( &allReduceSendBuffer[localHeight], 0, 
                      localHeight*sizeof(C) );
                    // Fill in the entries that we contribute to.
                    // We seek to find the minimum s in N such that
                    //   s*c = a0-b0 (mod r)
                    // where a0 is the column shift of MC, b0 is the row shift
                    // of MR, and s is our first local entry of MR that will 
                    // contribute to MC. I cannot think of an O(1) method, so
                    // I will instead use the worst-case O(lcm(c,r)/c) method.
                    const int sourcePeriod = g.LCM() / c;
                    const int targetPeriod = g.LCM() / r;
                    const int a0 = p01_MC_STAR.ColShift();
                    const int b0 = q01_MR_STAR.ColShift();

                    int sourceStart = 0;
                    const int f = (r+a0-b0) % r;
                    for( int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    const int globalShift = b0+sourceStart*c;
                    const int targetStart = (globalShift-a0)/r;
                    const int localLength = 
                        LocalLength(localHeight,targetStart,targetPeriod);
                    const C* q01_MR_STAR_LocalBuffer = 
                        q01_MR_STAR.LocalBuffer();
                    const int offset = localHeight + targetStart;
                    for( int i=0; i<localLength; ++i )
                        allReduceSendBuffer[offset+i*targetPeriod] = 
                            q01_MR_STAR_LocalBuffer[sourceStart+i*sourcePeriod];
                }
            }
            else
            {
                std::memset
                ( &allReduceSendBuffer[localHeight], 0, 
                  localHeight*sizeof(C) );
            }

            mpi::AllReduce
            ( &allReduceSendBuffer[0], 
              &allReduceRecvBuffer[0],
              2*localHeight, mpi::SUM, g.MRComm() );

            // Combine the second half into the first half        
            for( int i=0; i<localHeight; ++i )
                allReduceRecvBuffer[i] += allReduceRecvBuffer[i+localHeight];
 
            // Finish computing w01. During its computation, ensure that 
            // every process has a copy of the last element of the w01.
            // We know a priori that the last element of a01 is one.
            const C* a01_MC_STAR_LocalBuffer = a01_MC_STAR.LocalBuffer();
            C myDotProduct = blas::Dot
                ( localHeight, &allReduceRecvBuffer[0], 1, 
                               a01_MC_STAR_LocalBuffer, 1 );
            C dotProduct;
            mpi::AllReduce
            ( &myDotProduct, &dotProduct, 1, mpi::SUM, g.MCComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<C,MC,STAR> w01_MC_STAR(g);
            DistMatrix<C,MR,STAR> w01_MR_STAR(g);
            w01_MC_STAR.View( W_MC_STAR, 0, W00.Width(), w01.Height(), 1 );
            w01_MR_STAR.View( W_MR_STAR, 0, W00.Width(), w01.Height(), 1 );

            // Store w01[MC,* ]
            C scale = static_cast<C>(0.5)*dotProduct*Conj(tau);
            C* w01_MC_STAR_LocalBuffer = w01_MC_STAR.LocalBuffer();
            for( int i=0; i<localHeight; ++i )
                w01_MC_STAR_LocalBuffer[i] = tau*
                    ( allReduceRecvBuffer[i]-
                      scale*a01_MC_STAR_LocalBuffer[i] );

            // Form w01[MR,* ]
            w01_MR_STAR = w01_MC_STAR;
        }
        //--------------------------------------------------------------------//
        a01_MC_STAR.FreeAlignments();
        a01_MR_STAR.FreeAlignments();
        p01_MC_STAR.FreeAlignments();
        p01_MR_STAR.FreeAlignments();
        q01_MC_STAR.FreeAlignments();
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

    expandedABR.SetRealDiagonal( e, 1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elemental
