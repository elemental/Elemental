/*
   Copyright (c) 2009-2011, Jack Poulson
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

template<typename R> // representation of a real number
inline void
elemental::advanced::internal::HermitianPanelTridiagL
( DistMatrix<R,MC,MR  >& A,
  DistMatrix<R,MC,MR  >& W,
  DistMatrix<R,MC,STAR>& APan_MC_STAR, 
  DistMatrix<R,MR,STAR>& APan_MR_STAR,
  DistMatrix<R,MC,STAR>& W_MC_STAR,
  DistMatrix<R,MR,STAR>& W_MR_STAR )
{
    const int panelSize = W.Width();
    const int bottomSize = W.Height()-panelSize;
#ifndef RELEASE
    PushCallStack("advanced::internal::HermitianPanelTridiagL");
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
        throw std::logic_error("W and A must be aligned.");
#endif
    const Grid& g = A.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();

    // Create a distributed matrix for storing the subdiagonal
    DistMatrix<R,MD,STAR> e(g);
    e.AlignWithDiagonal( A, -1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<R,MC,MR> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g),  alpha21T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),            a21B(g),
                         A20(g), a21(g),     A22(g),  A20B(g);
    DistMatrix<R,MC,MR> 
        WTL(g), WTR(g),  W00(g), w01(g),     W02(g),  WCol(g),
        WBL(g), WBR(g),  w10(g), omega11(g), w12(g),
                         W20(g), w21(g),     W22(g),  W20B(g), w21Last(g);
    DistMatrix<R,MD,STAR> eT(g),  e0(g),
                          eB(g),  epsilon1(g),
                                  e2(g);

    // Temporary distributions
    std::vector<R> w21LastLocalBuffer(A.Height()/r+1);
    DistMatrix<R,MC,STAR> a21_MC_STAR(g);
    DistMatrix<R,MC,STAR> a21B_MC_STAR(g);
    DistMatrix<R,MR,STAR> a21_MR_STAR(g);
    DistMatrix<R,MC,STAR> p21_MC_STAR(g);
    DistMatrix<R,MC,STAR> p21B_MC_STAR(g);
    DistMatrix<R,MR,STAR> p21_MR_STAR(g);
    DistMatrix<R,MC,STAR> q21_MC_STAR(g);
    DistMatrix<R,MR,STAR> q21_MR_STAR(g);
    DistMatrix<R,MR,STAR> x01_MR_STAR(g);
    DistMatrix<R,MR,STAR> y01_MR_STAR(g);
    DistMatrix<R,MC,STAR> a21Last_MC_STAR(g);
    DistMatrix<R,MR,STAR> a21Last_MR_STAR(g);
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

        ACol.View2x1
        ( alpha11,
          a21 );

        WCol.View2x1
        ( omega11,
          w21 );

        // View the portions of A20 and W20 outside of this panel's square
        A20B.View( A, panelSize, 0, bottomSize, A20.Width() );
        W20B.View( W, panelSize, 0, bottomSize, W20.Width() );

        if( !firstIteration )
        {
            a21Last_MC_STAR.View
            ( APan_MC_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            a21Last_MR_STAR.View
            ( APan_MR_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            w21Last.View
            ( W, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
        }
            
        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        a21_MC_STAR.AlignWith( A22 );
        a21_MR_STAR.AlignWith( A22 );
        p21_MC_STAR.AlignWith( A22 );
        p21_MR_STAR.AlignWith( A22 );
        q21_MC_STAR.AlignWith( A22 );
        q21_MR_STAR.AlignWith( A22 );
        a21_MC_STAR.ResizeTo( a21.Height(), 1 );
        a21_MR_STAR.ResizeTo( a21.Height(), 1 );
        p21_MC_STAR.ResizeTo( a21.Height(), 1 );
        p21_MR_STAR.ResizeTo( a21.Height(), 1 );
        q21_MC_STAR.ResizeTo( a21.Height(), 1 );
        q21_MR_STAR.ResizeTo( a21.Height(), 1 );
        x01_MR_STAR.AlignWith( W20B );
        y01_MR_STAR.AlignWith( W20B );
        x01_MR_STAR.ResizeTo( W20B.Width(), 1 );
        y01_MR_STAR.ResizeTo( W20B.Width(), 1 );

        // View the portions of a21[MC,* ] and p21[MC,* ] below the current
        // panel's square
        a21B_MC_STAR.View
        ( a21_MC_STAR, a21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        p21B_MC_STAR.View
        ( p21_MC_STAR, p21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.MRRank() == alpha11.RowAlignment() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                int AColLocalHeight = ACol.LocalHeight();
                R* AColLocalBuffer = ACol.LocalBuffer();
                const R* a21Last_MC_STAR_LocalBuffer = 
                    a21Last_MC_STAR.LocalBuffer();
                for( int i=0; i<AColLocalHeight; ++i )
                    AColLocalBuffer[i] -=
                        w21LastLocalBuffer[i] + 
                        a21Last_MC_STAR_LocalBuffer[i]*w21LastFirstEntry;
            }
        }
        if( thisIsMyCol )
        {
            // Compute the Householder reflector
            tau = advanced::internal::ColReflector( alpha21T, a21B );
        }
        // Store the subdiagonal value and turn a21 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry
        alpha21T.GetDiagonal( epsilon1 );
        alpha21T.Set( 0, 0, (R)1 );

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
                std::memcpy
                ( &rowBroadcastBuffer[0], 
                  a21.LocalBuffer(), 
                  a21LocalHeight*sizeof(R) );
                rowBroadcastBuffer[a21LocalHeight] = tau;
            }
            // Broadcast a21 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+1, a21.RowAlignment(), g.MRComm() );
            // Store a21[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            std::memcpy
            ( a21_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0],
              a21LocalHeight*sizeof(R) );
            // Store a21[MC,* ] into APan[MC,* ]
            int APan_MC_STAR_Offset = APan_MC_STAR.LocalHeight()-a21LocalHeight;
            std::memcpy
            ( APan_MC_STAR.LocalBuffer(APan_MC_STAR_Offset,0), 
              &rowBroadcastBuffer[0],
              (APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset)*sizeof(R) );
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight];
            
            a21_MR_STAR = a21_MC_STAR;
            // Store a21[MR,* ]
            int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            std::memcpy
            ( APan_MR_STAR.LocalBuffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.LocalBuffer(),
              (APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset)*sizeof(R) );
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
                std::memcpy
                ( &rowBroadcastBuffer[0], 
                  a21.LocalBuffer(),
                  a21LocalHeight*sizeof(R) );
                std::memcpy
                ( &rowBroadcastBuffer[a21LocalHeight], 
                  &w21LastLocalBuffer[0],
                  w21LastLocalHeight*sizeof(R) );
                rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight] = tau;
            }
            // Broadcast a21, w21Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+w21LastLocalHeight+1, 
              a21.RowAlignment(), g.MRComm() );
            // Store a21[MC,* ] into its DistMatrix class 
            std::memcpy
            ( a21_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0],
              a21LocalHeight*sizeof(R) );
            // Store a21[MC,* ] into APan[MC,* ]
            int APan_MC_STAR_Offset = APan_MC_STAR.LocalHeight()-a21LocalHeight;
            std::memcpy
            ( APan_MC_STAR.LocalBuffer(APan_MC_STAR_Offset,A00.Width()), 
              &rowBroadcastBuffer[0],
              (APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset)*sizeof(R) );
            // Store w21Last[MC,* ] into its DistMatrix class
            w21Last_MC_STAR.AlignWith( alpha11 );
            w21Last_MC_STAR.ResizeTo( a21.Height()+1, 1 );
            std::memcpy
            ( w21Last_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[a21LocalHeight], 
              w21LastLocalHeight*sizeof(R) );
            // Store the bottom part of w21Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w21.
            int W_MC_STAR_Offset = W_MC_STAR.LocalHeight()-w21LastLocalHeight;
            std::memcpy
            ( W_MC_STAR.LocalBuffer(W_MC_STAR_Offset,A00.Width()-1),
              &rowBroadcastBuffer[a21LocalHeight],
              (W_MC_STAR.LocalHeight()-W_MC_STAR_Offset)*sizeof(R) );
            if( g.MRRank() == w21Last.RowAlignment() )
            {
                std::memcpy
                ( w21Last.LocalBuffer(),
                  &rowBroadcastBuffer[a21LocalHeight],
                  w21LastLocalHeight*sizeof(R) );
            }
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight];

            // Form a21[MR,* ] and w21Last[MR,* ] by combining the 
            // communications needed for taking a vector from 
            // [MC,* ] -> [MR,* ]: 
            //   local copy to [VC,* ], 
            //   Send/Recv to [VR,* ], 
            //   AllGather to [MR,* ]
            // We can combine the two by treating a21 as [0; a21] so that its 
            // alignments are equivalent to w21Last.

            // alpha11 is the top-left entry of the A22 from the last iteration,
            // so its alignments are identical
            int colAlignSource = alpha11.ColAlignment();
            int colAlignDest = alpha11.RowAlignment();
            int colShiftSource = alpha11.ColShift();
            int colShiftDest = alpha11.RowShift();

            int height = a21.Height()+1;
            int portionSize = 
                std::max(2*MaxLocalLength(height,p),mpi::MIN_COLL_MSG);

            int colShiftVRDest = Shift(g.VRRank(),colAlignDest,p);
            int colShiftVCSource = Shift(g.VCRank(),colAlignSource,p);
            int sendRankRM = (g.VRRank()+(p+colShiftVCSource-colShiftVRDest))%p;
            int recvRankCM = (g.VCRank()+(p+colShiftVRDest-colShiftVCSource))%p;
            int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

            std::vector<R> transposeBuffer( (r+1)*portionSize );
            R* sendBuf = &transposeBuffer[0];
            R* recvBuf = &transposeBuffer[r*portionSize];

            // (w21Last[VC,* ] <- w21Last[MC,* ]) and
            // ([0; a21][VC,* ] <- [0; a21][MC,* ])
            {
                // Pack the necessary portion of w21Last[MC,* ]
                int w21Shift = Shift(g.VCRank(),colAlignSource,p);
                int w21Offset = (w21Shift-colShiftSource)/r;
                int w21LocalHeight = LocalLength(height,w21Shift,p);
                const R* w21LastBuffer = 
                    w21Last_MC_STAR.LocalBuffer(w21Offset,0);
                for( int i=0; i<w21LocalHeight; ++i )
                    sendBuf[i] = w21LastBuffer[i*c];
                
                // Pack the necessary portion of a21[MC,* ]
                int a21Shift = (w21Shift+p-1) % p;
                int a21Offset = (a21Shift-((colShiftSource+r-1)%r))/r;
                int a21LocalHeight = LocalLength(height-1,a21Shift,p);
                const R* a21Buffer = a21_MC_STAR.LocalBuffer(a21Offset,0);
                for( int i=0; i<a21LocalHeight; ++i )
                    sendBuf[w21LocalHeight+i] = a21Buffer[i*c];
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
            w21Last_MR_STAR.AlignWith( alpha11 );
            w21Last_MR_STAR.ResizeTo( a21.Height()+1, 1 );
            for( int k=0; k<r; ++k )
            {
                // Unpack into w21Last[MR,* ]
                const R* w21Data = &sendBuf[k*portionSize];
                int w21Shift = Shift(g.MRRank()+c*k,colAlignDest,p);
                int w21Offset = (w21Shift-colShiftDest) / c;
                int w21LocalHeight = LocalLength(height,w21Shift,p);
                R* w21LastBuffer = w21Last_MR_STAR.LocalBuffer(w21Offset,0);
                for( int i=0; i<w21LocalHeight; ++i )
                    w21LastBuffer[i*r] = w21Data[i];

                // Unpack into a21[MR,* ]
                const R* a21Data = &sendBuf[k*portionSize+w21LocalHeight];
                int a21Shift = (w21Shift+p-1) % p;
                int a21Offset = (a21Shift-((colShiftDest+c-1)%c))/c;
                int a21LocalHeight = LocalLength(height-1,a21Shift,p);
                R* a21Buffer = a21_MR_STAR.LocalBuffer(a21Offset,0);
                for( int i=0; i<a21LocalHeight; ++i )
                    a21Buffer[i*r] = a21Data[i];
            }
            // Store w21Last[MR,* ]
            int W_MR_STAR_Offset = 
                W_MR_STAR.LocalHeight()-w21Last_MR_STAR.LocalHeight();
            std::memcpy
            ( W_MR_STAR.LocalBuffer(W_MR_STAR_Offset,A00.Width()-1),
              w21Last_MR_STAR.LocalBuffer(),
              (W_MR_STAR.LocalHeight()-W_MR_STAR_Offset)*sizeof(R) );
            // Store a21[MR,* ]
            int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            std::memcpy
            ( APan_MR_STAR.LocalBuffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.LocalBuffer(),
              (APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset)*sizeof(R) );

            // Update the portion of A22 that is in our current panel with 
            // w21Last and a21Last using two gers. We do not need their top 
            // entries. We trash the upper triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            DistMatrix<R,MC,STAR> a21Last_MC_STAR_Bottom(g);
            DistMatrix<R,MR,STAR> a21Last_MR_STAR_Bottom(g);
            DistMatrix<R,MC,STAR> w21Last_MC_STAR_Bottom(g);
            DistMatrix<R,MR,STAR> w21Last_MR_STAR_Bottom(g);
            a21Last_MC_STAR_Bottom.View
            ( a21Last_MC_STAR, 1, 0, a21Last_MC_STAR.Height()-1, 1 );
            a21Last_MR_STAR_Bottom.View
            ( a21Last_MR_STAR, 1, 0, a21Last_MR_STAR.Height()-1, 1 );
            w21Last_MC_STAR_Bottom.View
            ( w21Last_MC_STAR, 1, 0, w21Last_MC_STAR.Height()-1, 1 );
            w21Last_MR_STAR_Bottom.View
            ( w21Last_MR_STAR, 1, 0, w21Last_MR_STAR.Height()-1, 1 );
            const R* a21_MC_STAR_Buffer = a21Last_MC_STAR_Bottom.LocalBuffer();
            const R* a21_MR_STAR_Buffer = a21Last_MR_STAR_Bottom.LocalBuffer();
            const R* w21_MC_STAR_Buffer = w21Last_MC_STAR_Bottom.LocalBuffer();
            const R* w21_MR_STAR_Buffer = w21Last_MR_STAR_Bottom.LocalBuffer();
            R* A22Buffer = A22.LocalBuffer();
            int localHeight = W22.LocalHeight();
            int localWidth = W22.LocalWidth();
            int lDim = A22.LocalLDim();
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
        PopBlocksizeStack();
        p21_MC_STAR.SetToZero();
        q21_MR_STAR.SetToZero();
        basic::internal::LocalSymvColAccumulateL
        ( (R)1, A22, a21_MC_STAR, a21_MR_STAR, p21_MC_STAR, q21_MR_STAR );
        PushBlocksizeStack( 1 );

        basic::Gemv
        ( TRANSPOSE, 
          (R)1, W20B.LockedLocalMatrix(),
                a21B_MC_STAR.LockedLocalMatrix(),
          (R)0, x01_MR_STAR.LocalMatrix() );
        basic::Gemv
        ( TRANSPOSE, 
          (R)1, A20B.LockedLocalMatrix(),
                a21B_MC_STAR.LockedLocalMatrix(),
          (R)0, y01_MR_STAR.LocalMatrix() );
        // Combine the AllReduce column summations of x01[MR,* ], y01[MR,* ],
        // and q21[MR,* ]
        {
            int x01LocalHeight = x01_MR_STAR.LocalHeight();
            int q21LocalHeight = q21_MR_STAR.LocalHeight();
            std::vector<R> colSumSendBuffer(2*x01LocalHeight+q21LocalHeight);
            std::vector<R> colSumRecvBuffer(2*x01LocalHeight+q21LocalHeight);
            std::memcpy
            ( &colSumSendBuffer[0], 
              x01_MR_STAR.LocalBuffer(), 
              x01LocalHeight*sizeof(R) );
            std::memcpy
            ( &colSumSendBuffer[x01LocalHeight],
              y01_MR_STAR.LocalBuffer(), 
              x01LocalHeight*sizeof(R) );
            std::memcpy
            ( &colSumSendBuffer[2*x01LocalHeight],
              q21_MR_STAR.LocalBuffer(), 
              q21LocalHeight*sizeof(R) );
            mpi::AllReduce
            ( &colSumSendBuffer[0], 
              &colSumRecvBuffer[0],
              2*x01LocalHeight+q21LocalHeight, mpi::SUM, g.MCComm() );
            std::memcpy
            ( x01_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[0], 
              x01LocalHeight*sizeof(R) );
            std::memcpy
            ( y01_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[x01LocalHeight], 
              x01LocalHeight*sizeof(R) );
            std::memcpy
            ( q21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[2*x01LocalHeight], 
              q21LocalHeight*sizeof(R) );
        }

        basic::Gemv
        ( NORMAL, 
          (R)-1, A20B.LockedLocalMatrix(),
                 x01_MR_STAR.LockedLocalMatrix(),
          (R)+1, p21B_MC_STAR.LocalMatrix() );
        basic::Gemv
        ( NORMAL, 
          (R)-1, W20B.LockedLocalMatrix(),
                 y01_MR_STAR.LockedLocalMatrix(),
          (R)+1, p21B_MC_STAR.LocalMatrix() );

        if( W22.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // combine the Reduce to one of p21[MC,* ] with the redistribution 
            // of q21[MR,* ] -> q21[MC,MR] to the next process column.
            int localHeight = p21_MC_STAR.LocalHeight();
            std::vector<R> reduceToOneSendBuffer(2*localHeight);
            std::vector<R> reduceToOneRecvBuffer(2*localHeight);

            // Pack p21[MC,* ]
            std::memcpy
            ( &reduceToOneSendBuffer[0], 
              p21_MC_STAR.LocalBuffer(),
              localHeight*sizeof(R) );

            // Fill in contributions to q21[MC,MR] from q21[MR,* ]
            bool contributing = 
                ( q21_MR_STAR.ColShift() % g.GCD() ==
                  p21_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    std::memcpy
                    ( &reduceToOneSendBuffer[localHeight],
                      q21_MR_STAR.LocalBuffer(), 
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
                    int sourcePeriod = g.LCM() / c;
                    int targetPeriod = g.LCM() / r;
                    int a0 = p21_MC_STAR.ColShift();
                    int b0 = q21_MR_STAR.ColShift();

                    int sourceStart = 0;
                    int f = (r+a0-b0) % r;
                    for( int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    int globalShift = b0+sourceStart*c;
                    int targetStart = (globalShift-a0)/r;
                    int localLength =
                        LocalLength(localHeight,targetStart,targetPeriod);
                    const R* q21_MR_STAR_LocalBuffer = 
                        q21_MR_STAR.LocalBuffer();
                    int offset = localHeight + targetStart;
                    for( int i=0; i<localLength; ++i )                        
                        reduceToOneSendBuffer[offset+i*targetPeriod] = 
                            q21_MR_STAR_LocalBuffer[sourceStart+i*sourcePeriod];
                }
            }
            else
            {
                std::memset
                ( &reduceToOneSendBuffer[localHeight], 0, 
                  localHeight*sizeof(R) );
            }

            int nextProcessRow = (alpha11.ColAlignment()+1) % r;
            int nextProcessCol = (alpha11.RowAlignment()+1) % c;
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

                // Finish computing w21. During its computation, ensure that
                // every process has a copy of the first element of the w21.
                // We know a priori that the first element of a21 is one.
                const R* a21_MC_STAR_LocalBuffer = a21_MC_STAR.LocalBuffer();
                R myDotProduct = blas::Dot
                    ( localHeight, &reduceToOneRecvBuffer[0], 1, 
                                   a21_MC_STAR_LocalBuffer,   1 );
                R sendBuffer[2];
                R recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.MCRank()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[0] : 0 );
                mpi::AllReduce
                ( sendBuffer, recvBuffer, 2, mpi::SUM, g.MCComm() );
                R dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w21LastLocalBuffer
                // - w21LastFirstEntry
                R scale = 0.5*dotProduct*tau;
                for( int i=0; i<localHeight; ++i )
                    w21LastLocalBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a21_MC_STAR_LocalBuffer[i] );
                w21LastFirstEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w21[MC,* ] and w21[MR,* ] so that we may place them into W[MC,* ]
            // and W[MR,* ]
            int localHeight = p21_MC_STAR.LocalHeight();
            std::vector<R> allReduceSendBuffer(2*localHeight);
            std::vector<R> allReduceRecvBuffer(2*localHeight);

            // Pack p21[MC,* ]
            std::memcpy
            ( &allReduceSendBuffer[0], 
              p21_MC_STAR.LocalBuffer(),
              localHeight*sizeof(R) );

            // Fill in contributions to q21[MC,* ] from q21[MR,* ]
            bool contributing = 
                ( q21_MR_STAR.ColShift() % g.GCD() ==
                  p21_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    std::memcpy
                    ( &allReduceSendBuffer[localHeight],
                      q21_MR_STAR.LocalBuffer(), 
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
                    int sourcePeriod = g.LCM() / c;
                    int targetPeriod = g.LCM() / r;
                    int a0 = p21_MC_STAR.ColShift();
                    int b0 = q21_MR_STAR.ColShift();

                    int sourceStart = 0;
                    int f = (r+a0-b0) % r;
                    for( int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    int globalShift = b0+sourceStart*c;
                    int targetStart = (globalShift-a0)/r;
                    int localLength = 
                        LocalLength(localHeight,targetStart,targetPeriod);
                    const R* q21_MR_STAR_LocalBuffer = 
                        q21_MR_STAR.LocalBuffer();
                    int offset = localHeight + targetStart;
                    for( int i=0; i<localLength; ++i )
                        allReduceSendBuffer[offset+i*targetPeriod] = 
                            q21_MR_STAR_LocalBuffer[sourceStart+i*sourcePeriod];
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
 
            // Finish computing w21. 
            const R* a21_MC_STAR_LocalBuffer = a21_MC_STAR.LocalBuffer();
            R myDotProduct = blas::Dot
                ( localHeight, &allReduceRecvBuffer[0], 1, 
                               a21_MC_STAR_LocalBuffer, 1 );
            R dotProduct;
            mpi::AllReduce
            ( &myDotProduct, &dotProduct, 1, mpi::SUM, g.MCComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<R,MC,STAR> w21_MC_STAR(g);
            DistMatrix<R,MR,STAR> w21_MR_STAR(g);
            w21_MC_STAR.View
            ( W_MC_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );
            w21_MR_STAR.View
            ( W_MR_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );

            // Store w21[MC,* ]
            R scale = 0.5*dotProduct*tau;
            R* w21_MC_STAR_LocalBuffer = w21_MC_STAR.LocalBuffer();
            for( int i=0; i<localHeight; ++i )
                w21_MC_STAR_LocalBuffer[i] = tau*
                    ( allReduceRecvBuffer[i]-
                      scale*a21_MC_STAR_LocalBuffer[i] );

            // Form w21[MR,* ]
            w21_MR_STAR = w21_MC_STAR;
        }
        //--------------------------------------------------------------------//
        a21_MC_STAR.FreeAlignments();
        a21_MR_STAR.FreeAlignments();
        p21_MC_STAR.FreeAlignments();
        p21_MR_STAR.FreeAlignments();
        q21_MC_STAR.FreeAlignments();
        q21_MR_STAR.FreeAlignments();
        x01_MR_STAR.FreeAlignments();
        y01_MR_STAR.FreeAlignments();

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
        
        SlidePartitionDown
        ( eT,  e0,
               epsilon1,
         /**/ /********/
          eB,  e2 );
        firstIteration = false;
    }
    PopBlocksizeStack();

    // View the portion of A that e is the subdiagonal of, then place e into it
    DistMatrix<R,MC,MR> expandedATL(g);
    expandedATL.View( A, 0, 0, panelSize+1, panelSize+1 );
    expandedATL.SetDiagonal( e, -1 );
#ifndef RELEASE
    PopCallStack();
#endif
}

#ifndef WITHOUT_COMPLEX
template<typename R> // representation of a real number
inline void
elemental::advanced::internal::HermitianPanelTridiagL
( DistMatrix<std::complex<R>,MC,MR  >& A,
  DistMatrix<std::complex<R>,MC,MR  >& W,
  DistMatrix<std::complex<R>,MD,STAR>& t,
  DistMatrix<std::complex<R>,MC,STAR>& APan_MC_STAR, 
  DistMatrix<std::complex<R>,MR,STAR>& APan_MR_STAR,
  DistMatrix<std::complex<R>,MC,STAR>& W_MC_STAR,
  DistMatrix<std::complex<R>,MR,STAR>& W_MR_STAR )
{
    const int panelSize = W.Width();
    const int bottomSize = W.Height()-panelSize;
#ifndef RELEASE
    PushCallStack("advanced::internal::HermitianPanelTridiagL");
    if( A.Grid() != W.Grid() || W.Grid() != t.Grid() )
        throw std::logic_error
        ("A, W, and t must be distributed over the same grid.");
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
        throw std::logic_error("t is not aligned with A's subdiagonal.");
#endif
    typedef std::complex<R> C;

    const Grid& g = A.Grid();
    const int r = g.Height();
    const int c = g.Width();
    const int p = g.Size();

    // Create a distributed matrix for storing the subdiagonal
    DistMatrix<R,MD,STAR> e(g);
    e.AlignWithDiagonal( A, -1 );
    e.ResizeTo( panelSize, 1 );

    // Matrix views 
    DistMatrix<C,MC,MR> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  ACol(g),  alpha21T(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),            a21B(g),
                         A20(g), a21(g),     A22(g),  A20B(g);
    DistMatrix<C,MC,MR> 
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
    std::vector<C> w21LastLocalBuffer(A.Height()/r+1);
    DistMatrix<C,MC,STAR> a21_MC_STAR(g);
    DistMatrix<C,MC,STAR> a21B_MC_STAR(g);
    DistMatrix<C,MR,STAR> a21_MR_STAR(g);
    DistMatrix<C,MC,STAR> p21_MC_STAR(g);
    DistMatrix<C,MC,STAR> p21B_MC_STAR(g);
    DistMatrix<C,MR,STAR> p21_MR_STAR(g);
    DistMatrix<C,MC,STAR> q21_MC_STAR(g);
    DistMatrix<C,MR,STAR> q21_MR_STAR(g);
    DistMatrix<C,MR,STAR> x01_MR_STAR(g);
    DistMatrix<C,MR,STAR> y01_MR_STAR(g);
    DistMatrix<C,MC,STAR> a21Last_MC_STAR(g);
    DistMatrix<C,MR,STAR> a21Last_MR_STAR(g);
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

        ACol.View2x1
        ( alpha11,
          a21 );

        WCol.View2x1
        ( omega11,
          w21 );

        // View the portions of A20 and W20 outside of this panel's square
        A20B.View( A, panelSize, 0, bottomSize, A20.Width() );
        W20B.View( W, panelSize, 0, bottomSize, W20.Width() );

        if( !firstIteration )
        {
            a21Last_MC_STAR.View
            ( APan_MC_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            a21Last_MR_STAR.View
            ( APan_MR_STAR, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
            w21Last.View
            ( W, WTL.Height(), WTL.Width()-1, WBL.Height(), 1 );
        }
            
        PartitionDown
        ( a21, alpha21T,
               a21B,     1 );

        a21_MC_STAR.AlignWith( A22 );
        a21_MR_STAR.AlignWith( A22 );
        p21_MC_STAR.AlignWith( A22 );
        p21_MR_STAR.AlignWith( A22 );
        q21_MC_STAR.AlignWith( A22 );
        q21_MR_STAR.AlignWith( A22 );
        a21_MC_STAR.ResizeTo( a21.Height(), 1 );
        a21_MR_STAR.ResizeTo( a21.Height(), 1 );
        p21_MC_STAR.ResizeTo( a21.Height(), 1 );
        p21_MR_STAR.ResizeTo( a21.Height(), 1 );
        q21_MC_STAR.ResizeTo( a21.Height(), 1 );
        q21_MR_STAR.ResizeTo( a21.Height(), 1 );
        x01_MR_STAR.AlignWith( W20B );
        y01_MR_STAR.AlignWith( W20B );
        x01_MR_STAR.ResizeTo( W20B.Width(), 1 );
        y01_MR_STAR.ResizeTo( W20B.Width(), 1 );

        // View the portions of a21[MC,* ] and p21[MC,* ] below the current
        // panel's square
        a21B_MC_STAR.View
        ( a21_MC_STAR, a21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        p21B_MC_STAR.View
        ( p21_MC_STAR, p21_MC_STAR.Height()-bottomSize, 0, bottomSize, 1 );
        //--------------------------------------------------------------------//
        const bool thisIsMyCol = ( g.MRRank() == alpha11.RowAlignment() );
        if( thisIsMyCol )
        {
            if( !firstIteration )
            {
                // Finish updating the current column with two axpy's
                int AColLocalHeight = ACol.LocalHeight();
                C* AColLocalBuffer = ACol.LocalBuffer();
                const C* a21Last_MC_STAR_LocalBuffer = 
                    a21Last_MC_STAR.LocalBuffer();
                for( int i=0; i<AColLocalHeight; ++i )
                    AColLocalBuffer[i] -=
                        w21LastLocalBuffer[i] + 
                        a21Last_MC_STAR_LocalBuffer[i]*Conj(w21LastFirstEntry);
            }
        }
        if( thisIsMyCol )
        {
            // Compute the Householder reflector
            tau = advanced::internal::ColReflector( alpha21T, a21B );
            if( g.MCRank() == alpha21T.ColAlignment() )
                tau1.SetLocalEntry(0,0,tau);
        }
        // Store the subdiagonal value and turn a21 into a proper scaled 
        // reflector by explicitly placing the implicit one in its first entry.
        alpha21T.GetRealDiagonal( epsilon1 );
        alpha21T.Set( 0, 0, (C)1 );

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
                std::memcpy
                ( &rowBroadcastBuffer[0], 
                  a21.LocalBuffer(), 
                  a21LocalHeight*sizeof(C) );
                rowBroadcastBuffer[a21LocalHeight] = tau;
            }
            // Broadcast a21 and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+1, a21.RowAlignment(), g.MRComm() );
            // Store a21[MC,* ] into its DistMatrix class and also store a copy
            // for the next iteration
            std::memcpy
            ( a21_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0],
              a21LocalHeight*sizeof(C) );
            // Store a21[MC,* ] into APan[MC,* ]
            int APan_MC_STAR_Offset = APan_MC_STAR.LocalHeight()-a21LocalHeight;
            std::memcpy
            ( APan_MC_STAR.LocalBuffer(APan_MC_STAR_Offset,0), 
              &rowBroadcastBuffer[0],
              (APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset)*sizeof(C) );
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight];
            
            a21_MR_STAR = a21_MC_STAR;
            // Store a21[MR,* ]
            int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            std::memcpy
            ( APan_MR_STAR.LocalBuffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.LocalBuffer(),
              (APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset)*sizeof(C) );
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
                std::memcpy
                ( &rowBroadcastBuffer[0], 
                  a21.LocalBuffer(),
                  a21LocalHeight*sizeof(C) );
                std::memcpy
                ( &rowBroadcastBuffer[a21LocalHeight], 
                  &w21LastLocalBuffer[0],
                  w21LastLocalHeight*sizeof(C) );
                rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight] = tau;
            }
            // Broadcast a21, w21Last, and tau across the process row
            mpi::Broadcast
            ( &rowBroadcastBuffer[0], 
              a21LocalHeight+w21LastLocalHeight+1, 
              a21.RowAlignment(), g.MRComm() );
            // Store a21[MC,* ] into its DistMatrix class 
            std::memcpy
            ( a21_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[0],
              a21LocalHeight*sizeof(C) );
            // Store a21[MC,* ] into APan[MC,* ]
            int APan_MC_STAR_Offset = APan_MC_STAR.LocalHeight()-a21LocalHeight;
            std::memcpy
            ( APan_MC_STAR.LocalBuffer(APan_MC_STAR_Offset,A00.Width()), 
              &rowBroadcastBuffer[0],
              (APan_MC_STAR.LocalHeight()-APan_MC_STAR_Offset)*sizeof(C) );
            // Store w21Last[MC,* ] into its DistMatrix class
            w21Last_MC_STAR.AlignWith( alpha11 );
            w21Last_MC_STAR.ResizeTo( a21.Height()+1, 1 );
            std::memcpy
            ( w21Last_MC_STAR.LocalBuffer(), 
              &rowBroadcastBuffer[a21LocalHeight], 
              w21LastLocalHeight*sizeof(C) );
            // Store the bottom part of w21Last[MC,* ] into WB[MC,* ] and, 
            // if necessary, w21.
            int W_MC_STAR_Offset = W_MC_STAR.LocalHeight()-w21LastLocalHeight;
            std::memcpy
            ( W_MC_STAR.LocalBuffer(W_MC_STAR_Offset,A00.Width()-1),
              &rowBroadcastBuffer[a21LocalHeight],
              (W_MC_STAR.LocalHeight()-W_MC_STAR_Offset)*sizeof(C) );
            if( g.MRRank() == w21Last.RowAlignment() )
            {
                std::memcpy
                ( w21Last.LocalBuffer(),
                  &rowBroadcastBuffer[a21LocalHeight],
                  w21LastLocalHeight*sizeof(C) );
            }
            // Store tau
            tau = rowBroadcastBuffer[a21LocalHeight+w21LastLocalHeight];

            // Form a21[MR,* ] and w21Last[MR,* ] by combining the 
            // communications needed for taking a vector from 
            // [MC,* ] -> [MR,* ]: 
            //   local copy to [VC,* ], 
            //   Send/Recv to [VR,* ], 
            //   AllGather to [MR,* ]
            // We can combine the two by treating a21 as [0; a21] so that its 
            // alignments are equivalent to w21Last.

            // alpha11 is the top-left entry of the A22 from the last iteration,
            // so its alignments are identical
            int colAlignSource = alpha11.ColAlignment();
            int colAlignDest = alpha11.RowAlignment();
            int colShiftSource = alpha11.ColShift();
            int colShiftDest = alpha11.RowShift();

            int height = a21.Height()+1;
            int portionSize = 
                std::max(2*MaxLocalLength(height,p),mpi::MIN_COLL_MSG);

            int colShiftVRDest = Shift(g.VRRank(),colAlignDest,p);
            int colShiftVCSource = Shift(g.VCRank(),colAlignSource,p);
            int sendRankRM = (g.VRRank()+(p+colShiftVCSource-colShiftVRDest))%p;
            int recvRankCM = (g.VCRank()+(p+colShiftVRDest-colShiftVCSource))%p;
            int recvRankRM = (recvRankCM/r)+c*(recvRankCM%r);

            std::vector<C> transposeBuffer( (r+1)*portionSize );
            C* sendBuf = &transposeBuffer[0];
            C* recvBuf = &transposeBuffer[r*portionSize];

            // (w21Last[VC,* ] <- w21Last[MC,* ]) and
            // ([0; a21][VC,* ] <- [0; a21][MC,* ])
            {
                // Pack the necessary portion of w21Last[MC,* ]
                int w21Shift = Shift(g.VCRank(),colAlignSource,p);
                int w21Offset = (w21Shift-colShiftSource)/r;
                int w21LocalHeight = LocalLength(height,w21Shift,p);
                const C* w21LastBuffer = 
                    w21Last_MC_STAR.LocalBuffer(w21Offset,0);
                for( int i=0; i<w21LocalHeight; ++i )
                    sendBuf[i] = w21LastBuffer[i*c];
                
                // Pack the necessary portion of a21[MC,* ]
                int a21Shift = (w21Shift+p-1) % p;
                int a21Offset = (a21Shift-((colShiftSource+r-1)%r))/r;
                int a21LocalHeight = LocalLength(height-1,a21Shift,p);
                const C* a21Buffer = a21_MC_STAR.LocalBuffer(a21Offset,0);
                for( int i=0; i<a21LocalHeight; ++i )
                    sendBuf[w21LocalHeight+i] = a21Buffer[i*c];
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
            w21Last_MR_STAR.AlignWith( alpha11 );
            w21Last_MR_STAR.ResizeTo( a21.Height()+1, 1 );
            for( int k=0; k<r; ++k )
            {
                // Unpack into w21Last[MR,* ]
                const C* w21Data = &sendBuf[k*portionSize];
                int w21Shift = Shift(g.MRRank()+c*k,colAlignDest,p);
                int w21Offset = (w21Shift-colShiftDest) / c;
                int w21LocalHeight = LocalLength(height,w21Shift,p);
                C* w21LastBuffer = w21Last_MR_STAR.LocalBuffer(w21Offset,0);
                for( int i=0; i<w21LocalHeight; ++i )
                    w21LastBuffer[i*r] = w21Data[i];

                // Unpack into a21[MR,* ]
                const C* a21Data = &sendBuf[k*portionSize+w21LocalHeight];
                int a21Shift = (w21Shift+p-1) % p;
                int a21Offset = (a21Shift-((colShiftDest+c-1)%c))/c;
                int a21LocalHeight = LocalLength(height-1,a21Shift,p);
                C* a21Buffer = a21_MR_STAR.LocalBuffer(a21Offset,0);
                for( int i=0; i<a21LocalHeight; ++i )
                    a21Buffer[i*r] = a21Data[i];
            }
            // Store w21Last[MR,* ]
            int W_MR_STAR_Offset = 
                W_MR_STAR.LocalHeight()-w21Last_MR_STAR.LocalHeight();
            std::memcpy
            ( W_MR_STAR.LocalBuffer(W_MR_STAR_Offset,A00.Width()-1),
              w21Last_MR_STAR.LocalBuffer(),
              (W_MR_STAR.LocalHeight()-W_MR_STAR_Offset)*sizeof(C) );
            // Store a21[MR,* ]
            int APan_MR_STAR_Offset = 
                APan_MR_STAR.LocalHeight()-a21_MR_STAR.LocalHeight();
            std::memcpy
            ( APan_MR_STAR.LocalBuffer(APan_MR_STAR_Offset,A00.Width()),
              a21_MR_STAR.LocalBuffer(),
              (APan_MR_STAR.LocalHeight()-APan_MR_STAR_Offset)*sizeof(C) );

            // Update the portion of A22 that is in our current panel with 
            // w21Last and a21Last using two gers. We do not need their top 
            // entries. We trash the upper triangle of our panel of A since we 
            // are only doing slightly more work and we can replace it
            // afterwards.
            DistMatrix<C,MC,STAR> a21Last_MC_STAR_Bottom(g);
            DistMatrix<C,MR,STAR> a21Last_MR_STAR_Bottom(g);
            DistMatrix<C,MC,STAR> w21Last_MC_STAR_Bottom(g);
            DistMatrix<C,MR,STAR> w21Last_MR_STAR_Bottom(g);
            a21Last_MC_STAR_Bottom.View
            ( a21Last_MC_STAR, 1, 0, a21Last_MC_STAR.Height()-1, 1 );
            a21Last_MR_STAR_Bottom.View
            ( a21Last_MR_STAR, 1, 0, a21Last_MR_STAR.Height()-1, 1 );
            w21Last_MC_STAR_Bottom.View
            ( w21Last_MC_STAR, 1, 0, w21Last_MC_STAR.Height()-1, 1 );
            w21Last_MR_STAR_Bottom.View
            ( w21Last_MR_STAR, 1, 0, w21Last_MR_STAR.Height()-1, 1 );
            const C* a21_MC_STAR_Buffer = a21Last_MC_STAR_Bottom.LocalBuffer();
            const C* a21_MR_STAR_Buffer = a21Last_MR_STAR_Bottom.LocalBuffer();
            const C* w21_MC_STAR_Buffer = w21Last_MC_STAR_Bottom.LocalBuffer();
            const C* w21_MR_STAR_Buffer = w21Last_MR_STAR_Bottom.LocalBuffer();
            C* A22Buffer = A22.LocalBuffer();
            int localHeight = W22.LocalHeight();
            int localWidth = W22.LocalWidth();
            int lDim = A22.LocalLDim();
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
        PopBlocksizeStack();
        p21_MC_STAR.SetToZero();
        q21_MR_STAR.SetToZero();
        basic::internal::LocalHemvColAccumulateL
        ( (C)1, A22, a21_MC_STAR, a21_MR_STAR, p21_MC_STAR, q21_MR_STAR );
        PushBlocksizeStack( 1 );

        basic::Gemv
        ( ADJOINT, 
          (C)1, W20B.LockedLocalMatrix(),
                a21B_MC_STAR.LockedLocalMatrix(),
          (C)0, x01_MR_STAR.LocalMatrix() );
        basic::Gemv
        ( ADJOINT, 
          (C)1, A20B.LockedLocalMatrix(),
                a21B_MC_STAR.LockedLocalMatrix(),
          (C)0, y01_MR_STAR.LocalMatrix() );
        // Combine the AllReduce column summations of x01[MR,* ], y01[MR,* ],
        // and q21[MR,* ]
        {
            int x01LocalHeight = x01_MR_STAR.LocalHeight();
            int q21LocalHeight = q21_MR_STAR.LocalHeight();
            std::vector<C> colSumSendBuffer(2*x01LocalHeight+q21LocalHeight);
            std::vector<C> colSumRecvBuffer(2*x01LocalHeight+q21LocalHeight);
            std::memcpy
            ( &colSumSendBuffer[0], 
              x01_MR_STAR.LocalBuffer(), 
              x01LocalHeight*sizeof(C) );
            std::memcpy
            ( &colSumSendBuffer[x01LocalHeight],
              y01_MR_STAR.LocalBuffer(), 
              x01LocalHeight*sizeof(C) );
            std::memcpy
            ( &colSumSendBuffer[2*x01LocalHeight],
              q21_MR_STAR.LocalBuffer(), 
              q21LocalHeight*sizeof(C) );
            mpi::AllReduce
            ( &colSumSendBuffer[0], 
              &colSumRecvBuffer[0],
              2*x01LocalHeight+q21LocalHeight, mpi::SUM, g.MCComm() );
            std::memcpy
            ( x01_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[0], 
              x01LocalHeight*sizeof(C) );
            std::memcpy
            ( y01_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[x01LocalHeight], 
              x01LocalHeight*sizeof(C) );
            std::memcpy
            ( q21_MR_STAR.LocalBuffer(), 
              &colSumRecvBuffer[2*x01LocalHeight], 
              q21LocalHeight*sizeof(C) );
        }

        basic::Gemv
        ( NORMAL, 
          (C)-1, A20B.LockedLocalMatrix(),
                 x01_MR_STAR.LockedLocalMatrix(),
          (C)+1, p21B_MC_STAR.LocalMatrix() );
        basic::Gemv
        ( NORMAL, 
          (C)-1, W20B.LockedLocalMatrix(),
                 y01_MR_STAR.LockedLocalMatrix(),
          (C)+1, p21B_MC_STAR.LocalMatrix() );

        if( W22.Width() > 0 )
        {
            // This is not the last iteration of the panel factorization, 
            // combine the Reduce to one of p21[MC,* ] with the redistribution 
            // of q21[MR,* ] -> q21[MC,MR] to the next process column.
            int localHeight = p21_MC_STAR.LocalHeight();
            std::vector<C> reduceToOneSendBuffer(2*localHeight);
            std::vector<C> reduceToOneRecvBuffer(2*localHeight);

            // Pack p21[MC,* ]
            std::memcpy
            ( &reduceToOneSendBuffer[0], 
              p21_MC_STAR.LocalBuffer(),
              localHeight*sizeof(C) );

            // Fill in contributions to q21[MC,MR] from q21[MR,* ]
            bool contributing = 
                ( q21_MR_STAR.ColShift() % g.GCD() ==
                  p21_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    std::memcpy
                    ( &reduceToOneSendBuffer[localHeight],
                      q21_MR_STAR.LocalBuffer(), 
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
                    int sourcePeriod = g.LCM() / c;
                    int targetPeriod = g.LCM() / r;
                    int a0 = p21_MC_STAR.ColShift();
                    int b0 = q21_MR_STAR.ColShift();

                    int sourceStart = 0;
                    int f = (r+a0-b0) % r;
                    for( int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    int globalShift = b0+sourceStart*c;
                    int targetStart = (globalShift-a0)/r;
                    int localLength =
                        LocalLength(localHeight,targetStart,targetPeriod);
                    const C* q21_MR_STAR_LocalBuffer = 
                        q21_MR_STAR.LocalBuffer();
                    int offset = localHeight + targetStart;
                    for( int i=0; i<localLength; ++i )                        
                        reduceToOneSendBuffer[offset+i*targetPeriod] = 
                            q21_MR_STAR_LocalBuffer[sourceStart+i*sourcePeriod];
                }
            }
            else
            {
                std::memset
                ( &reduceToOneSendBuffer[localHeight], 0, 
                  localHeight*sizeof(C) );
            }

            int nextProcessRow = (alpha11.ColAlignment()+1) % r;
            int nextProcessCol = (alpha11.RowAlignment()+1) % c;
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

                // Finish computing w21. During its computation, ensure that 
                // every process has a copy of the first element of the w21.
                // We know a priori that the first element of a21 is one.
                const C* a21_MC_STAR_LocalBuffer = a21_MC_STAR.LocalBuffer();
                C myDotProduct = blas::Dot
                    ( localHeight, &reduceToOneRecvBuffer[0],   1, 
                                   &a21_MC_STAR_LocalBuffer[0], 1 );
                C sendBuffer[2];
                C recvBuffer[2];
                sendBuffer[0] = myDotProduct;
                sendBuffer[1] = ( g.MCRank()==nextProcessRow ? 
                                  reduceToOneRecvBuffer[0] : 0 );
                mpi::AllReduce
                ( sendBuffer, recvBuffer, 2, mpi::SUM, g.MCComm() );
                C dotProduct = recvBuffer[0];

                // Set up for the next iteration by filling in the values for:
                // - w21LastLocalBuffer
                // - w21LastFirstEntry
                C scale = static_cast<C>(0.5)*dotProduct*Conj(tau);
                for( int i=0; i<localHeight; ++i )
                    w21LastLocalBuffer[i] = tau*
                        ( reduceToOneRecvBuffer[i]-
                          scale*a21_MC_STAR_LocalBuffer[i] );
                w21LastFirstEntry = tau*( recvBuffer[1]-scale );
            }
        }
        else
        {
            // This is the last iteration, our last task is to finish forming
            // w21[MC,* ] and w21[MR,* ] so that we may place them into W[MC,* ]
            // and W[MR,* ]
            int localHeight = p21_MC_STAR.LocalHeight();
            std::vector<C> allReduceSendBuffer(2*localHeight);
            std::vector<C> allReduceRecvBuffer(2*localHeight);

            // Pack p21[MC,* ]
            std::memcpy
            ( &allReduceSendBuffer[0], 
              p21_MC_STAR.LocalBuffer(),
              localHeight*sizeof(C) );

            // Fill in contributions to q21[MC,* ] from q21[MR,* ]
            bool contributing = 
                ( q21_MR_STAR.ColShift() % g.GCD() ==
                  p21_MC_STAR.ColShift() % g.GCD() );
            if( contributing )
            {
                if( r == c )
                {
                    std::memcpy
                    ( &allReduceSendBuffer[localHeight],
                      q21_MR_STAR.LocalBuffer(), 
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
                    int sourcePeriod = g.LCM() / c;
                    int targetPeriod = g.LCM() / r;
                    int a0 = p21_MC_STAR.ColShift();
                    int b0 = q21_MR_STAR.ColShift();

                    int sourceStart = 0;
                    int f = (r+a0-b0) % r;
                    for( int s=0; s<sourcePeriod; ++s )
                    {
                        if( (s*c) % r == f )
                        {
                            sourceStart = s;
                            break;
                        }
                    }

                    int globalShift = b0+sourceStart*c;
                    int targetStart = (globalShift-a0)/r;
                    int localLength = 
                        LocalLength(localHeight,targetStart,targetPeriod);
                    const C* q21_MR_STAR_LocalBuffer = 
                        q21_MR_STAR.LocalBuffer();
                    int offset = localHeight + targetStart;
                    for( int i=0; i<localLength; ++i )
                        allReduceSendBuffer[offset+i*targetPeriod] = 
                            q21_MR_STAR_LocalBuffer[sourceStart+i*sourcePeriod];
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
 
            // Finish computing w21.
            const C* a21_MC_STAR_LocalBuffer = a21_MC_STAR.LocalBuffer();
            C myDotProduct = blas::Dot
                ( localHeight, &allReduceRecvBuffer[0], 1, 
                               a21_MC_STAR_LocalBuffer, 1 );
            C dotProduct;
            mpi::AllReduce
            ( &myDotProduct, &dotProduct, 1, mpi::SUM, g.MCComm() );

            // Grab views into W[MC,* ] and W[MR,* ]
            DistMatrix<C,MC,STAR> w21_MC_STAR(g);
            DistMatrix<C,MR,STAR> w21_MR_STAR(g);
            w21_MC_STAR.View
            ( W_MC_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );
            w21_MR_STAR.View
            ( W_MR_STAR, W00.Height()+1, W00.Width(), w21.Height(), 1 );

            // Store w21[MC,* ]
            C scale = static_cast<C>(0.5)*dotProduct*Conj(tau);
            C* w21_MC_STAR_LocalBuffer = w21_MC_STAR.LocalBuffer();
            for( int i=0; i<localHeight; ++i )
                w21_MC_STAR_LocalBuffer[i] = tau*
                    ( allReduceRecvBuffer[i]-
                      scale*a21_MC_STAR_LocalBuffer[i] );

            // Form w21[MR,* ]
            w21_MR_STAR = w21_MC_STAR;
        }
        //--------------------------------------------------------------------//
        a21_MC_STAR.FreeAlignments();
        a21_MR_STAR.FreeAlignments();
        p21_MC_STAR.FreeAlignments();
        p21_MR_STAR.FreeAlignments();
        q21_MC_STAR.FreeAlignments();
        q21_MR_STAR.FreeAlignments();
        x01_MR_STAR.FreeAlignments();
        y01_MR_STAR.FreeAlignments();

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
        firstIteration = false;
    }
    PopBlocksizeStack();

    // View the portion of A that e is the subdiagonal of, then place e into it
    DistMatrix<C,MC,MR> expandedATL(g);
    expandedATL.View( A, 0, 0, panelSize+1, panelSize+1 );
    expandedATL.SetRealDiagonal( e, -1 );
#ifndef RELEASE
    PopCallStack();
#endif
}
#endif // WITHOUT_COMPLEX

