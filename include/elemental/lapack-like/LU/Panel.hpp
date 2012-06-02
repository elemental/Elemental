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

template<typename F>
inline void
PanelLU
( DistMatrix<F,  STAR,STAR>& A, 
  DistMatrix<F,  MC,  STAR>& B, 
  DistMatrix<int,STAR,STAR>& p, 
  int pivotOffset )
{
#ifndef RELEASE
    PushCallStack("internal::PanelLU");
    if( A.Grid() != p.Grid() || p.Grid() != B.Grid() )
        throw std::logic_error
        ("Matrices must be distributed over the same grid");
    if( A.Width() != B.Width() )
        throw std::logic_error("A and B must be the same width");
    if( A.Height() != p.Height() || p.Width() != 1 )
        throw std::logic_error("p must be a vector that conforms with A");
#endif
    const Grid& g = A.Grid();
    const int r = g.Height();
    const int colShift = B.ColShift();
    const int colAlignment = B.ColAlignment();

    // Matrix views
    DistMatrix<F,STAR,STAR> 
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),  
                         A20(g), a21(g),     A22(g);

    DistMatrix<F,MC,STAR>
        BL(g), BR(g),
        B0(g), b1(g), B2(g);

    DistMatrix<int,STAR,STAR>
        pT(g),  p0(g),
        pB(g),  psi1(g),
                p2(g);

    const int width = A.Width();
    const int numBytes = (width+1)*sizeof(F)+sizeof(int);
    std::vector<byte> sendData(numBytes);
    std::vector<byte> recvData(numBytes);

    // Extract pointers to send and recv data
    F* sendBufFloat = (F*) &sendData[0];
    F* recvBufFloat = (F*) &recvData[0];
    int* sendBufInt = (int*) &sendData[(width+1)*sizeof(F)];
    int* recvBufInt = (int*) &recvData[(width+1)*sizeof(F)];

    // Start the algorithm
    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    PartitionRight( B, BL, BR, 0 );
    PartitionDown
    ( p, pT,
         pB, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        RepartitionRight
        ( BL, /**/ BR,  
          B0, /**/ b1, B2 );

        RepartitionDown
        ( pT,  p0,
         /**/ /****/
               psi1,
          pB,  p2 );

        //--------------------------------------------------------------------//
        
        // Store the index/value of the pivot candidate in A
        F pivotValue = alpha11.GetLocalEntry(0,0);
        int pivotIndex = a01.Height();
        for( int i=0; i<a21.Height(); ++i )
        {
            F value = a21.GetLocalEntry(i,0);
            if( FastAbs(value) > FastAbs(pivotValue) )
            {
                pivotValue = value;
                pivotIndex = a01.Height() + i + 1;
            }
        }

        // Update the pivot candidate to include local data from B
        for( int i=0; i<B.LocalHeight(); ++i )
        {
            F value = b1.GetLocalEntry(i,0);
            if( FastAbs(value) > FastAbs(pivotValue) )
            {
                pivotValue = value;
                pivotIndex = A.Height() + colShift + i*r;
            }
        }

        // Fill the send buffer with:
        // [ pivotValue | pivotRow | pivotIndex ]
        if( pivotIndex < A.Height() )
        {
            sendBufFloat[0] = A.GetLocalEntry(pivotIndex,a10.Width());

            const int ALDim = A.LocalLDim();
            const F* ABuffer = A.LocalBuffer(pivotIndex,0);
            for( int j=0; j<width; ++j )
                sendBufFloat[j+1] = ABuffer[j*ALDim];
        }
        else
        {
            const int localIndex = ((pivotIndex-A.Height())-colShift)/r;
            sendBufFloat[0] = b1.GetLocalEntry(localIndex,0);

            const int BLDim = B.LocalLDim();
            const F* BBuffer = B.LocalBuffer(localIndex,0);
            for( int j=0; j<width; ++j )
                sendBufFloat[j+1] = BBuffer[j*BLDim];
        }
        *sendBufInt = pivotIndex;

        // Communicate to establish the pivot information
        mpi::AllReduce
        ( &sendData[0], &recvData[0], numBytes, PivotOp<F>(), g.ColComm() );

        // Update the pivot vector
        const int maxIndex = *recvBufInt;
        p.SetLocalEntry(a01.Height(),0,maxIndex+pivotOffset);

        // Copy the current row into the pivot row
        if( maxIndex < A.Height() )
        {
            const int ALDim = A.LocalLDim();
            F* ASetBuffer = A.LocalBuffer(maxIndex,0);
            const F* AGetBuffer = A.LocalBuffer(A00.Height(),0);
            for( int j=0; j<width; ++j )
                ASetBuffer[j*ALDim] = AGetBuffer[j*ALDim];
        }
        else
        {
            const int ownerRank = (colAlignment+(maxIndex-A.Height())) % r;
            if( g.Row() == ownerRank )
            {
                const int localIndex = ((maxIndex-A.Height())-colShift) / r;

                const int ALDim = A.LocalLDim();
                const int BLDim = B.LocalLDim();
                F* BBuffer = B.LocalBuffer(localIndex,0);
                const F* ABuffer = A.LocalBuffer(A00.Height(),0);
                for( int j=0; j<width; ++j )
                    BBuffer[j*BLDim] = ABuffer[j*ALDim];
            }
        }

        // Copy the pivot row into the current row
        {
            F* ABuffer = A.LocalBuffer(A00.Height(),0);
            const int ALDim = A.LocalLDim();
            for( int j=0; j<width; ++j )
                ABuffer[j*ALDim] = recvBufFloat[j+1];
        }

        // Now we can perform the update of the current panel
        F alpha = alpha11.GetLocalEntry(0,0);
        if( alpha == (F)0 )
            throw SingularMatrixException();
        F alpha11Inv = ((F)1) / alpha;
        Scale( alpha11Inv, a21.LocalMatrix() );
        Scale( alpha11Inv, b1.LocalMatrix()  );
        Geru( (F)-1, a21.LocalMatrix(), a12.LocalMatrix(), A22.LocalMatrix() );
        Geru( (F)-1, b1.LocalMatrix(), a12.LocalMatrix(), B2.LocalMatrix() );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );

        SlidePartitionRight
        ( BL,     /**/ BR,  
          B0, b1, /**/ B2 );

        SlidePartitionDown
        ( pT,  p0,
               psi1,
         /**/ /****/
          pB,  p2 );
    }
    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
