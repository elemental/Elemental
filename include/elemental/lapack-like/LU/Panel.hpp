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
PanelLU( Matrix<F>& A, Matrix<int>& p, int pivotOffset )
{
#ifndef RELEASE
    PushCallStack("internal::PanelLU");
    if( A.Width() != p.Height() || p.Width() != 1 )
        throw std::logic_error("p must be a vector that conforms with A");
#endif
    // Matrix views
    Matrix<F> 
        ATL, ATR,  A00, a01,     A02,  
        ABL, ABR,  a10, alpha11, a12,  
                   A20, a21,     A22;

    const int width = A.Width();
    std::vector<F> buffer( width );

    // Start the algorithm
    PushBlocksizeStack( 1 );
    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22 );

        //--------------------------------------------------------------------//
        const int currentRow = A00.Height();
        
        // Find the index and value of the pivot candidate
        F pivot = alpha11.Get(0,0);
        int pivotRow = currentRow;
        for( int i=0; i<a21.Height(); ++i )
        {
            const F value = a21.Get(i,0);
            if( FastAbs(value) > FastAbs(pivot) )
            {
                pivot = value;
                pivotRow = currentRow + i + 1;
            }
        }
        p.Set( currentRow, 0, pivotRow+pivotOffset );

        // Swap the pivot row and current row
        for( int j=0; j<width; ++j )
        {
            buffer[j] = A.Get(currentRow,j);
            A.Set(currentRow,j,A.Get(pivotRow,j)); 
            A.Set(pivotRow,j,buffer[j]);
        }

        // Now we can perform the update of the current panel
        const F alpha = alpha11.Get(0,0);
        if( alpha == (F)0 )
            throw SingularMatrixException();
        const F alpha11Inv = ((F)1) / alpha;
        Scale( alpha11Inv, a21 );
        Geru( (F)-1, a21, a12, A22 );
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

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

        //--------------------------------------------------------------------//
        const int currentRow = a01.Height();
        
        // Store the index/value of the pivot candidate in A
        F pivot = alpha11.GetLocal(0,0);
        int pivotRow = currentRow;
        for( int i=0; i<a21.Height(); ++i )
        {
            F value = a21.GetLocal(i,0);
            if( FastAbs(value) > FastAbs(pivot) )
            {
                pivot = value;
                pivotRow = currentRow + i + 1;
            }
        }

        // Update the pivot candidate to include local data from B
        for( int i=0; i<B.LocalHeight(); ++i )
        {
            F value = b1.GetLocal(i,0);
            if( FastAbs(value) > FastAbs(pivot) )
            {
                pivot = value;
                pivotRow = A.Height() + colShift + i*r;
            }
        }

        // Fill the send buffer with:
        // [ pivotValue | pivot row data | pivotRow ]
        if( pivotRow < A.Height() )
        {
            sendBufFloat[0] = A.GetLocal(pivotRow,a10.Width());

            const int ALDim = A.LocalLDim();
            const F* ABuffer = A.LocalBuffer(pivotRow,0);
            for( int j=0; j<width; ++j )
                sendBufFloat[j+1] = ABuffer[j*ALDim];
        }
        else
        {
            const int localRow = ((pivotRow-A.Height())-colShift)/r;
            sendBufFloat[0] = b1.GetLocal(localRow,0);

            const int BLDim = B.LocalLDim();
            const F* BBuffer = B.LocalBuffer(localRow,0);
            for( int j=0; j<width; ++j )
                sendBufFloat[j+1] = BBuffer[j*BLDim];
        }
        *sendBufInt = pivotRow;

        // Communicate to establish the pivot information
        mpi::AllReduce
        ( &sendData[0], &recvData[0], numBytes, PivotOp<F>(), g.ColComm() );

        // Update the pivot vector
        pivotRow = *recvBufInt;
        p.SetLocal(currentRow,0,pivotRow+pivotOffset);

        // Copy the current row into the pivot row
        if( pivotRow < A.Height() )
        {
            const int ALDim = A.LocalLDim();
            F* ASetBuffer = A.LocalBuffer(pivotRow,0);
            const F* AGetBuffer = A.LocalBuffer(currentRow,0);
            for( int j=0; j<width; ++j )
                ASetBuffer[j*ALDim] = AGetBuffer[j*ALDim];
        }
        else
        {
            const int ownerRank = (colAlignment+(pivotRow-A.Height())) % r;
            if( g.Row() == ownerRank )
            {
                const int localRow = ((pivotRow-A.Height())-colShift) / r;

                const int ALDim = A.LocalLDim();
                const int BLDim = B.LocalLDim();
                F* BBuffer = B.LocalBuffer(localRow,0);
                const F* ABuffer = A.LocalBuffer(currentRow,0);
                for( int j=0; j<width; ++j )
                    BBuffer[j*BLDim] = ABuffer[j*ALDim];
            }
        }

        // Copy the pivot row into the current row
        {
            F* ABuffer = A.LocalBuffer(currentRow,0);
            const int ALDim = A.LocalLDim();
            for( int j=0; j<width; ++j )
                ABuffer[j*ALDim] = recvBufFloat[j+1];
        }

        // Now we can perform the update of the current panel
        const F alpha = alpha11.GetLocal(0,0);
        if( alpha == (F)0 )
            throw SingularMatrixException();
        const F alpha11Inv = ((F)1) / alpha;
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
    }
    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
