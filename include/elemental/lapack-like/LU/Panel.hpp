/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_LU_PANEL_HPP
#define LAPACK_LU_PANEL_HPP

#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level2/Geru.hpp"

namespace elem {
namespace lu {

template<typename F>
inline void
Panel( Matrix<F>& A, Matrix<int>& p, int pivotOffset=0 )
{
#ifndef RELEASE
    PushCallStack("lu::Panel");
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
    while( ATL.Width() < A.Width() )
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
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
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
Panel
( DistMatrix<F,  STAR,STAR>& A, 
  DistMatrix<F,  MC,  STAR>& B, 
  DistMatrix<int,STAR,STAR>& p, 
  int pivotOffset=0 )
{
#ifndef RELEASE
    PushCallStack("lu::Panel");
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
    // TODO: Think of how to make this safer with respect to alignment issues
    F* sendBufFloat = (F*)&sendData[0];
    F* recvBufFloat = (F*)&recvData[0];
    int* sendBufInt = (int*)&sendData[(width+1)*sizeof(F)];
    int* recvBufInt = (int*)&recvData[(width+1)*sizeof(F)];

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

            const int ALDim = A.LDim();
            const F* ABuffer = A.Buffer(pivotRow,0);
            for( int j=0; j<width; ++j )
                sendBufFloat[j+1] = ABuffer[j*ALDim];
        }
        else
        {
            const int localRow = ((pivotRow-A.Height())-colShift)/r;
            sendBufFloat[0] = b1.GetLocal(localRow,0);

            const int BLDim = B.LDim();
            const F* BBuffer = B.Buffer(localRow,0);
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
            const int ALDim = A.LDim();
            F* ASetBuffer = A.Buffer(pivotRow,0);
            const F* AGetBuffer = A.Buffer(currentRow,0);
            for( int j=0; j<width; ++j )
                ASetBuffer[j*ALDim] = AGetBuffer[j*ALDim];
        }
        else
        {
            const int ownerRank = (colAlignment+(pivotRow-A.Height())) % r;
            if( g.Row() == ownerRank )
            {
                const int localRow = ((pivotRow-A.Height())-colShift) / r;

                const int ALDim = A.LDim();
                const int BLDim = B.LDim();
                F* BBuffer = B.Buffer(localRow,0);
                const F* ABuffer = A.Buffer(currentRow,0);
                for( int j=0; j<width; ++j )
                    BBuffer[j*BLDim] = ABuffer[j*ALDim];
            }
        }

        // Copy the pivot row into the current row
        {
            F* ABuffer = A.Buffer(currentRow,0);
            const int ALDim = A.LDim();
            for( int j=0; j<width; ++j )
                ABuffer[j*ALDim] = recvBufFloat[j+1];
        }

        // Now we can perform the update of the current panel
        const F alpha = alpha11.GetLocal(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        Scale( alpha11Inv, a21.Matrix() );
        Scale( alpha11Inv, b1.Matrix()  );
        Geru( F(-1), a21.Matrix(), a12.Matrix(), A22.Matrix() );
        Geru( F(-1), b1.Matrix(), a12.Matrix(), B2.Matrix() );
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

} // namespace lu
} // namespace elem

#endif // ifndef LAPACK_LU_PANEL_HPP
