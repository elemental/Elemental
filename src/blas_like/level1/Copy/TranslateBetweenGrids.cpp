/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace copy {

// TODO: Generalize to more distributions
template<typename T,Dist U,Dist V>
void TranslateBetweenGrids
( const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B ) 
{
    DEBUG_ONLY(CSE cse("copy::TranslateBetweenGrids"))
    LogicError("General TranslateBetweenGrids not yet supported");
}

template<typename T>
void TranslateBetweenGrids
( const DistMatrix<T,MC,MR>& A, DistMatrix<T,MC,MR>& B ) 
{
    DEBUG_ONLY(CSE cse("copy::TranslateBetweenGrids [MC,MR]"))

    B.Resize( A.Height(), A.Width() );
    // Just need to ensure that each viewing comm contains the other team's
    // owning comm. Congruence is too strong.

    // Compute the number of process rows and columns that each process
    // needs to send to.
    const Int colStride = B.ColStride();
    const Int rowStride = B.RowStride();
    const Int colRank = B.ColRank();
    const Int rowRank = B.RowRank();
    const Int colStrideA = A.ColStride();
    const Int rowStrideA = A.RowStride();
    const Int colGCD = GCD( colStride, colStrideA );
    const Int rowGCD = GCD( rowStride, rowStrideA );
    const Int colLCM = colStride*colStrideA / colGCD;
    const Int rowLCM = rowStride*rowStrideA / rowGCD;
    const Int numColSends = colStride / colGCD;
    const Int numRowSends = rowStride / rowGCD;

    const Int colAlign = B.ColAlign();
    const Int rowAlign = B.RowAlign();
    const Int colAlignA = A.ColAlign();
    const Int rowAlignA = A.RowAlign();

    const bool inBGrid = B.Grid().InGrid();
    const bool inAGrid = A.Grid().InGrid();
    if( !inBGrid && !inAGrid )
        return;

    const Int maxSendSize =
        (A.Height()/(colStrideA*numColSends)+1) *
        (A.Width()/(rowStrideA*numRowSends)+1);

    // Translate the ranks from A's VC communicator to B's viewing so that
    // we can match send/recv communicators. Since A's VC communicator is not
    // necessarily defined on every process, we instead work with A's owning
    // group and account for row-major ordering if necessary.
    const int sizeA = A.Grid().Size();
    vector<int> rankMap(sizeA), ranks(sizeA);
    if( A.Grid().Order() == COLUMN_MAJOR )
    {
        for( int j=0; j<sizeA; ++j )
            ranks[j] = j;
    }
    else
    {
        // The (i,j) = i + j*colStrideA rank in the column-major ordering is
        // equal to the j + i*rowStrideA rank in a row-major ordering.
        // Since we desire rankMap[i+j*colStrideA] to correspond to process
        // (i,j) in A's grid's rank in this viewing group, ranks[i+j*colStrideA]
        // should correspond to process (i,j) in A's owning group. Since the
        // owning group is ordered row-major in this case, its rank is
        // j+i*rowStrideA. Note that setting
        // ranks[j+i*rowStrideA] = i+j*colStrideA is *NOT* valid.
        for( int i=0; i<colStrideA; ++i )
            for( int j=0; j<rowStrideA; ++j )
                ranks[i+j*colStrideA] = j+i*rowStrideA;
    }
    mpi::Translate
    ( A.Grid().OwningGroup(), sizeA, &ranks[0],
      B.Grid().ViewingComm(), &rankMap[0] );

    // Have each member of A's grid individually send to all numRow x numCol
    // processes in order, while the members of this grid receive from all
    // necessary processes at each step.
    Int requiredMemory = 0;
    if( inAGrid )
        requiredMemory += maxSendSize;
    if( inBGrid )
        requiredMemory += maxSendSize;
    vector<T> auxBuf( requiredMemory );
    Int offset = 0;
    T* sendBuf = &auxBuf[offset];
    if( inAGrid )
        offset += maxSendSize;
    T* recvBuf = &auxBuf[offset];

    Int recvRow = 0; // avoid compiler warnings...
    if( inAGrid )
        recvRow = Mod(Mod(A.ColRank()-colAlignA,colStrideA)+colAlign,colStride);
    for( Int colSend=0; colSend<numColSends; ++colSend )
    {
        Int recvCol = 0; // avoid compiler warnings...
        if( inAGrid )
            recvCol=Mod(Mod(A.RowRank()-rowAlignA,rowStrideA)+rowAlign,
                        rowStride);
        for( Int rowSend=0; rowSend<numRowSends; ++rowSend )
        {
            mpi::Request sendRequest;
            // Fire off this round of non-blocking sends
            if( inAGrid )
            {
                // Pack the data
                Int sendHeight = Length(A.LocalHeight(),colSend,numColSends);
                Int sendWidth = Length(A.LocalWidth(),rowSend,numRowSends);
                copy::util::InterleaveMatrix
                ( sendHeight, sendWidth,
                  A.LockedBuffer(colSend,rowSend),
                 numColSends, numRowSends*A.LDim(),
                  sendBuf, 1, sendHeight );
                // Send data
                const Int recvVCRank = recvRow + recvCol*colStride;
                const Int recvViewingRank = B.Grid().VCToViewing( recvVCRank );
                mpi::ISend
                ( sendBuf, sendHeight*sendWidth, recvViewingRank,
                  B.Grid().ViewingComm(), sendRequest );
            }
            // Perform this round of recv's
            if( inBGrid )
            {
                const Int sendColOffset = colAlignA;
                const Int recvColOffset =
                    (colSend*colStrideA+colAlign) % colStride;
                const Int sendRowOffset = rowAlignA;
                const Int recvRowOffset =
                    (rowSend*rowStrideA+rowAlign) % rowStride;

                const Int firstSendRow =
                    Mod( Mod(colRank-recvColOffset,colStride)+sendColOffset,
                         colStrideA );
                const Int firstSendCol =
                    Mod( Mod(rowRank-recvRowOffset,rowStride)+sendRowOffset,
                         rowStrideA );

                const Int colShift = Mod( colRank-recvColOffset, colStride );
                const Int rowShift = Mod( rowRank-recvRowOffset, rowStride );
                const Int numColRecvs = Length( colStrideA, colShift, colStride );
                const Int numRowRecvs = Length( rowStrideA, rowShift, rowStride );

                // Recv data
                // For now, simply receive sequentially. Until we switch to
                // nonblocking recv's, we won't be using much of the
                // recvBuf
                Int sendRow = firstSendRow;
                for( Int colRecv=0; colRecv<numColRecvs; ++colRecv )
                {
                    const Int sendColShift = Shift( sendRow, colAlignA, colStrideA ) + colSend*colStrideA;
                    const Int sendHeight = Length( A.Height(), sendColShift, colLCM );
                    const Int localColOffset = (sendColShift-B.ColShift()) / colStride;

                    Int sendCol = firstSendCol;
                    for( Int rowRecv=0; rowRecv<numRowRecvs; ++rowRecv )
                    {
                        const Int sendRowShift = Shift( sendCol, rowAlignA, rowStrideA ) + rowSend*rowStrideA;
                        const Int sendWidth = Length( A.Width(), sendRowShift, rowLCM );
                        const Int localRowOffset = (sendRowShift-B.RowShift()) / rowStride;

                        const Int sendVCRank = sendRow+sendCol*colStrideA;
                        mpi::Recv
                        ( recvBuf, sendHeight*sendWidth, rankMap[sendVCRank],
                          B.Grid().ViewingComm() );

                        // Unpack the data
                        copy::util::InterleaveMatrix
                        ( sendHeight, sendWidth,
                          recvBuf, 1, sendHeight,
                          B.Buffer(localColOffset,localRowOffset),
                          colLCM/colStride, (rowLCM/rowStride)*B.LDim() );

                        // Set up the next send col
                        sendCol = (sendCol + rowStride) % rowStrideA;
                    }
                    // Set up the next send row
                    sendRow = (sendRow + colStride) % colStrideA;
                }
            }
            // Ensure that this round of non-blocking sends completes
            if( inAGrid )
            {
                mpi::Wait( sendRequest );
                recvCol = (recvCol + rowStrideA) % rowStride;
            }
        }
        if( inAGrid )
            recvRow = (recvRow + colStrideA) % colStride;
    }
}

template<typename T>
void TranslateBetweenGrids
( const DistMatrix<T,STAR,STAR>& A, DistMatrix<T,STAR,STAR>& B ) 
{
    DEBUG_ONLY(CSE cse("copy::TranslateBetweenGrids [STAR,STAR]"))
    const Int height = A.Height();
    const Int width = A.Width();
    B.Resize( height, width );

    // TODO:Decide whether this condition can be lifted or simplified.
    mpi::Comm viewingCommA = A.Grid().ViewingComm();
    mpi::Comm viewingCommB = B.Grid().ViewingComm();
    if( !mpi::Congruent( viewingCommA, viewingCommB ) )
        LogicError
        ("Redistributing between nonmatching grids currently requires"
         " the viewing communicators to match.");

    const Int rankA = A.RedundantRank();
    const Int rankB = B.RedundantRank();

    // Compute and allocate the amount of required memory
    Int requiredMemory = 0;
    if( rankA == 0 ) 
        requiredMemory += height*width;
    if( B.Participating() )
        requiredMemory += height*width;
    vector<T> buffer( requiredMemory );
    Int offset = 0;
    T* sendBuf = &buffer[offset];
    if( rankA == 0 ) 
        offset += height*width;
    T* bcastBuffer = &buffer[offset];

    // Send from the root of A to the root of B's matrix's grid
    mpi::Request sendRequest;
    if( rankA == 0 ) 
    {
        util::InterleaveMatrix
        ( height, width,
          A.LockedBuffer(), 1, A.LDim(),
          sendBuf,          1, height );
        // TODO: Use mpi::Translate instead?
        const Int recvViewingRank = B.Grid().VCToViewing(0);
        mpi::ISend
        ( sendBuf, height*width, recvViewingRank,
          viewingCommB, sendRequest );
    }

    // Receive on the root of B's matrix's grid and then broadcast
    // over the owning communicator
    if( B.Participating() )
    {
        if( rankB == 0 ) 
        {
            // TODO: Use mpi::Translate instead?
            const Int sendViewingRank = A.Grid().VCToViewing(0);
            mpi::Recv
            ( bcastBuffer, height*width, sendViewingRank,
              viewingCommB );
        }

        mpi::Broadcast( bcastBuffer, height*width, 0, B.RedundantComm() );

        util::InterleaveMatrix
        ( height, width,
          bcastBuffer, 1, height,
          B.Buffer(),  1, B.LDim() );
    }

    if( rankA == 0 )
        mpi::Wait( sendRequest );
}

#define PROTO_DIST(T,U,V) \
  template void TranslateBetweenGrids \
  ( const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B );

#define PROTO(T) \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,MR) \
  PROTO_DIST(T,MC,STAR) \
  PROTO_DIST(T,MD,STAR) \
  PROTO_DIST(T,MR,MC) \
  PROTO_DIST(T,MR,STAR) \
  PROTO_DIST(T,STAR,MC) \
  PROTO_DIST(T,STAR,MD) \
  PROTO_DIST(T,STAR,MR) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC) \
  PROTO_DIST(T,STAR,VR) \
  PROTO_DIST(T,VC,STAR) \
  PROTO_DIST(T,VR,STAR) 

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
