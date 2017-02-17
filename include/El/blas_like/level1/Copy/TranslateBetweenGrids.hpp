/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_TRANSLATEBETWEENGRIDS_HPP
#define EL_BLAS_COPY_TRANSLATEBETWEENGRIDS_HPP

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void TranslateBetweenGrids
( const DistMatrix<T,U,V>& A,
        DistMatrix<T,U,V>& B )
{
    EL_DEBUG_CSE
    GeneralPurpose( A, B );
}

// TODO(poulson): Compare against copy::GeneralPurpose
template<typename T>
void TranslateBetweenGrids
( const DistMatrix<T,MC,MR>& A,
        DistMatrix<T,MC,MR>& B )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int mLocA = A.LocalHeight();
    const Int nLocA = A.LocalWidth();
    B.Resize( m, n );
    mpi::Comm viewingCommB = B.Grid().ViewingComm();
    mpi::Group owningGroupA = A.Grid().OwningGroup();

    // Just need to ensure that each viewing comm contains the other team's
    // owning comm. Congruence is too strong.

    // Compute the number of process rows and columns that each process
    // needs to send to.
    const Int colStride = B.ColStride();
    const Int rowStride = B.RowStride();
    const Int colShiftB = B.ColShift();
    const Int rowShiftB = B.RowShift();
    const Int colRank = B.ColRank();
    const Int rowRank = B.RowRank();
    const Int colRankA = A.ColRank();
    const Int rowRankA = A.RowRank();
    const Int colStrideA = A.ColStride();
    const Int rowStrideA = A.RowStride();
    const Int colGCD = GCD( colStride, colStrideA );
    const Int rowGCD = GCD( rowStride, rowStrideA );
    const Int colLCM = colStride*colStrideA / colGCD;
    const Int rowLCM = rowStride*rowStrideA / rowGCD;
    const Int numColSends = colStride / colGCD;
    const Int numRowSends = rowStride / rowGCD;

    const Int colAlignA = A.ColAlign();
    const Int rowAlignA = A.RowAlign();
    const Int colAlignB = B.ColAlign();
    const Int rowAlignB = B.RowAlign();

    const bool inBGrid = B.Participating();
    const bool inAGrid = A.Participating();
    if( !inBGrid && !inAGrid )
        return;

    const Int maxSendSize =
      (m/(colStrideA*numColSends)+1) * (n/(rowStrideA*numRowSends)+1);

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
    ( owningGroupA, sizeA, ranks.data(), viewingCommB, rankMap.data() );

    // Have each member of A's grid individually send to all numRow x numCol
    // processes in order, while the members of this grid receive from all
    // necessary processes at each step.
    Int requiredMemory = 0;
    if( inAGrid )
        requiredMemory += maxSendSize;
    if( inBGrid )
        requiredMemory += maxSendSize;
    vector<T> auxBuf;
    FastResize( auxBuf, requiredMemory );
    Int offset = 0;
    T* sendBuf = &auxBuf[offset];
    if( inAGrid )
        offset += maxSendSize;
    T* recvBuf = &auxBuf[offset];

    Int recvRow = 0; // avoid compiler warnings...
    if( inAGrid )
        recvRow = Mod(Mod(colRankA-colAlignA,colStrideA)+colAlignB,colStride);
    for( Int colSend=0; colSend<numColSends; ++colSend )
    {
        Int recvCol = 0; // avoid compiler warnings...
        if( inAGrid )
            recvCol=Mod(Mod(rowRankA-rowAlignA,rowStrideA)+rowAlignB,rowStride);
        for( Int rowSend=0; rowSend<numRowSends; ++rowSend )
        {
            mpi::Request<T> sendRequest;
            // Fire off this round of non-blocking sends
            if( inAGrid )
            {
                // Pack the data
                Int sendHeight = Length(mLocA,colSend,numColSends);
                Int sendWidth = Length(nLocA,rowSend,numRowSends);
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
                  viewingCommB, sendRequest );
            }
            // Perform this round of recv's
            if( inBGrid )
            {
                const Int sendColOffset = colAlignA;
                const Int recvColOffset =
                  Mod(colSend*colStrideA+colAlignB,colStride);
                const Int sendRowOffset = rowAlignA;
                const Int recvRowOffset =
                  Mod(rowSend*rowStrideA+rowAlignB,rowStride);

                const Int colShift = Mod( colRank-recvColOffset, colStride );
                const Int rowShift = Mod( rowRank-recvRowOffset, rowStride );

                const Int firstSendRow = Mod(colShift+sendColOffset,colStrideA);
                const Int firstSendCol = Mod(rowShift+sendRowOffset,rowStrideA);

                const Int numColRecvs = Length(colStrideA,colShift,colStride);
                const Int numRowRecvs = Length(rowStrideA,rowShift,rowStride);

                // Recv data
                // For now, simply receive sequentially. Until we switch to
                // nonblocking recv's, we won't be using much of the
                // recvBuf
                Int sendRow = firstSendRow;
                for( Int colRecv=0; colRecv<numColRecvs; ++colRecv )
                {
                    const Int sendColShift =
                      Shift( sendRow, colAlignA, colStrideA ) +
                      colSend*colStrideA;
                    const Int sendHeight = Length( m, sendColShift, colLCM );
                    const Int localColOffset =
                      (sendColShift-colShiftB) / colStride;

                    Int sendCol = firstSendCol;
                    for( Int rowRecv=0; rowRecv<numRowRecvs; ++rowRecv )
                    {
                        const Int sendRowShift =
                          Shift( sendCol, rowAlignA, rowStrideA ) +
                          rowSend*rowStrideA;
                        const Int sendWidth = Length( n, sendRowShift, rowLCM );
                        const Int localRowOffset =
                          (sendRowShift-rowShiftB) / rowStride;

                        const Int sendVCRank = sendRow+sendCol*colStrideA;
                        mpi::Recv
                        ( recvBuf, sendHeight*sendWidth, rankMap[sendVCRank],
                          viewingCommB );

                        // Unpack the data
                        copy::util::InterleaveMatrix
                        ( sendHeight, sendWidth,
                          recvBuf, 1, sendHeight,
                          B.Buffer(localColOffset,localRowOffset),
                          colLCM/colStride, (rowLCM/rowStride)*B.LDim() );

                        // Set up the next send col
                        sendCol = Mod(sendCol+rowStride,rowStrideA);
                    }
                    // Set up the next send row
                    sendRow = Mod(sendRow+colStride,colStrideA);
                }
            }
            // Ensure that this round of non-blocking sends completes
            if( inAGrid )
            {
                mpi::Wait( sendRequest );
                recvCol = Mod(recvCol+rowStrideA,rowStride);
            }
        }
        if( inAGrid )
            recvRow = Mod(recvRow+colStrideA,colStride);
    }
}

template<typename T>
void TranslateBetweenGrids
( const DistMatrix<T,STAR,STAR>& A,
        DistMatrix<T,STAR,STAR>& B )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    B.Resize( height, width );

    // Attempt to distinguish between the owning groups of A and B both being
    // subsets of the same viewing communicator, the owning group of A being
    // the same as the viewing communicator of B (A is the *parent* of B),
    // and the viewing communicator of A being the owning communicator of B
    // (B is the *parent* of A).
    //
    // TODO(poulson): Decide whether these condition can be simplified.
    mpi::Comm commA = A.Grid().VCComm();
    mpi::Comm commB = B.Grid().VCComm();
    mpi::Comm viewingCommA = A.Grid().ViewingComm();
    mpi::Comm viewingCommB = B.Grid().ViewingComm();
    const int commSizeA = mpi::Size(commA);
    const int commSizeB = mpi::Size(commB);
    const int viewingCommSizeA = mpi::Size(viewingCommA);
    const int viewingCommSizeB = mpi::Size(viewingCommB);
    bool usingViewingA=false, usingViewingB=false;
    mpi::Comm activeCommA, activeCommB;
    if( viewingCommSizeA == viewingCommSizeB )
    {
        usingViewingA = true;
        usingViewingB = true;
        activeCommA = viewingCommA;
        activeCommB = viewingCommB;
        if( !mpi::Congruent( viewingCommA, viewingCommB ) )
            LogicError("Viewing communicators were not congruent");
    }
    else if( viewingCommSizeA == commSizeB )
    {
        usingViewingA = true;
        usingViewingB = false;
        activeCommA = viewingCommA;
        activeCommB = commB;
        if( !mpi::Congruent( viewingCommA, commB ) )
            LogicError("Communicators were not congruent");
    }
    else if( commSizeA == viewingCommSizeB )
    {
        usingViewingA = false;
        usingViewingB = true;
        activeCommA = commA;
        activeCommB = viewingCommB;
        if( !mpi::Congruent( commA, viewingCommB ) )
            LogicError("Communicators were not congruent");
    }
    else
    {
        usingViewingA = false;
        usingViewingB = false;
        activeCommA = commA;
        activeCommB = commB;
        LogicError("Unsupported TranslateBetweenGrids instance");
    }

    const Int rankA = A.RedundantRank();
    const Int rankB = B.RedundantRank();

    // Compute and allocate the amount of required memory
    Int requiredMemory = 0;
    if( rankA == 0 )
        requiredMemory += height*width;
    if( B.Participating() )
        requiredMemory += height*width;
    vector<T> buffer;
    FastResize( buffer, requiredMemory );
    Int offset = 0;
    T* sendBuf = &buffer[offset];
    if( rankA == 0 )
        offset += height*width;
    T* bcastBuffer = &buffer[offset];

    // Send from the root of A to the root of B's matrix's grid
    mpi::Request<T> sendRequest;
    if( rankA == 0 )
    {
        util::InterleaveMatrix
        ( height, width,
          A.LockedBuffer(), 1, A.LDim(),
          sendBuf,          1, height );
        // TODO(poulson): Use mpi::Translate instead?
        const Int recvRank = ( usingViewingB ? B.Grid().VCToViewing(0) : 0 );
        mpi::ISend
        ( sendBuf, height*width, recvRank, activeCommB, sendRequest );
    }

    // Receive on the root of B's matrix's grid and then broadcast
    // over the owning communicator
    if( B.Participating() )
    {
        if( rankB == 0 )
        {
            // TODO(poulson): Use mpi::Translate instead?
            const Int sendRank =
              ( usingViewingA ? A.Grid().VCToViewing(0) : 0 );
            mpi::Recv( bcastBuffer, height*width, sendRank, activeCommB );
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

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_TRANSLATEBETWEENGRIDS_HPP
