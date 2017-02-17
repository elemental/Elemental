/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_ROWALLGATHER_HPP
#define EL_BLAS_COPY_ROWALLGATHER_HPP

namespace El {
namespace copy {

// (U,V) |-> (U,Collect(V))
template<typename T>
void RowAllGather( const ElementalMatrix<T>& A, ElementalMatrix<T>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.ColDist() != B.ColDist() ||
          Collect(A.RowDist()) != B.RowDist() )
          LogicError("Incompatible distributions");
    )
    AssertSameGrids( A, B );
    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignColsAndResize( A.ColAlign(), height, width, false, false );

    if( A.Participating() )
    {
        const Int colDiff = B.ColAlign() - A.ColAlign();
        if( colDiff == 0 )
        {
            if( A.RowStride() == 1 )
            {
                Copy( A.LockedMatrix(), B.Matrix() );
            }
            else if( width == 1 )
            {
                if( A.RowRank() == A.RowAlign() )
                    B.Matrix() = A.LockedMatrix();
                mpi::Broadcast
                ( B.Buffer(), B.LocalHeight(), A.RowAlign(), A.RowComm() );
            }
            else
            {
                const Int rowStride = A.RowStride();
                const Int localHeight = A.LocalHeight();
                const Int maxLocalWidth = MaxLength(width,rowStride);

                const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
                vector<T> buffer;
                FastResize( buffer, (rowStride+1)*portionSize );
                T* sendBuf = &buffer[0];
                T* recvBuf = &buffer[portionSize];

                // Pack
                util::InterleaveMatrix
                ( localHeight, A.LocalWidth(),
                  A.LockedBuffer(), 1, A.LDim(),
                  sendBuf,          1, localHeight );

                // Communicate
                mpi::AllGather
                ( sendBuf, portionSize, recvBuf, portionSize, A.RowComm() );

                // Unpack
                util::RowStridedUnpack
                ( localHeight, width, A.RowAlign(), rowStride,
                  recvBuf, portionSize,
                  B.Buffer(), B.LDim() );
            }
        }
        else
        {
#ifdef EL_UNALIGNED_WARNINGS
            if( A.Grid().Rank() == 0 )
                Output("Unaligned RowAllGather");
#endif
            const Int sendColRank = Mod( A.ColRank()+colDiff, A.ColStride() );
            const Int recvColRank = Mod( A.ColRank()-colDiff, A.ColStride() );

            if( width == 1 )
            {
                if( A.RowRank() == A.RowAlign() )
                    mpi::SendRecv
                    ( A.LockedBuffer(), A.LocalHeight(), sendColRank,
                      B.Buffer(),       B.LocalHeight(), recvColRank,
                      A.ColComm() );

                // Perform the row broadcast
                mpi::Broadcast
                ( B.Buffer(), B.LocalHeight(), A.RowAlign(), A.RowComm() );
            }
            else
            {
                const Int rowStride = A.RowStride();
                const Int localHeight = A.LocalHeight();
                const Int localWidthA = A.LocalWidth();
                const Int localHeightB = B.LocalHeight();
                const Int maxLocalHeight = MaxLength(height,A.ColStride());
                const Int maxLocalWidth = MaxLength(width,rowStride);

                const Int portionSize = mpi::Pad(maxLocalHeight*maxLocalWidth);
                vector<T> buffer;
                FastResize( buffer, (rowStride+1)*portionSize );
                T* firstBuf = &buffer[0];
                T* secondBuf = &buffer[portionSize];

                // Pack
                util::InterleaveMatrix
                ( localHeight, localWidthA,
                  A.LockedBuffer(), 1, A.LDim(),
                  secondBuf,        1, localHeight );

                // Realign
                mpi::SendRecv
                ( secondBuf, portionSize, sendColRank,
                  firstBuf,  portionSize, recvColRank, A.ColComm() );

                // Perform the row AllGather
                mpi::AllGather
                ( firstBuf,  portionSize,
                  secondBuf, portionSize, A.RowComm() );

                // Unpack
                util::RowStridedUnpack
                ( localHeightB, width, A.RowAlign(), rowStride,
                  secondBuf, portionSize,
                  B.Buffer(), B.LDim() );
            }
        }
    }
    // Consider a A[STAR,MD] -> B[STAR,STAR] redistribution, where only the
    // owning team of the MD distribution of A participates in the initial phase
    // and the second phase broadcasts over the cross communicator.
    if( A.Grid().InGrid() && A.CrossComm() != mpi::COMM_SELF )
        El::Broadcast( B, A.CrossComm(), A.Root() );
}

template<typename T>
void RowAllGather( const BlockMatrix<T>& A, BlockMatrix<T>& B )
{
    EL_DEBUG_CSE
    AssertSameGrids( A, B );

    EL_DEBUG_ONLY(
      if( A.ColDist() != B.ColDist() ||
          Collect(A.RowDist()) != B.RowDist() )
          LogicError("Incompatible distributions");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    const Int colCut = A.ColCut();
    const Int rowCut = A.RowCut();
    const Int blockHeight = A.BlockHeight();
    const Int blockWidth = A.BlockWidth();

    B.AlignAndResize
    ( blockHeight, blockWidth, A.ColAlign(), 0, colCut, 0,
      height, width, false, false );
    // TODO(poulson): Realign if the cuts are different
    if( A.BlockHeight() != B.BlockHeight() || A.ColCut() != B.ColCut() )
    {
        EL_DEBUG_ONLY(
          Output("Performing expensive GeneralPurpose RowAllGather");
        )
        GeneralPurpose( A, B );
        return;
    }

    if( A.Participating() )
    {
        const Int colDiff = B.ColAlign() - A.ColAlign();
        const Int firstBlockWidth = blockWidth - rowCut;
        if( colDiff == 0 )
        {
            if( A.RowStride() == 1 )
            {
                Copy( A.LockedMatrix(), B.Matrix() );
            }
            else if( width <= firstBlockWidth )
            {
                if( A.RowRank() == A.RowAlign() )
                    B.Matrix() = A.LockedMatrix();
                El::Broadcast( B, A.RowComm(), A.RowAlign() );
            }
            else
            {
                const Int rowStride = A.RowStride();
                const Int localHeight = A.LocalHeight();
                const Int maxLocalWidth =
                  MaxBlockedLength(width,blockWidth,rowCut,rowStride);

                const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
                vector<T> buffer;
                FastResize( buffer, (rowStride+1)*portionSize );
                T* sendBuf = &buffer[0];
                T* recvBuf = &buffer[portionSize];

                // Pack
                util::InterleaveMatrix
                ( localHeight, A.LocalWidth(),
                  A.LockedBuffer(), 1, A.LDim(),
                  sendBuf,          1, localHeight );

                // Communicate
                mpi::AllGather
                ( sendBuf, portionSize, recvBuf, portionSize, A.RowComm() );

                // Unpack
                util::BlockedRowStridedUnpack
                ( localHeight, width, A.RowAlign(), rowStride,
                  A.BlockWidth(), A.RowCut(),
                  recvBuf, portionSize,
                  B.Buffer(), B.LDim() );
            }
        }
        else
        {
#ifdef EL_UNALIGNED_WARNINGS
            if( A.Grid().Rank() == 0 )
                Output("Unaligned RowAllGather");
#endif
            const Int sendColRank = Mod( A.ColRank()+colDiff, A.ColStride() );
            const Int recvColRank = Mod( A.ColRank()-colDiff, A.ColStride() );

            if( width <= firstBlockWidth )
            {
                if( A.RowRank() == A.RowAlign() )
                    El::SendRecv
                    ( A.LockedMatrix(), B.Matrix(),
                      A.ColComm(), sendColRank, recvColRank );
                El::Broadcast( B, A.RowComm(), A.RowAlign() );
            }
            else
            {
                const Int rowStride = A.RowStride();
                const Int localHeight = A.LocalHeight();
                const Int localWidthA = A.LocalWidth();
                const Int localHeightB = B.LocalHeight();
                const Int maxLocalHeight =
                  MaxBlockedLength(height,blockHeight,colCut,A.ColStride());
                const Int maxLocalWidth =
                  MaxBlockedLength(width,blockWidth,rowCut,rowStride);

                const Int portionSize = mpi::Pad(maxLocalHeight*maxLocalWidth);
                vector<T> buffer;
                FastResize( buffer, (rowStride+1)*portionSize );
                T* firstBuf = &buffer[0];
                T* secondBuf = &buffer[portionSize];

                // Pack
                util::InterleaveMatrix
                ( localHeight, localWidthA,
                  A.LockedBuffer(), 1, A.LDim(),
                  secondBuf,        1, localHeight );

                // Realign
                mpi::SendRecv
                ( secondBuf, portionSize, sendColRank,
                  firstBuf,  portionSize, recvColRank, A.ColComm() );

                // Perform the row AllGather
                mpi::AllGather
                ( firstBuf,  portionSize,
                  secondBuf, portionSize, A.RowComm() );

                // Unpack
                util::BlockedRowStridedUnpack
                ( localHeightB, width, A.RowAlign(), rowStride,
                  blockWidth, rowCut, secondBuf, portionSize,
                  B.Buffer(), B.LDim() );
            }
        }
    }
    // Consider a A[STAR,MD] -> B[STAR,STAR] redistribution, where only the
    // owning team of the MD distribution of A participates in the initial phase
    // and the second phase broadcasts over the cross communicator.
    if( A.Grid().InGrid() && A.CrossComm() != mpi::COMM_SELF )
        El::Broadcast( B, A.CrossComm(), A.Root() );
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_ROWALLGATHER_HPP
