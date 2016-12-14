/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_PARTIALROWFILTER_HPP
#define EL_BLAS_COPY_PARTIALROWFILTER_HPP

namespace El {
namespace copy {

template<typename T>
void PartialRowFilter
( const ElementalMatrix<T>& A, ElementalMatrix<T>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.ColDist() != B.ColDist() ||
          A.RowDist() != Partial(B.RowDist()) )
          LogicError("Incompatible distributions");
    )
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignRowsAndResize( A.RowAlign(), height, width, false, false );
    if( !B.Participating() )
        return;

    const Int rowAlign = B.RowAlign();
    const Int rowStride = B.RowStride();
    const Int rowStridePart = B.PartialRowStride();
    const Int rowStrideUnion = B.PartialUnionRowStride();
    const Int rowShiftA = A.RowShift();
    const Int rowDiff = (rowAlign%rowStridePart) - A.RowAlign();

    const Int localWidth = B.LocalWidth();

    if( rowDiff == 0 )
    {
        const Int rowShift = B.RowShift();
        const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
        util::InterleaveMatrix
        ( height, localWidth,
          A.LockedBuffer(0,rowOffset), 1, rowStrideUnion*A.LDim(),
          B.Buffer(),                  1, B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            cerr << "Unaligned PartialRowFilter" << endl;
#endif
        const Int rowRankPart = B.PartialRowRank();
        const Int rowRankUnion = B.PartialUnionRowRank();

        // Realign
        // -------
        const Int sendRowRankPart = Mod( rowRankPart+rowDiff, rowStridePart );
        const Int recvRowRankPart = Mod( rowRankPart-rowDiff, rowStridePart );
        const Int sendRowRank = sendRowRankPart + rowStridePart*rowRankUnion;
        const Int sendRowShift = Shift( sendRowRank, rowAlign, rowStride );
        const Int sendRowOffset = (sendRowShift-rowShiftA) / rowStridePart;
        const Int localWidthSend = Length( width, sendRowShift, rowStride );
        const Int sendSize = height*localWidthSend;
        const Int recvSize = height*localWidth;
        vector<T> buffer;
        FastResize( buffer, sendSize+recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];
        // Pack
        util::InterleaveMatrix
        ( height, localWidthSend,
          A.LockedBuffer(0,sendRowOffset), 1, rowStrideUnion*A.LDim(),
          sendBuf,                         1, height );
        // Change the column alignment
        mpi::SendRecv
        ( sendBuf, sendSize, sendRowRankPart,
          recvBuf, recvSize, recvRowRankPart, B.PartialRowComm() );

        // Unpack
        // ------
        util::InterleaveMatrix
        ( height, localWidth,
          recvBuf,    1, height,
          B.Buffer(), 1, B.LDim() );
    }
}

template<typename T>
void PartialRowFilter
( const BlockMatrix<T>& A, BlockMatrix<T>& B )
{
    EL_DEBUG_CSE
    AssertSameGrids( A, B );
    // TODO: More efficient implementation
    GeneralPurpose( A, B );
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_PARTIALROWFILTER_HPP
