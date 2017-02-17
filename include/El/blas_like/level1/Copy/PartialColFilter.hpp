/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_PARTIALCOLFILTER_HPP
#define EL_BLAS_COPY_PARTIALCOLFILTER_HPP

namespace El {
namespace copy {

template<typename T>
void PartialColFilter
( const ElementalMatrix<T>& A, ElementalMatrix<T>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.ColDist() != Partial(B.ColDist()) ||
          A.RowDist() != B.RowDist() )
          LogicError("Incompatible distributions");
    )
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignColsAndResize( A.ColAlign(), height, width, false, false );
    if( !B.Participating() )
        return;

    const Int colAlign = B.ColAlign();
    const Int colStride = B.ColStride();
    const Int colStridePart = B.PartialColStride();
    const Int colStrideUnion = B.PartialUnionColStride();
    const Int colShiftA = A.ColShift();
    const Int colDiff = Mod(colAlign,colStridePart)-A.ColAlign();

    const Int localHeight = B.LocalHeight();

    if( colDiff == 0 )
    {
        const Int colShift = B.ColShift();
        const Int colOffset = (colShift-colShiftA) / colStridePart;
        util::InterleaveMatrix
        ( localHeight, width,
          A.LockedBuffer(colOffset,0), colStrideUnion, A.LDim(),
          B.Buffer(),                  1,              B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            cerr << "Unaligned PartialColFilter" << endl;
#endif
        const Int colRankPart = B.PartialColRank();
        const Int colRankUnion = B.PartialUnionColRank();

        // Realign
        // -------
        const Int sendColRankPart = Mod( colRankPart+colDiff, colStridePart );
        const Int recvColRankPart = Mod( colRankPart-colDiff, colStridePart );
        const Int sendColRank = sendColRankPart + colStridePart*colRankUnion;
        const Int sendColShift = Shift( sendColRank, colAlign, colStride );
        const Int sendColOffset = (sendColShift-colShiftA) / colStridePart;
        const Int localHeightSend = Length( height, sendColShift, colStride );
        const Int sendSize = localHeightSend*width;
        const Int recvSize = localHeight    *width;
        vector<T> buffer;
        FastResize( buffer, sendSize+recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];
        // Pack
        util::InterleaveMatrix
        ( localHeightSend, width,
          A.LockedBuffer(sendColOffset,0), colStrideUnion, A.LDim(),
          sendBuf,                         1,              localHeightSend );
        // Change the column alignment
        mpi::SendRecv
        ( sendBuf, sendSize, sendColRankPart,
          recvBuf, recvSize, recvColRankPart, B.PartialColComm() );

        // Unpack
        // ------
        util::InterleaveMatrix
        ( localHeight, width,
          recvBuf,    1, localHeight,
          B.Buffer(), 1, B.LDim() );
    }
}

template<typename T>
void PartialColFilter
( const BlockMatrix<T>& A, BlockMatrix<T>& B )
{
    EL_DEBUG_CSE
    AssertSameGrids( A, B );
    // TODO(poulson): More efficient implementation
    GeneralPurpose( A, B );
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_PARTIALCOLFILTER_HPP
