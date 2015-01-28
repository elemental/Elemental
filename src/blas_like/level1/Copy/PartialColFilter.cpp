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

template<typename T>
void PartialColFilter
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("copy::PartialColFilter");
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
    const Int colDiff = (colAlign%colStridePart)-A.ColAlign();

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
            std::cerr << "Unaligned PartialColFilter" << std::endl;
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
        std::vector<T> buffer( sendSize+recvSize );
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
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::PartialColFilter"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO(T) \
  template void PartialColFilter \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void PartialColFilter \
  ( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
