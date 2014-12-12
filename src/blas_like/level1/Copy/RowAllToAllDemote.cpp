/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void RowAllToAllDemote
  ( const DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A, 
          DistMatrix<T,                U,             V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowAllToAllDemote"))
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignRowsAndResize( A.RowAlign(), height, width, false, false );
    if( !B.Participating() )
        return;

    const Int rowAlign = B.RowAlign();
    const Int colAlignA = A.ColAlign();

    const Int rowStride = B.RowStride();
    const Int rowStridePart = B.PartialRowStride();
    const Int rowStrideUnion = B.PartialUnionRowStride();
    const Int rowRankPart = B.PartialRowRank();
    const Int rowDiff = (rowAlign%rowStridePart) - A.RowAlign();

    const Int rowShiftA = A.RowShift();

    const Int localWidthB = B.LocalWidth();
    const Int localHeightA = A.LocalHeight();
    const Int maxLocalHeight = MaxLength(height,rowStrideUnion);
    const Int maxLocalWidth = MaxLength(width,rowStride);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    std::vector<T> buffer( 2*rowStrideUnion*portionSize );
    T* firstBuf  = &buffer[0];
    T* secondBuf = &buffer[rowStrideUnion*portionSize];

    if( rowDiff == 0 )
    {
        // Pack            
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const Int rowRank = rowRankPart + k*rowStridePart;
            const Int rowShift = Shift_( rowRank, rowAlign, rowStride );
            const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
            const Int localWidth = Length_( width, rowShift, rowStride );
            InterleaveMatrix
            ( localHeightA, localWidth,
              A.LockedBuffer(0,rowOffset), 1, rowStrideUnion*A.LDim(),
              &firstBuf[k*portionSize],    1, localHeightA );
        }

        // Simultaneously Scatter in rows and Gather in columns
        mpi::AllToAll
        ( firstBuf,  portionSize,
          secondBuf, portionSize, B.PartialUnionRowComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const Int colShift = Shift_( k, colAlignA, rowStrideUnion );
            const Int localHeight = Length_( height, colShift, rowStrideUnion );
            InterleaveMatrix
            ( localHeight, localWidthB,
              &secondBuf[k*portionSize],  1,          localHeight,
              B.Buffer(colShift,0),   rowStrideUnion, B.LDim() );
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            std::cerr << "Unaligned RowAllToAllDemote" << std::endl;
#endif
        const Int sendRowRankPart = Mod( rowRankPart+rowDiff, rowStridePart );
        const Int recvRowRankPart = Mod( rowRankPart-rowDiff, rowStridePart );

        // Pack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const Int rowRank = sendRowRankPart + k*rowStridePart;
            const Int rowShift = Shift_( rowRank, rowAlign, rowStride );
            const Int rowOffset = (rowShift-rowShiftA) / rowStridePart;
            const Int localWidth = Length_( width, rowShift, rowStride );
            InterleaveMatrix
            ( localHeightA, localWidth,
              A.LockedBuffer(0,rowOffset), 1, rowStrideUnion*A.LDim(),
              &secondBuf[k*portionSize],   1, localHeightA );
        }

        // Simultaneously Scatter in rows and Gather in columns
        mpi::AllToAll
        ( secondBuf, portionSize,
          firstBuf,  portionSize, B.PartialUnionRowComm() );

        // Realign the result
        mpi::SendRecv
        ( firstBuf,  rowStrideUnion*portionSize, sendRowRankPart,
          secondBuf, rowStrideUnion*portionSize, recvRowRankPart,
          B.PartialRowComm() );

        // Unpack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<rowStrideUnion; ++k )
        {
            const Int colShift = Shift_( k, colAlignA, rowStrideUnion );
            const Int localHeight = Length_( height, colShift, rowStrideUnion );
            InterleaveMatrix
            ( localHeight, localWidthB,
              &secondBuf[k*portionSize], 1,          localHeight,
              B.Buffer(colShift,0),  rowStrideUnion, B.LDim() );
        }
    }
}

#define PROTO_DIST(T,U,V) \
  template void RowAllToAllDemote \
  ( const DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A, \
          DistMatrix<T,                U,             V   >& B );

#define PROTO(T) \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,  MR  ) \
  PROTO_DIST(T,MC,  STAR) \
  PROTO_DIST(T,MD,  STAR) \
  PROTO_DIST(T,MR,  MC  ) \
  PROTO_DIST(T,MR,  STAR) \
  PROTO_DIST(T,STAR,MC  ) \
  PROTO_DIST(T,STAR,MD  ) \
  PROTO_DIST(T,STAR,MR  ) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC  ) \
  PROTO_DIST(T,STAR,VR  ) \
  PROTO_DIST(T,VC,  STAR) \
  PROTO_DIST(T,VR,  STAR) 

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
