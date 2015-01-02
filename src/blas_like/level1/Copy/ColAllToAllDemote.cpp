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

template<typename T,Dist U,Dist V>
void ColAllToAllDemote
( const DistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A,
        DistMatrix<T,        U,                     V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::ColAllToAllDemote"))
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignColsAndResize( A.ColAlign(), height, width, false, false );
    if( !B.Participating() )
        return;

    const Int colAlign = B.ColAlign();
    const Int rowAlignA = A.RowAlign();

    const Int colStride = B.ColStride();
    const Int colStridePart = B.PartialColStride();
    const Int colStrideUnion = B.PartialUnionColStride();
    const Int colRankPart = B.PartialColRank();
    const Int colDiff = (colAlign%colStridePart) - A.ColAlign();

    const Int colShiftA = A.ColShift();

    const Int localHeightB = B.LocalHeight();
    const Int localWidthA = A.LocalWidth();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int maxLocalWidth = MaxLength(width,colStrideUnion);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    std::vector<T> buffer( 2*colStrideUnion*portionSize );
    T* firstBuf  = &buffer[0];
    T* secondBuf = &buffer[colStrideUnion*portionSize];

    if( colDiff == 0 )
    {
        // Pack            
        util::PartialColStridedPack
        ( height, localWidthA,
          colAlign, colStride, 
          colStrideUnion, colStridePart, colRankPart,
          colShiftA,
          A.LockedBuffer(), A.LDim(),
          firstBuf,         portionSize );

        // Simultaneously Scatter in columns and Gather in rows
        mpi::AllToAll
        ( firstBuf,  portionSize,
          secondBuf, portionSize, B.PartialUnionColComm() );

        // Unpack
        util::RowStridedUnpack
        ( localHeightB, width,
          rowAlignA, colStrideUnion,
          secondBuf, portionSize,
          B.Buffer(), B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            std::cerr << "Unaligned ColAllToAllDemote" << std::endl;
#endif
        const Int sendColRankPart = Mod( colRankPart+colDiff, colStridePart );
        const Int recvColRankPart = Mod( colRankPart-colDiff, colStridePart );

        // Pack
        util::PartialColStridedPack
        ( height, localWidthA,
          colAlign, colStride, 
          colStrideUnion, colStridePart, sendColRankPart,
          colShiftA,
          A.LockedBuffer(), A.LDim(),
          secondBuf,        portionSize );

        // Simultaneously Scatter in columns and Gather in rows
        mpi::AllToAll
        ( secondBuf, portionSize,
          firstBuf,  portionSize, B.PartialUnionColComm() );

        // Realign the result
        mpi::SendRecv
        ( firstBuf,  colStrideUnion*portionSize, sendColRankPart,
          secondBuf, colStrideUnion*portionSize, recvColRankPart,
          B.PartialColComm() );

        // Unpack
        util::RowStridedUnpack
        ( localHeightB, width,
          rowAlignA, colStrideUnion,
          secondBuf, portionSize,
          B.Buffer(), B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void ColAllToAllDemote
( const BlockDistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A,
        BlockDistMatrix<T,        U,                     V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::ColAllToAllDemote"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void ColAllToAllDemote \
  ( const DistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A, \
          DistMatrix<T,        U,                     V   >& B ); \
  template void ColAllToAllDemote \
  ( const BlockDistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A, \
          BlockDistMatrix<T,        U,                     V   >& B );

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
