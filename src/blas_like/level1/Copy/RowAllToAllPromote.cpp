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
void RowAllToAllPromote
( const DistMatrix<T,                U,             V   >& A,
        DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowAllToAllPromote"))
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignRowsAndResize
    ( A.RowAlign()%B.RowStride(), height, width, false, false );
    if( !B.Participating() )
        return;

    const Int rowAlign = A.RowAlign();

    const Int rowStride = A.RowStride();
    const Int rowStridePart = A.PartialRowStride();
    const Int rowStrideUnion = A.PartialUnionRowStride();
    const Int rowRankPart = A.PartialRowRank();
    const Int rowDiff = B.RowAlign() - (rowAlign%rowStridePart);

    const Int maxLocalWidth = MaxLength(width,rowStride);
    const Int maxLocalHeight = MaxLength(height,rowStrideUnion);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    std::vector<T> buffer( 2*rowStrideUnion*portionSize );
    T* firstBuf  = &buffer[0];
    T* secondBuf = &buffer[rowStrideUnion*portionSize];

    if( rowDiff == 0 )
    {
        // Pack            
        util::ColStridedPack
        ( height, A.LocalWidth(),
          B.ColAlign(), rowStrideUnion,
          A.LockedBuffer(), A.LDim(),
          firstBuf,         portionSize );

        // Simultaneously Gather in rows and Scatter in columns
        mpi::AllToAll
        ( firstBuf,  portionSize,
          secondBuf, portionSize, A.PartialUnionRowComm() );

        // Unpack
        util::PartialRowStridedUnpack
        ( B.LocalHeight(), width,
          rowAlign, rowStride,
          rowStrideUnion, rowStridePart, rowRankPart,
          B.RowShift(),
          secondBuf, portionSize,
          B.Buffer(), B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( A.Grid().Rank() == 0 )
            std::cerr << "Unaligned RowAllToAllPromote" << std::endl;
#endif
        const Int sendRowRankPart = Mod( rowRankPart+rowDiff, rowStridePart );
        const Int recvRowRankPart = Mod( rowRankPart-rowDiff, rowStridePart );

        // Pack
        util::ColStridedPack
        ( height, A.LocalWidth(),
          B.ColAlign(), rowStrideUnion,
          A.LockedBuffer(), A.LDim(),
          secondBuf,        portionSize );

        // Realign the input
        mpi::SendRecv
        ( secondBuf, rowStrideUnion*portionSize, sendRowRankPart,
          firstBuf,  rowStrideUnion*portionSize, recvRowRankPart,
          A.PartialRowComm() );

        // Simultaneously Scatter in rows and Gather in columns
        mpi::AllToAll
        ( firstBuf,  portionSize,
          secondBuf, portionSize, A.PartialUnionRowComm() );

        // Unpack
        util::PartialRowStridedUnpack
        ( B.LocalHeight(), width,
          rowAlign, rowStride,
          rowStrideUnion, rowStridePart, recvRowRankPart,
          B.RowShift(),
          secondBuf, portionSize,
          B.Buffer(), B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void RowAllToAllPromote
( const BlockDistMatrix<T,                U,             V   >& A,
        BlockDistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowAllToAllPromote"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void RowAllToAllPromote \
  ( const DistMatrix<T,                U,             V   >& A, \
          DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& B ); \
  template void RowAllToAllPromote \
  ( const BlockDistMatrix<T,                U,             V   >& A, \
          BlockDistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& B );

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
