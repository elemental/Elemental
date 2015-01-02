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

    const Int rowStride = B.RowStride();
    const Int rowStridePart = B.PartialRowStride();
    const Int rowStrideUnion = B.PartialUnionRowStride();
    const Int rowRankPart = B.PartialRowRank();
    const Int rowDiff = (rowAlign%rowStridePart) - A.RowAlign();

    const Int maxLocalHeight = MaxLength(height,rowStrideUnion);
    const Int maxLocalWidth = MaxLength(width,rowStride);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    std::vector<T> buffer( 2*rowStrideUnion*portionSize );
    T* firstBuf  = &buffer[0];
    T* secondBuf = &buffer[rowStrideUnion*portionSize];

    if( rowDiff == 0 )
    {
        // Pack            
        util::PartialRowStridedPack
        ( A.LocalHeight(), width,
          rowAlign, rowStride,
          rowStrideUnion, rowStridePart, rowRankPart,
          A.RowShift(),
          A.LockedBuffer(), A.LDim(),
          firstBuf,         portionSize );

        // Simultaneously Scatter in rows and Gather in columns
        mpi::AllToAll
        ( firstBuf,  portionSize,
          secondBuf, portionSize, B.PartialUnionRowComm() );

        // Unpack
        util::ColStridedUnpack
        ( height, B.LocalWidth(), 
          A.ColAlign(), rowStrideUnion,
          secondBuf, portionSize,
          B.Buffer(), B.LDim() );
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
        util::PartialRowStridedPack
        ( A.LocalHeight(), width,
          rowAlign, rowStride,
          rowStrideUnion, rowStridePart, sendRowRankPart,
          A.RowShift(),
          A.LockedBuffer(), A.LDim(),
          secondBuf,        portionSize );

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
        util::ColStridedUnpack
        ( height, B.LocalWidth(), 
          A.ColAlign(), rowStrideUnion,
          secondBuf,  portionSize,
          B.Buffer(), B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void RowAllToAllDemote
  ( const BlockDistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A, 
          BlockDistMatrix<T,                U,             V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowAllToAllDemote"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void RowAllToAllDemote \
  ( const DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A, \
          DistMatrix<T,                U,             V   >& B ); \
  template void RowAllToAllDemote \
  ( const BlockDistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A, \
          BlockDistMatrix<T,                U,             V   >& B );

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
