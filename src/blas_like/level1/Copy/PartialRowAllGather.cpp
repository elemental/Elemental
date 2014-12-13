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
void PartialRowAllGather
( const DistMatrix<T,U,        V   >& A, 
        DistMatrix<T,U,Partial<V>()>& B ) 
{
    DEBUG_ONLY(CallStackEntry cse("copy::PartialRowAllGather"))
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignRowsAndResize
    ( A.RowAlign()%B.RowStride(), height, width, false, false );
    if( !A.Participating() )
        return;

    DEBUG_ONLY(
        if( A.LocalHeight() != height )
            LogicError("This routine assumes columns are not distributed");
    )
    const Int rowStride = A.RowStride();
    const Int rowStrideUnion = A.PartialUnionRowStride();
    const Int rowStridePart = A.PartialRowStride();
    const Int rowRankPart = A.PartialRowRank();
    const Int rowDiff = B.RowAlign() - (A.RowAlign()%rowStridePart);

    const Int maxLocalWidth = MaxLength(width,rowStride);
    const Int portionSize = mpi::Pad( height*maxLocalWidth );
    std::vector<T> buffer( (rowStrideUnion+1)*portionSize );
    T* firstBuf = &buffer[0];
    T* secondBuf = &buffer[portionSize];

    if( rowDiff == 0 )
    {
        // Pack
        util::InterleaveMatrix
        ( height, A.LocalWidth(),
          A.LockedBuffer(), 1, A.LDim(),
          firstBuf,         1, height );

        // Communicate
        mpi::AllGather
        ( firstBuf, portionSize, secondBuf, portionSize,
          A.PartialUnionRowComm() );

        // Unpack
        util::PartialRowStridedUnpack
        ( height, width,
          A.RowAlign(), rowStride,
          rowStrideUnion, rowStridePart, rowRankPart,
          B.RowShift(),
          secondBuf, portionSize,
          B.Buffer(), B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( A.Grid().Rank() == 0 )
            std::cerr << "Unaligned PartialRowAllGather" << std::endl;
#endif
        // Perform a SendRecv to match the row alignments
        util::InterleaveMatrix
        ( height, A.LocalWidth(),
          A.LockedBuffer(), 1, A.LDim(),
          secondBuf,        1, height );
        const Int sendRowRank = Mod( A.RowRank()+rowDiff, rowStride );
        const Int recvRowRank = Mod( A.RowRank()-rowDiff, rowStride );
        mpi::SendRecv
        ( secondBuf, portionSize, sendRowRank,
          firstBuf,  portionSize, recvRowRank, A.RowComm() );

        // Use the SendRecv as an input to the partial union AllGather
        mpi::AllGather
        ( firstBuf,  portionSize,
          secondBuf, portionSize, A.PartialUnionRowComm() );

        // Unpack
        util::PartialRowStridedUnpack
        ( height, width,
          B.RowAlign(), rowStride,
          rowStrideUnion, rowStridePart, rowRankPart,
          B.RowShift(),
          secondBuf, portionSize,
          B.Buffer(), B.LDim() );
    }
}

#define PROTO_DIST(T,U,V) \
  template void PartialRowAllGather \
  ( const DistMatrix<T,U,        V   >& A, \
          DistMatrix<T,U,Partial<V>()>& B );

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
