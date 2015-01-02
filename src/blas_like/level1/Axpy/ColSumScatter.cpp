/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace axpy {

template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha,
  const DistMatrix<T,Collect<U>(),V>& A,
        DistMatrix<T,        U,   V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::ColSumScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("A and B must be the same size");
#ifdef EL_VECTOR_WARNINGS
    if( A.Width() == 1 && B.Grid().Rank() == 0 )
    {
        std::cerr <<
          "The vector version of axpy::ColSumScatter does not"
          " yet have a vector version implemented, but it would only "
          "require a modification of the vector version of RowSumScatter"
          << std::endl;
    }
#endif
#ifdef EL_CACHE_WARNINGS
    if( A.Width() != 1 && B.Grid().Rank() == 0 )
    {
        std::cerr <<
          "axpy::ColSumScatter potentially causes a large "
          "amount of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [* ,V] matrix instead." << std::endl;
    }
#endif
    if( !B.Participating() )
        return;
    const Int height = B.Height();
    const Int localHeight = B.LocalHeight();
    const Int localWidth = B.LocalWidth();

    const Int colAlign = B.ColAlign();
    const Int colStride = B.ColStride();

    const Int rowDiff = B.RowAlign()-A.RowAlign();
    if( rowDiff == 0 )
    {
        const Int maxLocalHeight = MaxLength(height,colStride);

        const Int recvSize = mpi::Pad( maxLocalHeight*localWidth );
        const Int sendSize = colStride*recvSize;
        std::vector<T> buffer( sendSize );

        // Pack 
        copy::util::ColStridedPack
        ( height, localWidth,
          colAlign, colStride,
          A.LockedBuffer(), A.LDim(),
          buffer.data(),    recvSize );
    
        // Communicate
        mpi::ReduceScatter( buffer.data(), recvSize, B.ColComm() );

        // Update with our received data
        util::InterleaveMatrixUpdate
        ( alpha, localHeight, localWidth,
          buffer.data(), 1, localHeight,
          B.Buffer(),    1, B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            std::cerr << "Unaligned axpy::ColSumScatter" << std::endl;
#endif
        const Int localWidthA = A.LocalWidth();
        const Int maxLocalHeight = MaxLength(height,colStride);

        const Int recvSize_RS = mpi::Pad( maxLocalHeight*localWidthA );
        const Int sendSize_RS = colStride*recvSize_RS;
        const Int recvSize_SR = localHeight*localWidth;

        std::vector<T> buffer( recvSize_RS + Max(sendSize_RS,recvSize_SR) );
        T* firstBuf = &buffer[0];
        T* secondBuf = &buffer[recvSize_RS];

        // Pack
        copy::util::ColStridedPack
        ( height, localWidth,
          colAlign, colStride,
          A.LockedBuffer(), A.LDim(),
          secondBuf,        recvSize_RS );

        // Reduce-scatter over each col
        mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, B.ColComm() );

        // Trade reduced data with the appropriate col
        const Int sendCol = Mod( B.RowRank()+rowDiff, B.RowStride() );
        const Int recvCol = Mod( B.RowRank()-rowDiff, B.RowStride() );
        mpi::SendRecv
        ( firstBuf,  localHeight*localWidthA, sendCol,
          secondBuf, localHeight*localWidth,  recvCol, B.RowComm() );

        // Update with our received data
        util::InterleaveMatrixUpdate
        ( alpha, localHeight, localWidth,
          secondBuf,  1, localHeight,
          B.Buffer(), 1, B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void ColSumScatter
( T alpha,
  const BlockDistMatrix<T,Collect<U>(),V>& A,
        BlockDistMatrix<T,        U,   V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::ColSumScatter"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void ColSumScatter \
  ( T alpha, \
    const DistMatrix<T,Collect<U>(),V>& A, \
          DistMatrix<T,        U,   V>& B ); \
  template void ColSumScatter \
  ( T alpha, \
    const BlockDistMatrix<T,Collect<U>(),V>& A, \
          BlockDistMatrix<T,        U,   V>& B );

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

} // namespace axpy
} // namespace El
