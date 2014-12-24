/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace axpy {

template<typename T,Dist U,Dist V>
void PartialRowSumScatter
( T alpha,
  const DistMatrix<T,U,Partial<V>()>& A,
        DistMatrix<T,U,        V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::PartialRowSumScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrix sizes did not match");
    if( !B.Participating() )
        return;

    if( B.RowAlign() % A.RowStride() == A.RowAlign() )
    {
        const Int rowStride = B.RowStride();
        const Int rowStridePart = B.PartialRowStride();
        const Int rowStrideUnion = B.PartialUnionRowStride();
        const Int rowRankPart = B.PartialRowRank();

        const Int height = B.Height();
        const Int width = B.Width();
        const Int maxLocalWidth = MaxLength( width, rowStride );
        const Int recvSize = mpi::Pad( height*maxLocalWidth );
        const Int sendSize = rowStrideUnion*recvSize;

        std::vector<T> buffer( sendSize );

        // Pack
        copy::util::PartialRowStridedPack
        ( height, width,
          B.RowAlign(), rowStride, 
          rowStrideUnion, rowStridePart, rowRankPart,
          A.RowShift(),
          A.LockedBuffer(), A.LDim(),
          buffer.data(),    recvSize );

        // Communicate
        mpi::ReduceScatter( buffer.data(), recvSize, B.PartialUnionRowComm() );

        // Unpack our received data
        util::InterleaveMatrixUpdate
        ( alpha, height, B.LocalWidth(),
          buffer.data(), 1, height,
          B.Buffer(),    1, B.LDim() );
    }
    else
        LogicError("Unaligned axpy::PartialRowSumScatter not implemented");
}

template<typename T,Dist U,Dist V>
void PartialRowSumScatter
( T alpha,
  const BlockDistMatrix<T,U,Partial<V>()>& A,
        BlockDistMatrix<T,U,        V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::PartialRowSumScatter"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void PartialRowSumScatter \
  ( T alpha, \
    const DistMatrix<T,U,Partial<V>()>& A, \
          DistMatrix<T,U,        V   >& B ); \
  template void PartialRowSumScatter \
  ( T alpha, \
    const BlockDistMatrix<T,U,Partial<V>()>& A, \
          BlockDistMatrix<T,U,        V   >& B );

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
