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
void SumScatter
( T alpha,
  const DistMatrix<T,Collect<U>(),Collect<V>()>& A,
        DistMatrix<T,        U,           V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::SumScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Sizes of A and B must match");
    if( !B.Participating() )
        return;

    const Int colStride = B.ColStride();
    const Int rowStride = B.RowStride();
    const Int colAlign = B.ColAlign();
    const Int rowAlign = B.RowAlign();

    const Int height = B.Height();
    const Int width = B.Width();
    const Int localHeight = B.LocalHeight();
    const Int localWidth = B.LocalWidth();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int maxLocalWidth = MaxLength(width,rowStride);

    const Int recvSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
    const Int sendSize = colStride*rowStride*recvSize;

    std::vector<T> buffer( sendSize );

    // Pack 
    copy::util::StridedPack
    ( height, width,
      colAlign, colStride,
      rowAlign, rowStride,
      A.LockedBuffer(), A.LDim(),
      buffer.data(),    recvSize );

    // Communicate
    mpi::ReduceScatter( buffer.data(), recvSize, B.DistComm() );

    // Unpack our received data
    util::InterleaveMatrixUpdate
    ( alpha, localHeight, localWidth,
      buffer.data(), 1, localHeight,
      B.Buffer(),    1, B.LDim() );
}

template<typename T,Dist U,Dist V>
void SumScatter
( T alpha,
  const BlockDistMatrix<T,Collect<U>(),Collect<V>()>& A,
        BlockDistMatrix<T,        U,           V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::SumScatter"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void SumScatter \
  ( T alpha, \
    const DistMatrix<T,Collect<U>(),Collect<V>()>& A, \
          DistMatrix<T,        U,           V   >& B ); \
  template void SumScatter \
  ( T alpha, \
    const BlockDistMatrix<T,Collect<U>(),Collect<V>()>& A, \
          BlockDistMatrix<T,        U,           V   >& B );

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
