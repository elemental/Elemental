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
void PartialColSumScatter
( T alpha,
  const DistMatrix<T,Partial<U>(),V>& A,
        DistMatrix<T,        U,   V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::PartialColSumScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("A and B must be the same size");

#ifdef EL_CACHE_WARNINGS
    if( A.Width() != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "PartialColSumScatterUpdate potentially causes a large amount"
          " of cache-thrashing. If possible, avoid it by forming the "
          "(conjugate-)transpose of the [UGath,* ] matrix instead."
          << std::endl;
    }
#endif
    if( B.ColAlign() % A.ColStride() == A.ColAlign() )
    {
        const Int colStride = B.ColStride();
        const Int colStridePart = B.PartialColStride();
        const Int colStrideUnion = B.PartialUnionColStride();
        const Int colRankPart = B.PartialColRank();
        const Int colAlign = B.ColAlign();
        const Int colShiftOfA = A.ColShift();

        const Int height = B.Height();
        const Int width = B.Width();
        const Int localHeight = B.LocalHeight();
        const Int maxLocalHeight = MaxLength( height, colStride );
        const Int recvSize = mpi::Pad( maxLocalHeight*width );
        const Int sendSize = colStrideUnion*recvSize;

        std::vector<T> buffer( sendSize );

        // Pack
        EL_OUTER_PARALLEL_FOR
        for( Int k=0; k<colStrideUnion; ++k )
        {
            const Int thisRank = colRankPart+k*colStridePart;
            const Int thisColShift = Shift_( thisRank, colAlign, colStride );
            const Int thisColOffset =
                (thisColShift-colShiftOfA) / colStridePart;
            const Int thisLocalHeight =
                Length_( height, thisColShift, colStride );
            InterleaveMatrix
            ( thisLocalHeight, width,
              A.LockedBuffer(thisColOffset,0), colStrideUnion, A.LDim(),
              &buffer[k*recvSize], 1, thisLocalHeight );
        }

        // Communicate
        mpi::ReduceScatter( buffer.data(), recvSize, B.PartialUnionColComm() );

        // Unpack our received data
        InterleaveMatrixUpdate
        ( alpha, localHeight, width,
          buffer.data(), 1, localHeight,
          B.Buffer(),    1, B.LDim() );
    }
    else
        LogicError("Unaligned axpy::PartialColSumScatter not implemented");
}

#define PROTO_DIST(T,U,V) \
  template void PartialColSumScatter \
  ( T alpha, \
    const DistMatrix<T,Partial<U>(),V>& A, \
          DistMatrix<T,        U,   V>& B );

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
