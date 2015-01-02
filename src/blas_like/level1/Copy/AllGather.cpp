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
void AllGather
( const DistMatrix<T,        U,           V   >& A, 
        DistMatrix<T,Collect<U>(),Collect<V>()>& B ) 
{
    DEBUG_ONLY(CallStackEntry cse("copy::AllGather"))
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.SetGrid( A.Grid() );
    B.Resize( height, width );

    if( A.Participating() )
    {
        const Int colStride = A.ColStride();
        const Int rowStride = A.RowStride();
        const Int distStride = colStride*rowStride;
        const Int maxLocalHeight = MaxLength(height,colStride);
        const Int maxLocalWidth = MaxLength(width,rowStride);
        const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );
        std::vector<T> buf( (distStride+1)*portionSize );
        T* sendBuf = &buf[0];
        T* recvBuf = &buf[portionSize];

        // Pack
        util::InterleaveMatrix
        ( A.LocalHeight(), A.LocalWidth(),
          A.LockedBuffer(), 1, A.LDim(),
          sendBuf,          1, A.LocalHeight() );

        // Communicate
        mpi::AllGather
        ( sendBuf, portionSize, recvBuf, portionSize, A.DistComm() );

        // Unpack
        util::StridedUnpack
        ( height, width,
          A.ColAlign(), colStride,
          A.RowAlign(), rowStride,
          recvBuf, portionSize,
          B.Buffer(), B.LDim() );
    }
    if( A.Grid().InGrid() && A.CrossComm() != mpi::COMM_SELF )
    {
        // Pack from the root
        const Int BLocalHeight = B.LocalHeight();
        const Int BLocalWidth = B.LocalWidth();
        std::vector<T> buf(BLocalHeight*BLocalWidth);
        if( A.CrossRank() == A.Root() )
            util::InterleaveMatrix
            ( BLocalHeight, BLocalWidth,
              B.LockedBuffer(), 1, B.LDim(),
              buf.data(),       1, BLocalHeight ); 

        // Broadcast from the root
        mpi::Broadcast
        ( buf.data(), BLocalHeight*BLocalWidth, A.Root(), A.CrossComm() );

        // Unpack if not the root
        if( A.CrossRank() != A.Root() )
            util::InterleaveMatrix
            ( BLocalHeight, BLocalWidth,
              buf.data(), 1, BLocalHeight,
              B.Buffer(), 1, B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void AllGather
( const BlockDistMatrix<T,        U,           V   >& A, 
        BlockDistMatrix<T,Collect<U>(),Collect<V>()>& B ) 
{
    DEBUG_ONLY(CallStackEntry cse("copy::AllGather"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void AllGather \
  ( const DistMatrix<T,        U,           V   >& A, \
          DistMatrix<T,Collect<U>(),Collect<V>()>& B ); \
  template void AllGather \
  ( const BlockDistMatrix<T,        U,           V   >& A, \
          BlockDistMatrix<T,Collect<U>(),Collect<V>()>& B );

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
