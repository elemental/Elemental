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
void ColFilter
( const DistMatrix<T,Collect<U>(),V>& A,
        DistMatrix<T,        U,   V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::ColFilter"))
    AssertSameGrids( A, B );

    B.AlignRowsAndResize
    ( A.RowAlign(), A.Height(), A.Width(), false, false );
    if( !B.Participating() )
        return;

    const Int colStride = B.ColStride();
    const Int colShift = B.ColShift();
    const Int rowDiff = B.RowAlign() - A.RowAlign();

    const Int localHeight = B.LocalHeight();
    const Int localWidth = B.LocalWidth();

    if( rowDiff == 0 )
    {
        InterleaveMatrix
        ( localHeight, localWidth,
          A.LockedBuffer(colShift,0), colStride, A.LDim(),
          B.Buffer(),                 1,         B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            std::cerr << "Unaligned ColFilter" << std::endl;
#endif
        const Int rowStride = B.RowStride();
        const Int sendRowRank = Mod( B.RowRank()+rowDiff, rowStride );
        const Int recvRowRank = Mod( B.RowRank()-rowDiff, rowStride );
        const Int localWidthA = A.LocalWidth();
        const Int sendSize = localHeight*localWidthA;
        const Int recvSize = localHeight*localWidth;
        std::vector<T> buffer( sendSize+recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        InterleaveMatrix
        ( localHeight, localWidthA,
          A.LockedBuffer(colShift,0), colStride, A.LDim(),
          sendBuf,                    1,         localHeight );

        // Realign
        mpi::SendRecv
        ( sendBuf, sendSize, sendRowRank,
          recvBuf, recvSize, recvRowRank, B.RowComm() );

        // Unpack
        InterleaveMatrix
        ( localHeight, localWidth,
          recvBuf,    1, localHeight,
          B.Buffer(), 1, B.LDim() );
    }
}

#define PROTO_DIST(T,U,V) \
  template void ColFilter \
  ( const DistMatrix<T,Collect<U>(),V>& A, \
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

} // namespace copy
} // namespace El
