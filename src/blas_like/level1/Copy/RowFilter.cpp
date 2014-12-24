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
void RowFilter
( const DistMatrix<T,U,Collect<V>()>& A,
        DistMatrix<T,U,        V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowFilter"))
    AssertSameGrids( A, B );

    B.AlignColsAndResize
    ( A.ColAlign(), A.Height(), A.Width(), false, false );
    if( !B.Participating() )
        return;

    const Int colDiff = B.ColAlign() - A.ColAlign();
    const Int rowStride = B.RowStride();
    const Int rowShift = B.RowShift();

    const Int localHeight = B.LocalHeight();
    const Int localWidth = B.LocalWidth();

    if( colDiff == 0 )
    {
        util::InterleaveMatrix
        ( localHeight, localWidth,
          A.LockedBuffer(0,rowShift), 1, rowStride*A.LDim(),
          B.Buffer(),                 1, B.LDim() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            std::cerr << "Unaligned RowFilter" << std::endl;
#endif
        const Int colStride = B.ColStride();
        const Int sendColRank = Mod( B.ColRank()+colDiff, colStride );
        const Int recvColRank = Mod( B.ColRank()-colDiff, colStride );
        const Int localHeightA = A.LocalHeight();
        const Int sendSize = localHeightA*localWidth;
        const Int recvSize = localHeight *localWidth;

        std::vector<T> buffer( sendSize+recvSize );
        T* sendBuf = &buffer[0];
        T* recvBuf = &buffer[sendSize];

        // Pack
        util::InterleaveMatrix
        ( localHeightA, localWidth,
          A.LockedBuffer(0,rowShift), 1, rowStride*A.LDim(),
          sendBuf,                    1, localHeightA );

        // Realign
        mpi::SendRecv
        ( sendBuf, sendSize, sendColRank,
          recvBuf, recvSize, recvColRank, B.ColComm() );

        // Unpack
        util::InterleaveMatrix
        ( localHeight, localWidth,
          recvBuf,    1, localHeight,
          B.Buffer(), 1, B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void RowFilter
( const BlockDistMatrix<T,U,Collect<V>()>& A,
        BlockDistMatrix<T,U,        V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowFilter"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void RowFilter \
  ( const DistMatrix<T,U,Collect<V>()>& A, \
          DistMatrix<T,U,        V   >& B ); \
  template void RowFilter \
  ( const BlockDistMatrix<T,U,Collect<V>()>& A, \
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

} // namespace copy
} // namespace El
