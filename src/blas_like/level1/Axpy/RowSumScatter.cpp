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
void RowSumScatter
( T alpha,
  const DistMatrix<T,U,Collect<V>()>& A,
        DistMatrix<T,U,        V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::RowSumScatter"))
    AssertSameGrids( A, B );
    if( A.Height() != B.Height() || A.Width() != B.Width() )
        LogicError("Matrix sizes did not match");
    if( !B.Participating() )
        return;

    const Int width = B.Width();
    const Int colDiff = B.ColAlign()-A.ColAlign();
    if( colDiff == 0 )
    {
        if( width == 1 )
        {
            const Int localHeight = B.LocalHeight();
            const Int portionSize = mpi::Pad( localHeight );
            std::vector<T> buffer( portionSize );

            // Reduce to rowAlign
            const Int rowAlign = B.RowAlign();
            mpi::Reduce
            ( A.LockedBuffer(), buffer.data(), portionSize, 
              rowAlign, B.RowComm() );

            if( B.RowRank() == rowAlign )
                util::InterleaveMatrixUpdate
                ( alpha, localHeight, 1,
                  buffer.data(), 1, localHeight,
                  B.Buffer(),    1, B.LDim() );
        }
        else
        {
            const Int rowStride = B.RowStride();
            const Int rowAlign = B.RowAlign();

            const Int localHeight = B.LocalHeight();
            const Int localWidth = B.LocalWidth();
            const Int maxLocalWidth = MaxLength(width,rowStride);

            const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
            const Int sendSize = rowStride*portionSize;

            // Pack 
            std::vector<T> buffer( sendSize );
            copy::util::RowStridedPack
            ( localHeight, width,
              rowAlign, rowStride,
              A.LockedBuffer(), A.LDim(),
              buffer.data(), portionSize );

            // Communicate
            mpi::ReduceScatter( buffer.data(), portionSize, B.RowComm() );

            // Update with our received data
            util::InterleaveMatrixUpdate
            ( alpha, localHeight, localWidth,
              buffer.data(), 1, localHeight,
              B.Buffer(),    1, B.LDim() );
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            std::cerr << "Unaligned axpy::RowSumScatter" << std::endl;
#endif
        const Int colRank = B.ColRank();
        const Int colStride = B.ColStride();

        const Int sendRow = Mod( colRank+colDiff, colStride );
        const Int recvRow = Mod( colRank-colDiff, colStride );

        const Int localHeight = B.LocalHeight();
        const Int localHeightA = A.LocalHeight();

        if( width == 1 )
        {
            std::vector<T> buffer( localHeight+localHeightA );
            T* sendBuf = &buffer[0];
            T* recvBuf = &buffer[localHeightA];

            // Reduce to rowAlign
            const Int rowAlign = B.RowAlign();
            mpi::Reduce
            ( A.LockedBuffer(), sendBuf, localHeightA, rowAlign, B.RowComm() );

            if( B.RowRank() == rowAlign )
            {
                // Perform the realignment
                mpi::SendRecv
                ( sendBuf, localHeightA, sendRow,
                  recvBuf, localHeight,  recvRow, B.ColComm() );

                util::InterleaveMatrixUpdate
                ( alpha, localHeight, 1,
                  recvBuf,    1, localHeight,
                  B.Buffer(), 1, B.LDim() );
            }
        }
        else
        {
            const Int rowStride = B.RowStride();
            const Int rowAlign = B.RowAlign();

            const Int localWidth = B.LocalWidth();
            const Int maxLocalWidth = MaxLength(width,rowStride);

            const Int recvSize_RS = mpi::Pad( localHeightA*maxLocalWidth );
            const Int sendSize_RS = rowStride * recvSize_RS;
            const Int recvSize_SR = localHeight * localWidth;

            std::vector<T> buffer( recvSize_RS + Max(sendSize_RS,recvSize_SR) );
            T* firstBuf = &buffer[0];
            T* secondBuf = &buffer[recvSize_RS];

            // Pack 
            copy::util::RowStridedPack
            ( localHeightA, width,
              rowAlign, rowStride,
              A.LockedBuffer(), A.LDim(),
              secondBuf,        recvSize_RS );

            // Reduce-scatter over each process row
            mpi::ReduceScatter( secondBuf, firstBuf, recvSize_RS, B.RowComm() );

            // Trade reduced data with the appropriate process row
            mpi::SendRecv
            ( firstBuf,  localHeightA*localWidth, sendRow,
              secondBuf, localHeight*localWidth,  recvRow, B.ColComm() );

            // Update with our received data
            util::InterleaveMatrixUpdate
            ( alpha, localHeight, localWidth,
              secondBuf,  1, localHeight,
              B.Buffer(), 1, B.LDim() );
        }
    }    
}

template<typename T,Dist U,Dist V>
void RowSumScatter
( T alpha,
  const BlockDistMatrix<T,U,Collect<V>()>& A,
        BlockDistMatrix<T,U,        V   >& B )
{
    DEBUG_ONLY(CallStackEntry cse("axpy::RowSumScatter"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO_DIST(T,U,V) \
  template void RowSumScatter \
  ( T alpha, \
    const DistMatrix<T,U,Collect<V>()>& A, \
          DistMatrix<T,U,        V   >& B ); \
  template void RowSumScatter \
  ( T alpha, \
    const BlockDistMatrix<T,U,Collect<V>()>& A, \
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
