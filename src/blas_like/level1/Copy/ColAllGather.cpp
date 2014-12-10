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
void ColAllGather
( const DistMatrix<T,        U,   V>& A, 
        DistMatrix<T,Collect<U>(),V>& B ) 
{
    DEBUG_ONLY(CallStackEntry cse("copy::ColAllGather"))
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
#ifdef EL_CACHE_WARNINGS
    if( height != 1 && A.Grid().Rank() == 0 )
    {
        std::cerr <<
          "The matrix redistribution [* ,V] <- [U,V] potentially causes a "
          "large amount of cache-thrashing. If possible, avoid it by "
          "performing the redistribution with a (conjugate)-transpose"
          << std::endl;
    }
#endif
    B.AlignRowsAndResize( A.RowAlign(), height, width, false, false );

    if( A.Participating() )
    {
        if( A.RowAlign() == B.RowAlign() )
        {
            if( height == 1 )
            {
                const Int localWidthB = B.LocalWidth();
                std::vector<T> bcastBuf(localWidthB);

                if( A.ColRank() == A.ColAlign() )
                {
                    B.Matrix() = A.LockedMatrix();
                    StridedMemCopy
                    ( bcastBuf.data(),  1,
                      B.LockedBuffer(), B.LDim(), localWidthB );
                }

                // Broadcast within the column comm
                mpi::Broadcast
                ( bcastBuf.data(), localWidthB, A.ColAlign(), A.ColComm() );

                // Unpack
                StridedMemCopy
                ( B.Buffer(),      B.LDim(), 
                  bcastBuf.data(), 1,        localWidthB );
            }
            else
            {
                const Int colStride = A.ColStride();
                const Int maxLocalHeight = MaxLength(height,colStride);
                const Int localWidth = A.LocalWidth();
                const Int portionSize = mpi::Pad( maxLocalHeight*localWidth );

                std::vector<T> buffer( (colStride+1)*portionSize );
                T* sendBuf = &buffer[0];
                T* recvBuf = &buffer[portionSize];

                // Pack
                InterleaveMatrix
                ( A.LocalHeight(), localWidth,
                  A.LockedBuffer(), 1, A.LDim(),
                  sendBuf,          1, A.LocalHeight() );

                // Communicate
                mpi::AllGather
                ( sendBuf, portionSize, recvBuf, portionSize, A.ColComm() );

                // Unpack
                const Int colAlign = A.ColAlign();
                EL_OUTER_PARALLEL_FOR
                for( Int k=0; k<colStride; ++k )
                {
                    const T* data = &recvBuf[k*portionSize];
                    const Int colShift = Shift_( k, colAlign, colStride );
                    const Int localHeight =
                        Length_( height, colShift, colStride );
                    InterleaveMatrix
                    ( localHeight, localWidth,
                      data,                 1,         localHeight,
                      B.Buffer(colShift,0), colStride, B.LDim() );
                }
            }
        }
        else
        {
#ifdef EL_UNALIGNED_WARNINGS
            if( A.Grid().Rank() == 0 )
                std::cerr << "Unaligned [U,V] -> [* ,V]." << std::endl;
#endif
            const Int rowDiff = B.RowAlign()-A.RowAlign();
            const Int rowStride = A.RowStride();
            const Int sendRowRank = Mod( A.RowRank()+rowDiff, rowStride );
            const Int recvRowRank = Mod( A.RowRank()-rowDiff, rowStride );

            if( height == 1 )
            {
                const Int localWidthB = B.LocalWidth();
                T* bcastBuf;

                if( A.ColRank() == A.ColAlign() )
                {
                    const Int localWidth = A.LocalWidth();

                    std::vector<T> buffer( localWidth+localWidthB );
                    T* sendBuf = &buffer[0];
                    bcastBuf   = &buffer[localWidth];

                    // Pack
                    StridedMemCopy
                    ( sendBuf, 1, A.LockedBuffer(), A.LDim(), localWidth );

                    // Communicate
                    mpi::SendRecv
                    ( sendBuf,  localWidth,  sendRowRank,
                      bcastBuf, localWidthB, recvRowRank, A.RowComm() );
                }
                else
                {
                    buffer.resize( localWidthB );
                    bcastBuf = buffer.data();
                }

                // Communicate
                mpi::Broadcast
                ( bcastBuf, localWidthB, A.ColAlign(), A.ColComm() );

                // Unpack
                StridedMemCopy
                ( B.Buffer(), B.LDim(), bcastBuf, 1, localWidthB );
            }
            else
            {
                const Int localWidthA = A.LocalWidth();
                const Int localWidthB = B.LocalWidth();
                const Int localHeightA = A.LocalHeight();
                const Int colStride = A.ColStride();
                const Int maxLocalHeight = MaxLength(height,colStride);
                const Int maxLocalWidth = MaxLength(width,rowStride);
                const Int portionSize =
                    mpi::Pad( maxLocalHeight*maxLocalWidth );

                std::vector<T> buffer( (colStride+1)*portionSize );
                T* firstBuf  = &buffer[0];
                T* secondBuf = &buffer[portionSize];

                // Pack
                InterleaveMatrix
                ( A.LocalHeight(), A.LocalWidth(),
                  A.LockedBuffer(), 1, A.LDim(),
                  secondBuf,        1, A.LocalHeight() );

                // Realign
                mpi::SendRecv
                ( secondBuf, portionSize, sendRowRank,
                  firstBuf,  portionSize, recvRowRank, A.RowComm() );

                // AllGather the aligned data
                mpi::AllGather
                ( firstBuf,  portionSize,
                  secondBuf, portionSize, A.ColComm() );

                // Unpack the contents of each member of the column team
                const Int colAlign = A.ColAlign();
                EL_OUTER_PARALLEL_FOR
                for( Int k=0; k<colStride; ++k )
                {
                    const T* data = &secondBuf[k*portionSize];
                    const Int colShift = Shift_( k, colAlign, colStride );
                    const Int localHeight =
                        Length_( height, colShift, colStride );
                    InterleaveMatrix
                    ( localHeight, localWidthB,
                      data,                 1,         localHeight,
                      B.Buffer(colShift,0), colStride, B.LDim() );
                }
            }
        }
    }
    if( A.Grid().InGrid() && A.CrossComm() != mpi::COMM_SELF )
    {
        // Pack from the root
        const Int localHeight = B.LocalHeight();
        const Int localWidth = B.LocalWidth();
        std::vector<T> buf( localHeight*localWidth );
        if( A.CrossRank() == A.Root() )
            InterleaveMatrix
            ( localHeight, localWidth,
              B.LockedBuffer(), 1, B.LDim(),
              buf.data(),       1, localHeight );

        // Broadcast from the root
        mpi::Broadcast
        ( buf.data(), localHeight*localWidth, A.Root(), A.CrossComm() );

        // Unpack if not the root
        if( A.CrossRank() != A.Root() )
            InterleaveMatrix
            ( localHeight, localWidth,
              buf.data(), 1, localHeight,
              B.Buffer(), 1, B.LDim() );
    }
}

#define PROTO_DIST(T,U,V) \
  template void ColAllGather \
  ( const DistMatrix<T,        U,   V>& A, \
          DistMatrix<T,Collect<U>(),V>& B );

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
