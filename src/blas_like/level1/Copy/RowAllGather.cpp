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
void RowAllGather
( const DistMatrix<T,U,        V   >& A, 
        DistMatrix<T,U,Collect<V>()>& B ) 
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowAllGather"))
    AssertSameGrids( A, B );
    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignColsAndResize( A.ColAlign(), height, width, false, false );

    if( A.Participating() )
    {
        const Int colDiff = B.ColAlign() - A.ColAlign();
        if( colDiff == 0 )
        {
            if( width == 1 )
            {
                if( A.RowRank() == A.RowAlign() )
                    B.Matrix() = A.LockedMatrix();
                mpi::Broadcast
                ( B.Buffer(), B.LocalHeight(), A.RowAlign(), A.RowComm() );
            }
            else
            {
                const Int rowStride = A.RowStride();
                const Int localHeight = A.LocalHeight();
                const Int maxLocalWidth = MaxLength(width,rowStride);

                const Int portionSize = mpi::Pad( localHeight*maxLocalWidth );
                std::vector<T> buffer( (rowStride+1)*portionSize );
                T* sendBuf = &buffer[0];
                T* recvBuf = &buffer[portionSize];

                // Pack
                InterleaveMatrix
                ( localHeight, A.LocalWidth(),
                  A.LockedBuffer(), 1, A.LDim(),
                  sendBuf,          1, localHeight );

                // Communicate
                mpi::AllGather
                ( sendBuf, portionSize, recvBuf, portionSize, A.RowComm() );

                // Unpack
                const Int rowAlign = A.RowAlign();
                EL_OUTER_PARALLEL_FOR
                for( Int k=0; k<rowStride; ++k )
                {
                    const Int rowShift = Shift_( k, rowAlign, rowStride );
                    const Int localWidth =
                        Length_( width, rowShift, rowStride );
                    InterleaveMatrix
                    ( localHeight, localWidth,
                      &recvBuf[k*portionSize], 1, localHeight,
                      B.Buffer(0,rowShift),    1, rowStride*B.LDim() );
                }
            }
        }
        else
        {
#ifdef EL_UNALIGNED_WARNINGS
            if( A.Grid().Rank() == 0 )
                std::cerr << "Unaligned RowAllGather." << std::endl;
#endif
            const Int sendColRank = Mod( A.ColRank()+colDiff, A.ColStride() );
            const Int recvColRank = Mod( A.ColRank()-colDiff, A.ColStride() );

            if( width == 1 )
            {
                if( A.RowRank() == A.RowAlign() )
                    mpi::SendRecv
                    ( A.LockedBuffer(), A.LocalHeight(), sendColRank,
                      B.Buffer(),       B.LocalHeight(), recvColRank,
                      A.ColComm() );

                // Perform the row broadcast
                mpi::Broadcast
                ( B.Buffer(), B.LocalHeight(), A.RowAlign(), A.RowComm() );
            }
            else
            {
                const Int rowStride = A.RowStride();
                const Int localHeight = A.LocalHeight();
                const Int localWidthA = A.LocalWidth();
                const Int localHeightB = B.LocalHeight();
                const Int maxLocalHeight = MaxLength(height,A.ColStride());
                const Int maxLocalWidth = MaxLength(width,rowStride);

                const Int portionSize = mpi::Pad(maxLocalHeight*maxLocalWidth);
                std::vector<T> buffer( (rowStride+1)*portionSize );
                T* firstBuf = &buffer[0];
                T* secondBuf = &buffer[portionSize];

                // Pack
                InterleaveMatrix
                ( localHeight, localWidthA,
                  A.LockedBuffer(), 1, A.LDim(),
                  secondBuf,        1, localHeight );

                // Realign
                mpi::SendRecv
                ( secondBuf, portionSize, sendColRank,
                  firstBuf,  portionSize, recvColRank, A.ColComm() );

                // Perform the row AllGather
                mpi::AllGather
                ( firstBuf,  portionSize,
                  secondBuf, portionSize, A.RowComm() );

                // Unpack
                const Int rowAlign = A.RowAlign();
                EL_OUTER_PARALLEL_FOR
                for( Int k=0; k<rowStride; ++k )
                {
                    const Int rowShift = Shift_( k, rowAlign, rowStride );
                    const Int localWidth =
                        Length_( width, rowShift, rowStride );
                    InterleaveMatrix
                    ( localHeightB, localWidth,
                      &secondBuf[k*portionSize], 1, localHeightB,
                      B.Buffer(0,rowShift),      1, rowStride*B.LDim() );
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
  template void RowAllGather \
  ( const DistMatrix<T,U,        V   >& A, \
          DistMatrix<T,U,Collect<V>()>& B );

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
