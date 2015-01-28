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

// (U,V) |-> (Collect(U),V)
template<typename T>
void ColAllGather( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("copy::ColAllGather");
        if( B.ColDist() != Collect(A.ColDist()) ||
            B.RowDist() != A.RowDist() )
            LogicError("Incompatible distributions");
    )
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
        const Int rowDiff = B.RowAlign()-A.RowAlign();
        if( rowDiff == 0 )
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
                util::InterleaveMatrix
                ( A.LocalHeight(), localWidth,
                  A.LockedBuffer(), 1, A.LDim(),
                  sendBuf,          1, A.LocalHeight() );

                // Communicate
                mpi::AllGather
                ( sendBuf, portionSize, recvBuf, portionSize, A.ColComm() );

                // Unpack
                util::ColStridedUnpack
                ( height, localWidth, A.ColAlign(), colStride,
                  recvBuf,    portionSize,
                  B.Buffer(), B.LDim() );
            }
        }
        else
        {
#ifdef EL_UNALIGNED_WARNINGS
            if( A.Grid().Rank() == 0 )
                std::cerr << "Unaligned [U,V] -> [* ,V]." << std::endl;
#endif
            const Int sendRowRank = Mod( A.RowRank()+rowDiff, A.RowStride() );
            const Int recvRowRank = Mod( A.RowRank()-rowDiff, A.RowStride() );

            if( height == 1 )
            {
                const Int localWidthB = B.LocalWidth();
                std::vector<T> buffer;
                T* bcastBuf;

                if( A.ColRank() == A.ColAlign() )
                {
                    const Int localWidth = A.LocalWidth();
                    buffer.resize( localWidth+localWidthB );
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
                const Int colStride = A.ColStride();
                const Int maxLocalHeight = MaxLength(height,colStride);
                const Int maxLocalWidth = MaxLength(width,A.RowStride());
                const Int portionSize =
                    mpi::Pad( maxLocalHeight*maxLocalWidth );

                std::vector<T> buffer( (colStride+1)*portionSize );
                T* firstBuf  = &buffer[0];
                T* secondBuf = &buffer[portionSize];

                // Pack
                util::InterleaveMatrix
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
                util::ColStridedUnpack
                ( height, B.LocalWidth(), A.ColAlign(), colStride,
                  secondBuf,  portionSize, 
                  B.Buffer(), B.LDim() );
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
            util::InterleaveMatrix
            ( localHeight, localWidth,
              B.LockedBuffer(), 1, B.LDim(),
              buf.data(),       1, localHeight );

        // Broadcast from the root
        mpi::Broadcast
        ( buf.data(), localHeight*localWidth, A.Root(), A.CrossComm() );

        // Unpack if not the root
        if( A.CrossRank() != A.Root() )
            util::InterleaveMatrix
            ( localHeight, localWidth,
              buf.data(), 1, localHeight,
              B.Buffer(), 1, B.LDim() );
    }
}

template<typename T>
void ColAllGather
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B ) 
{
    DEBUG_ONLY(CallStackEntry cse("copy::ColAllGather"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO(T) \
  template void ColAllGather \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void ColAllGather \
  ( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
