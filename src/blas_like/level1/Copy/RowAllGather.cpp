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

// (U,V) |-> (U,Collect(V))
template<typename T>
void RowAllGather( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ) 
{
    DEBUG_ONLY(
        CallStackEntry cse("copy::RowAllGather");
        if( A.ColDist() != B.ColDist() || 
            Collect(A.RowDist()) != B.RowDist() )
            LogicError("Incompatible distributions");
    )
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
                util::InterleaveMatrix
                ( localHeight, A.LocalWidth(),
                  A.LockedBuffer(), 1, A.LDim(),
                  sendBuf,          1, localHeight );

                // Communicate
                mpi::AllGather
                ( sendBuf, portionSize, recvBuf, portionSize, A.RowComm() );

                // Unpack
                util::RowStridedUnpack
                ( localHeight, width, A.RowAlign(), rowStride,
                  recvBuf, portionSize,
                  B.Buffer(), B.LDim() );
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
                util::InterleaveMatrix
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
                util::RowStridedUnpack
                ( localHeightB, width, A.RowAlign(), rowStride,
                  secondBuf, portionSize,
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
void RowAllGather
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B ) 
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowAllGather"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO(T) \
  template void RowAllGather \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void RowAllGather \
  ( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
