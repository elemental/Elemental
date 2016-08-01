/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_ROWALLGATHER_HPP
#define EL_BLAS_COPY_ROWALLGATHER_HPP

namespace El {
namespace copy {

// (U,V) |-> (U,Collect(V))
template<typename T>
void RowAllGather( const ElementalMatrix<T>& A, ElementalMatrix<T>& B ) 
{
    DEBUG_CSE
    DEBUG_ONLY(
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
            if( A.RowStride() == 1 )
            {
                Copy( A.LockedMatrix(), B.Matrix() );
            }
            else if( width == 1 )
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
                vector<T> buffer;
                FastResize( buffer, (rowStride+1)*portionSize );
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
                cerr << "Unaligned RowAllGather." << endl;
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
                vector<T> buffer;
                FastResize( buffer, (rowStride+1)*portionSize );
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
        vector<T> buf;
        FastResize( buf, localHeight*localWidth );
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
( const BlockMatrix<T>& A, BlockMatrix<T>& B ) 
{
    DEBUG_CSE
    AssertSameGrids( A, B );
    // TODO: More efficient implementation
    GeneralPurpose( A, B );
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_ROWALLGATHER_HPP
