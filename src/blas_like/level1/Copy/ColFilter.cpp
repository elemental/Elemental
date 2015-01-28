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

// (Collect(U),V) |-> (U,V)
template<typename T>
void ColFilter( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("copy::ColFilter");
        if( A.ColDist() != Collect(B.ColDist()) ||
            A.RowDist() != B.RowDist() )
            LogicError("Incompatible distributions");
    )
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
        util::InterleaveMatrix
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
        util::InterleaveMatrix
        ( localHeight, localWidthA,
          A.LockedBuffer(colShift,0), colStride, A.LDim(),
          sendBuf,                    1,         localHeight );

        // Realign
        mpi::SendRecv
        ( sendBuf, sendSize, sendRowRank,
          recvBuf, recvSize, recvRowRank, B.RowComm() );

        // Unpack
        util::InterleaveMatrix
        ( localHeight, localWidth,
          recvBuf,    1, localHeight,
          B.Buffer(), 1, B.LDim() );
    }
}

template<typename T>
void ColFilter
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::ColFilter"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO(T) \
  template void ColFilter \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void ColFilter \
  ( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
