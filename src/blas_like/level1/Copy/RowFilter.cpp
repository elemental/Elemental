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

// (U,Collect(V)) |-> (U,V)
template<typename T>
void RowFilter
( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("copy::RowFilter");
        if( A.ColDist() != B.ColDist() ||
            A.RowDist() != Collect(B.RowDist()) )
            LogicError("Incompatible distributions");
    )
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

template<typename T>
void RowFilter
( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::RowFilter"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO(T) \
  template void RowFilter \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template void RowFilter \
  ( const AbstractBlockDistMatrix<T>& A, AbstractBlockDistMatrix<T>& B );

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
