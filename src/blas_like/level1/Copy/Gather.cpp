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

template<typename T>
void Gather
( const AbstractDistMatrix<T>& A,
        DistMatrix<T,CIRC,CIRC>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::Gather"))
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.SetGrid( A.Grid() );
    B.Resize( height, width );

    // Gather the colShifts and rowShifts
    // ==================================
    Int myShifts[2];
    myShifts[0] = A.ColShift();
    myShifts[1] = A.RowShift();
    std::vector<Int> shifts;
    const Int crossSize = B.CrossSize();
    if( B.CrossRank() == B.Root() )
        shifts.resize( 2*crossSize );
    mpi::Gather( myShifts, 2, shifts.data(), 2, B.Root(), B.CrossComm() );

    // Gather the payload data
    // =======================
    const bool irrelevant = ( A.RedundantRank()!=0 || A.CrossRank()!=A.Root() );
    int totalSend = ( irrelevant ? 0 : A.LocalHeight()*A.LocalWidth() );
    std::vector<int> recvCounts, recvOffsets;
    if( B.CrossRank() == B.Root() )
        recvCounts.resize( crossSize );
    mpi::Gather( &totalSend, 1, recvCounts.data(), 1, B.Root(), B.CrossComm() );
    int totalRecv = Scan( recvCounts, recvOffsets );
    std::vector<T> sendBuf(totalSend), recvBuf(totalRecv);
    if( !irrelevant )
        copy::util::InterleaveMatrix
        ( A.LocalHeight(), A.LocalWidth(),
          A.LockedBuffer(), 1, A.LDim(),
          sendBuf.data(),   1, A.LocalHeight() );
    mpi::Gather
    ( sendBuf.data(), totalSend,
      recvBuf.data(), recvCounts.data(), recvOffsets.data(), 
      B.Root(), B.CrossComm() );

    // Unpack
    // ======
    if( B.Root() == B.CrossRank() )
    {
        for( Int q=0; q<crossSize; ++q )
        {
            if( recvCounts[q] == 0 )
                continue;
            const Int colShift = shifts[2*q+0];
            const Int rowShift = shifts[2*q+1];
            const Int colStride = A.ColStride();
            const Int rowStride = A.RowStride();
            const Int localHeight = Length( height, colShift, colStride );
            const Int localWidth = Length( width, rowShift, rowStride );
            copy::util::InterleaveMatrix
            ( localHeight, localWidth,
              &recvBuf[recvOffsets[q]],    1,         localHeight,
              B.Buffer(colShift,rowShift), colStride, rowStride*B.LDim() );
        }
    }
}

template<typename T>
void Gather
( const AbstractBlockDistMatrix<T>& A,
        BlockDistMatrix<T,CIRC,CIRC>& B )
{
    DEBUG_ONLY(CallStackEntry cse("copy::Gather"))
    AssertSameGrids( A, B );
    LogicError("This routine is not yet written");
}

#define PROTO(T) \
  template void Gather \
  ( const AbstractDistMatrix<T>& A, \
          DistMatrix<T,CIRC,CIRC>& B ); \
  template void Gather \
  ( const AbstractBlockDistMatrix<T>& A, \
          BlockDistMatrix<T,CIRC,CIRC>& B );

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
