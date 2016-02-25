/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_GATHER_HPP
#define EL_BLAS_COPY_GATHER_HPP

namespace El {
namespace copy {

template<typename T>
void Gather
( const ElementalMatrix<T>& A,
        DistMatrix<T,CIRC,CIRC>& B )
{
    DEBUG_ONLY(CSE cse("copy::Gather"))
    AssertSameGrids( A, B );
    if( A.DistSize() == 1 && A.CrossSize() == 1 )
    {
        B.Resize( A.Height(), A.Width() );
        if( B.CrossRank() == B.Root() )
            Copy( A.LockedMatrix(), B.Matrix() );
        return;
    }

    const Int height = A.Height();
    const Int width = A.Width();
    B.SetGrid( A.Grid() );
    B.Resize( height, width );

    // Gather the colShifts and rowShifts
    // ==================================
    Int myShifts[2];
    myShifts[0] = A.ColShift();
    myShifts[1] = A.RowShift();
    vector<Int> shifts;
    const Int crossSize = B.CrossSize();
    if( B.CrossRank() == B.Root() )
        shifts.resize( 2*crossSize );
    mpi::Gather( myShifts, 2, shifts.data(), 2, B.Root(), B.CrossComm() );

    // Gather the payload data
    // =======================
    const bool irrelevant = ( A.RedundantRank()!=0 || A.CrossRank()!=A.Root() );
    int totalSend = ( irrelevant ? 0 : A.LocalHeight()*A.LocalWidth() );
    vector<int> recvCounts, recvOffsets;
    if( B.CrossRank() == B.Root() )
        recvCounts.resize( crossSize );
    mpi::Gather( &totalSend, 1, recvCounts.data(), 1, B.Root(), B.CrossComm() );
    int totalRecv = Scan( recvCounts, recvOffsets );
    vector<T> sendBuf, recvBuf;
    FastResize( sendBuf, totalSend );
    FastResize( recvBuf, totalRecv );
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
( const BlockMatrix<T>& A,
        DistMatrix<T,CIRC,CIRC,BLOCK>& B )
{
    DEBUG_ONLY(CSE cse("copy::Gather"))
    AssertSameGrids( A, B );
    if( A.DistSize() == 1 && A.CrossSize() == 1 )
    {
        B.Resize( A.Height(), A.Width() );
        if( B.CrossRank() == B.Root() )
            Copy( A.LockedMatrix(), B.Matrix() );
        return;
    }

    const Int height = A.Height();
    const Int width = A.Width();
    B.SetGrid( A.Grid() );
    B.Resize( height, width );

    // Gather the colShifts and rowShifts
    // ==================================
    Int myShifts[2];
    myShifts[0] = A.ColShift();
    myShifts[1] = A.RowShift();
    vector<Int> shifts;
    const Int crossSize = B.CrossSize();
    if( B.CrossRank() == B.Root() )
        shifts.resize( 2*crossSize );
    mpi::Gather( myShifts, 2, shifts.data(), 2, B.Root(), B.CrossComm() );

    // Gather the payload data
    // =======================
    const bool irrelevant = ( A.RedundantRank()!=0 || A.CrossRank()!=A.Root() );
    int totalSend = ( irrelevant ? 0 : A.LocalHeight()*A.LocalWidth() );
    vector<int> recvCounts, recvOffsets;
    if( B.CrossRank() == B.Root() )
        recvCounts.resize( crossSize );
    mpi::Gather( &totalSend, 1, recvCounts.data(), 1, B.Root(), B.CrossComm() );
    int totalRecv = Scan( recvCounts, recvOffsets );
    vector<T> sendBuf, recvBuf;
    FastResize( sendBuf, totalSend );
    FastResize( recvBuf, totalRecv );
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
    const Int mb = A.BlockHeight();
    const Int nb = A.BlockWidth();
    const Int colCut = A.ColCut();
    const Int rowCut = A.RowCut();
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
            const Int localHeight =
              BlockedLength( height, colShift, mb, colCut, colStride );
            const Int localWidth = 
              BlockedLength( width, rowShift, nb, rowCut, rowStride );
            const T* data = &recvBuf[recvOffsets[q]];
            for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            {
                const Int jBefore = rowShift*nb - rowCut;
                const Int jLocAdj = ( rowShift==0 ? jLoc+rowCut : jLoc );
                const Int numFilledLocalBlocks = jLocAdj / nb;
                const Int jMid = numFilledLocalBlocks*nb*rowStride;
                const Int jPost = jLocAdj-numFilledLocalBlocks*nb;
                const Int j = jBefore + jMid + jPost;
                const T* sourceCol = &data[jLoc*localHeight];
                for( Int iLoc=0; iLoc<localHeight; ++iLoc )
                {
                    const Int iBefore = colShift*mb - colCut;
                    const Int iLocAdj = (colShift==0 ? iLoc+colCut : iLoc);
                    const Int numFilledLocalBlocks = iLocAdj / mb;
                    const Int iMid = numFilledLocalBlocks*mb*colStride;
                    const Int iPost = iLocAdj-numFilledLocalBlocks*mb;
                    const Int i = iBefore + iMid + iPost;
                    B.SetLocal(i,j,sourceCol[iLoc]);
                }
            }
        }
    }
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_GATHER_HPP
