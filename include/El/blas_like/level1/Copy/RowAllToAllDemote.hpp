/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_ROWALLTOALLDEMOTE_HPP
#define EL_BLAS_COPY_ROWALLTOALLDEMOTE_HPP

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void RowAllToAllDemote
  ( const DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>()>& A,
          DistMatrix<T,                U,             V   >& B )
{
    EL_DEBUG_CSE
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignRowsAndResize( A.RowAlign(), height, width, false, false );
    if( !B.Participating() )
        return;

    const Int rowAlign = B.RowAlign();

    const Int rowStride = B.RowStride();
    const Int rowStridePart = B.PartialRowStride();
    const Int rowStrideUnion = B.PartialUnionRowStride();
    const Int rowRankPart = B.PartialRowRank();
    const Int rowDiff = Mod(rowAlign,rowStridePart) - A.RowAlign();

    const Int maxLocalHeight = MaxLength(height,rowStrideUnion);
    const Int maxLocalWidth = MaxLength(width,rowStride);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    if( rowDiff == 0 )
    {
        if( B.PartialUnionRowStride() == 1 )
        {
            Copy( A.LockedMatrix(), B.Matrix() );
        }
        else
        {
            vector<T> buffer;
            FastResize( buffer, 2*rowStrideUnion*portionSize );
            T* firstBuf  = &buffer[0];
            T* secondBuf = &buffer[rowStrideUnion*portionSize];

            // Pack
            util::PartialRowStridedPack
            ( A.LocalHeight(), width,
              rowAlign, rowStride,
              rowStrideUnion, rowStridePart, rowRankPart,
              A.RowShift(),
              A.LockedBuffer(), A.LDim(),
              firstBuf,         portionSize );

            // Simultaneously Scatter in rows and Gather in columns
            mpi::AllToAll
            ( firstBuf,  portionSize,
              secondBuf, portionSize, B.PartialUnionRowComm() );

            // Unpack
            util::ColStridedUnpack
            ( height, B.LocalWidth(),
              A.ColAlign(), rowStrideUnion,
              secondBuf, portionSize,
              B.Buffer(), B.LDim() );
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            cerr << "Unaligned RowAllToAllDemote" << endl;
#endif
        const Int sendRowRankPart = Mod( rowRankPart+rowDiff, rowStridePart );
        const Int recvRowRankPart = Mod( rowRankPart-rowDiff, rowStridePart );

        vector<T> buffer;
        FastResize( buffer, 2*rowStrideUnion*portionSize );
        T* firstBuf  = &buffer[0];
        T* secondBuf = &buffer[rowStrideUnion*portionSize];

        // Pack
        util::PartialRowStridedPack
        ( A.LocalHeight(), width,
          rowAlign, rowStride,
          rowStrideUnion, rowStridePart, sendRowRankPart,
          A.RowShift(),
          A.LockedBuffer(), A.LDim(),
          secondBuf,        portionSize );

        // Simultaneously Scatter in rows and Gather in columns
        mpi::AllToAll
        ( secondBuf, portionSize,
          firstBuf,  portionSize, B.PartialUnionRowComm() );

        // Realign the result
        mpi::SendRecv
        ( firstBuf,  rowStrideUnion*portionSize, sendRowRankPart,
          secondBuf, rowStrideUnion*portionSize, recvRowRankPart,
          B.PartialRowComm() );

        // Unpack
        util::ColStridedUnpack
        ( height, B.LocalWidth(),
          A.ColAlign(), rowStrideUnion,
          secondBuf,  portionSize,
          B.Buffer(), B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void RowAllToAllDemote
  ( const DistMatrix<T,PartialUnionCol<U,V>(),Partial<V>(),BLOCK>& A,
          DistMatrix<T,                U,             V   ,BLOCK>& B )
{
    EL_DEBUG_CSE
    AssertSameGrids( A, B );
    // TODO(poulson): More efficient implementation
    GeneralPurpose( A, B );
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_ROWALLTOALLDEMOTE_HPP
