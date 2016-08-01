/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_COLALLTOALLDEMOTE_HPP
#define EL_BLAS_COPY_COLALLTOALLDEMOTE_HPP

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void ColAllToAllDemote
( const DistMatrix<T,Partial<U>(),PartialUnionRow<U,V>()>& A,
        DistMatrix<T,        U,                     V   >& B )
{
    DEBUG_CSE
    AssertSameGrids( A, B );

    const Int height = A.Height();
    const Int width = A.Width();
    B.AlignColsAndResize( A.ColAlign(), height, width, false, false );
    if( !B.Participating() )
        return;

    const Int colAlign = B.ColAlign();
    const Int rowAlignA = A.RowAlign();

    const Int colStride = B.ColStride();
    const Int colStridePart = B.PartialColStride();
    const Int colStrideUnion = B.PartialUnionColStride();
    const Int colRankPart = B.PartialColRank();
    const Int colDiff = Mod(colAlign,colStridePart) - A.ColAlign();

    const Int colShiftA = A.ColShift();

    const Int localHeightB = B.LocalHeight();
    const Int localWidthA = A.LocalWidth();
    const Int maxLocalHeight = MaxLength(height,colStride);
    const Int maxLocalWidth = MaxLength(width,colStrideUnion);
    const Int portionSize = mpi::Pad( maxLocalHeight*maxLocalWidth );

    if( colDiff == 0 )
    {
        if( B.PartialUnionColStride() == 1 )
        {
            Copy( A.LockedMatrix(), B.Matrix() );
        }
        else
        {
            vector<T> buffer;
            FastResize( buffer, 2*colStrideUnion*portionSize );
            T* firstBuf  = &buffer[0];
            T* secondBuf = &buffer[colStrideUnion*portionSize];

            // Pack            
            util::PartialColStridedPack
            ( height, localWidthA,
              colAlign, colStride, 
              colStrideUnion, colStridePart, colRankPart,
              colShiftA,
              A.LockedBuffer(), A.LDim(),
              firstBuf,         portionSize );

            // Simultaneously Scatter in columns and Gather in rows
            mpi::AllToAll
            ( firstBuf,  portionSize,
              secondBuf, portionSize, B.PartialUnionColComm() );

            // Unpack
            util::RowStridedUnpack
            ( localHeightB, width,
              rowAlignA, colStrideUnion,
              secondBuf, portionSize,
              B.Buffer(), B.LDim() );
        }
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( B.Grid().Rank() == 0 )
            cerr << "Unaligned ColAllToAllDemote" << endl;
#endif
        const Int sendColRankPart = Mod( colRankPart+colDiff, colStridePart );
        const Int recvColRankPart = Mod( colRankPart-colDiff, colStridePart );

        vector<T> buffer;
        FastResize( buffer, 2*colStrideUnion*portionSize );
        T* firstBuf  = &buffer[0];
        T* secondBuf = &buffer[colStrideUnion*portionSize];

        // Pack
        util::PartialColStridedPack
        ( height, localWidthA,
          colAlign, colStride, 
          colStrideUnion, colStridePart, sendColRankPart,
          colShiftA,
          A.LockedBuffer(), A.LDim(),
          secondBuf,        portionSize );

        // Simultaneously Scatter in columns and Gather in rows
        mpi::AllToAll
        ( secondBuf, portionSize,
          firstBuf,  portionSize, B.PartialUnionColComm() );

        // Realign the result
        mpi::SendRecv
        ( firstBuf,  colStrideUnion*portionSize, sendColRankPart,
          secondBuf, colStrideUnion*portionSize, recvColRankPart,
          B.PartialColComm() );

        // Unpack
        util::RowStridedUnpack
        ( localHeightB, width,
          rowAlignA, colStrideUnion,
          secondBuf, portionSize,
          B.Buffer(), B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void ColAllToAllDemote
( const DistMatrix<T,Partial<U>(),PartialUnionRow<U,V>(),BLOCK>& A,
        DistMatrix<T,        U,                     V   ,BLOCK>& B )
{
    DEBUG_CSE
    AssertSameGrids( A, B );
    // TODO: More efficient implementation
    GeneralPurpose( A, B );
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_COLALLTOALLDEMOTE_HPP
