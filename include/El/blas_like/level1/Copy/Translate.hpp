/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_TRANSLATE_HPP
#define EL_BLAS_COPY_TRANSLATE_HPP

namespace El {
namespace copy {

template<typename T,Dist U,Dist V>
void Translate
( const DistMatrix<T,U,V>& A,
        DistMatrix<T,U,V>& B )
{
    EL_DEBUG_CSE
    if( A.Grid() != B.Grid() )
    {
        copy::TranslateBetweenGrids( A, B );
        return;
    }

    const Grid& g = A.Grid();
    const Int height = A.Height();
    const Int width = A.Width();
    const Int colAlign = A.ColAlign();
    const Int rowAlign = A.RowAlign();
    const Int root = A.Root();
    B.SetGrid( g );
    if( !B.RootConstrained() )
        B.SetRoot( root );
    if( !B.ColConstrained() )
        B.AlignCols( colAlign, false );
    if( !B.RowConstrained() )
        B.AlignRows( rowAlign, false );
    B.Resize( height, width );
    if( !g.InGrid() )
        return;

    const bool aligned = colAlign == B.ColAlign() && rowAlign == B.RowAlign();
    if( aligned && root == B.Root() )
    {
        Copy( A.LockedMatrix(), B.Matrix() );
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            cerr << "Unaligned [U,V] <- [U,V]" << endl;
#endif
        const Int colRank = A.ColRank();
        const Int rowRank = A.RowRank();
        const Int crossRank = A.CrossRank();
        const Int colStride = A.ColStride();
        const Int rowStride = A.RowStride();
        const Int maxHeight = MaxLength( height, colStride );
        const Int maxWidth  = MaxLength( width,  rowStride );
        const Int pkgSize = mpi::Pad( maxHeight*maxWidth );
        vector<T> buffer;
        if( crossRank == root || crossRank == B.Root() )
            FastResize( buffer, pkgSize );

        const Int colAlignB = B.ColAlign();
        const Int rowAlignB = B.RowAlign();
        const Int localHeightB =
            Length( height, colRank, colAlignB, colStride );
        const Int localWidthB = Length( width, rowRank, rowAlignB, rowStride );
        const Int recvSize = mpi::Pad( localHeightB*localWidthB );

        if( crossRank == root )
        {
            // Pack the local data
            util::InterleaveMatrix
            ( A.LocalHeight(), A.LocalWidth(),
              A.LockedBuffer(), 1, A.LDim(),
              buffer.data(),    1, A.LocalHeight() );

            if( !aligned )
            {
                // If we were not aligned, then SendRecv over the DistComm
                const Int colDiff = colAlignB-colAlign;
                const Int rowDiff = rowAlignB-rowAlign;

                const Int toRow = Mod(colRank+colDiff,colStride);
                const Int toCol = Mod(rowRank+rowDiff,rowStride);
                const Int toRank = toRow + toCol*colStride;

                const Int fromRow = Mod(colRank-colDiff,colStride);
                const Int fromCol = Mod(rowRank-rowDiff,rowStride);
                const Int fromRank = fromRow + fromCol*colStride;

                mpi::SendRecv
                ( buffer.data(), pkgSize, toRank, fromRank, A.DistComm() );
            }
        }
        if( root != B.Root() )
        {
            // Send to the correct new root over the cross communicator
            if( crossRank == root )
                mpi::Send( buffer.data(), recvSize, B.Root(), B.CrossComm() );
            else if( crossRank == B.Root() )
                mpi::Recv( buffer.data(), recvSize, root, B.CrossComm() );
        }
        // Unpack
        if( crossRank == B.Root() )
            util::InterleaveMatrix
            ( localHeightB, localWidthB,
              buffer.data(), 1, localHeightB,
              B.Buffer(),    1, B.LDim() );
    }
}

template<typename T,Dist U,Dist V>
void Translate
( const DistMatrix<T,U,V,BLOCK>& A,
        DistMatrix<T,U,V,BLOCK>& B )
{
    EL_DEBUG_CSE
    const Int height = A.Height();
    const Int width = A.Width();
    const Int blockHeight = A.BlockHeight();
    const Int blockWidth = A.BlockWidth();
    const Int colAlign = A.ColAlign();
    const Int rowAlign = A.RowAlign();
    const Int colCut = A.ColCut();
    const Int rowCut = A.RowCut();
    const Int root = A.Root();
    B.SetGrid( A.Grid() );
    if( !B.RootConstrained() )
        B.SetRoot( root, false );
    // TODO(poulson): Clarify under what conditions blocksizes are modified
    // (perhaps via a BlockSizeConstrained() function?)
    if( !B.ColConstrained() && B.BlockHeight() == blockHeight )
        B.AlignCols( blockHeight, colAlign, colCut, false );
    if( !B.RowConstrained() && B.BlockWidth() == blockWidth )
        B.AlignRows( blockWidth, rowAlign, rowCut, false );
    B.Resize( height, width );
    const bool aligned =
        blockHeight == B.BlockHeight() && blockWidth == B.BlockWidth() &&
        colAlign    == B.ColAlign()    && rowAlign   == B.RowAlign() &&
        colCut      == B.ColCut()      && rowCut     == B.RowCut();
    if( A.Grid().Size() == 1 || (aligned && root == B.Root()) )
    {
        Copy( A.LockedMatrix(), B.Matrix() );
    }
    else
    {
        GeneralPurpose( A, B );
    }
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_TRANSLATE_HPP
