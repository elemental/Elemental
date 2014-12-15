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

template<typename T,Dist U,Dist V>
void Translate( const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B ) 
{
    DEBUG_ONLY(CallStackEntry cse("copy::Translate"))
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
        B.Matrix() = A.LockedMatrix();
    }
    else
    {
#ifdef EL_UNALIGNED_WARNINGS
        if( g.Rank() == 0 )
            std::cerr << "Unaligned [U,V] <- [U,V]" << std::endl;
#endif
        const Int colRank = A.ColRank();
        const Int rowRank = A.RowRank();
        const Int crossRank = A.CrossRank();
        const Int colStride = A.ColStride();
        const Int rowStride = A.RowStride();
        const Int maxHeight = MaxLength( height, colStride );
        const Int maxWidth  = MaxLength( width,  rowStride );
        const Int pkgSize = mpi::Pad( maxHeight*maxWidth );
        std::vector<T> buffer;
        if( crossRank == root || crossRank == B.Root() )
            buffer.resize( pkgSize ); 

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

#define PROTO_DIST(T,U,V) \
  template void Translate( const DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B );

#define PROTO(T) \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,MR) \
  PROTO_DIST(T,MC,STAR) \
  PROTO_DIST(T,MD,STAR) \
  PROTO_DIST(T,MR,MC) \
  PROTO_DIST(T,MR,STAR) \
  PROTO_DIST(T,STAR,MC) \
  PROTO_DIST(T,STAR,MD) \
  PROTO_DIST(T,STAR,MR) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC) \
  PROTO_DIST(T,STAR,VR) \
  PROTO_DIST(T,VC,STAR) \
  PROTO_DIST(T,VR,STAR) 

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
