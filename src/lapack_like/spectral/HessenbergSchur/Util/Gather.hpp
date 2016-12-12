/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_UTIL_GATHER_HPP
#define EL_HESS_SCHUR_UTIL_GATHER_HPP

namespace El {
namespace hess_schur {
namespace util {

template<typename Field>
void GatherSubdiagonal
( const DistMatrix<Field,MC,MR,BLOCK>& H,
  const IR& winInd,
        DistMatrix<Field,STAR,STAR>& hSubWin )
{
    EL_DEBUG_CSE
    const Int winSize = winInd.end - winInd.beg;
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();
    const auto& HLoc = H.LockedMatrix();
    EL_DEBUG_ONLY(
      if( H.BlockHeight() != H.BlockWidth() )
          LogicError("Assumed square distribution blocks");
      if( H.ColCut() != H.RowCut() )
          LogicError("Assumed symmetric cuts");
      if( blockSize < 2 )
          LogicError("Assumed blocks of size at least two");
    )
    Zeros( hSubWin, winSize-1, 1 );

    // TODO(poulson): Optimize these redistributions
    hSubWin.Reserve( winSize-1 );

    Int localRowOffset = H.LocalRowOffset( winInd.beg );
    Int localColOffset = H.LocalColOffset( winInd.beg );
    int rowOwner = H.RowOwner( winInd.beg );
    int colOwner = H.ColOwner( winInd.beg );

    Int thisCut = Mod( H.ColCut()+winInd.beg, H.BlockHeight() );
    Int winIndex = 0;
    while( winIndex < winSize )
    {
        const Int thisBlockSize = Min( blockSize-thisCut, winSize-winIndex );
        const int nextRowOwner = Mod( rowOwner+1, grid.Height() );
        const int nextColOwner = Mod( colOwner+1, grid.Width() );

        //  ---------------------------
        // | x x x x x x | x x x x x x |
        // | s x x x x x | x x x x x x |
        // |   s x x x x | x x x x x x |
        // |     s x x x | x x x x x x |
        // |       s x x | x x x x x x |
        // |         s x | x x x x x x |
        // |-------------|-------------|
        // |           s | x x x x x x |
        // |             | s x x x x x |
        // |             |   s x x x x |
        // |             |     s x x x |
        // |             |       s x x |
        // |             |         s x |
        //  ---------------------------
        const bool haveLastOff = (winSize > winIndex+thisBlockSize);

        // Handle this main diagonal block
        if( colOwner == grid.Col() && rowOwner == grid.Row() )
        {
            for( Int offset=0; offset<thisBlockSize; ++offset )
            {
                if( offset < thisBlockSize-1 )
                    hSubWin.QueueUpdate
                    ( winIndex+offset, 0,
                      HLoc(localRowOffset+offset+1,localColOffset+offset) );
            }
        }
        // Handle the last subdiagonal (if it exists)
        if( haveLastOff &&
            nextRowOwner == grid.Row() && colOwner == grid.Col() )
        {
            const Int subLocalRowOffset = localRowOffset +
              ( grid.Height() == 1 ? thisBlockSize : 0 );
            hSubWin.QueueUpdate
            ( winIndex+thisBlockSize-1, 0,
              HLoc(subLocalRowOffset,localColOffset+thisBlockSize-1) );
        }

        winIndex += thisBlockSize;
        thisCut = 0;
        if( colOwner == grid.Col() )
            localColOffset += thisBlockSize;
        if( rowOwner == grid.Row() )
            localRowOffset += thisBlockSize;
        rowOwner = nextRowOwner;
        colOwner = nextColOwner;
    }
    hSubWin.ProcessQueues();
}

template<typename Field>
void GatherBidiagonal
( const DistMatrix<Field,MC,MR,BLOCK>& H,
  const IR& winInd,
        DistMatrix<Field,STAR,STAR>& hMainWin,
        DistMatrix<Field,STAR,STAR>& hSubWin )
{
    EL_DEBUG_CSE
    const Int winSize = winInd.end - winInd.beg;
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();
    const auto& HLoc = H.LockedMatrix();
    EL_DEBUG_ONLY(
      if( H.BlockHeight() != H.BlockWidth() )
          LogicError("Assumed square distribution blocks");
      if( H.ColCut() != H.RowCut() )
          LogicError("Assumed symmetric cuts");
      if( blockSize < 2 )
          LogicError("Assumed blocks of size at least two");
    )
    Zeros( hMainWin, winSize, 1 );
    Zeros( hSubWin, winSize-1, 1 );

    // TODO(poulson): Optimize these redistributions
    hMainWin.Reserve( winSize );
    hSubWin.Reserve( winSize-1 );

    Int localRowOffset = H.LocalRowOffset( winInd.beg );
    Int localColOffset = H.LocalColOffset( winInd.beg );
    int rowOwner = H.RowOwner( winInd.beg );
    int colOwner = H.ColOwner( winInd.beg );

    Int thisCut = Mod( H.ColCut()+winInd.beg, H.BlockHeight() );
    Int winIndex = 0;
    while( winIndex < winSize )
    {
        const Int thisBlockSize = Min( blockSize-thisCut, winSize-winIndex );
        const int nextRowOwner = Mod( rowOwner+1, grid.Height() );
        const int nextColOwner = Mod( colOwner+1, grid.Width() );

        //  ---------------------------
        // | m x x x x x | x x x x x x |
        // | s m x x x x | x x x x x x |
        // |   s m x x x | x x x x x x |
        // |     s m x x | x x x x x x |
        // |       s m x | x x x x x x |
        // |         s m | x x x x x x |
        // |-------------|-------------|
        // |           s | m x x x x x |
        // |             | s m x x x x |
        // |             |   s m x x x |
        // |             |     s m x x |
        // |             |       s m x |
        // |             |         s m |
        //  ---------------------------
        const bool haveLastOff = (winSize > winIndex+thisBlockSize);

        // Handle this main diagonal block
        if( colOwner == grid.Col() && rowOwner == grid.Row() )
        {
            for( Int offset=0; offset<thisBlockSize; ++offset )
            {
                hMainWin.QueueUpdate
                ( winIndex+offset, 0,
                  HLoc(localRowOffset+offset,localColOffset+offset) );
                if( offset < thisBlockSize-1 )
                    hSubWin.QueueUpdate
                    ( winIndex+offset, 0,
                      HLoc(localRowOffset+offset+1,localColOffset+offset) );
            }
        }
        // Handle the last subdiagonal (if it exists)
        if( haveLastOff &&
            nextRowOwner == grid.Row() && colOwner == grid.Col() )
        {
            const Int subLocalRowOffset = localRowOffset +
              ( grid.Height() == 1 ? thisBlockSize : 0 );
            hSubWin.QueueUpdate
            ( winIndex+thisBlockSize-1, 0,
              HLoc(subLocalRowOffset,localColOffset+thisBlockSize-1) );
        }

        winIndex += thisBlockSize;
        thisCut = 0;
        if( colOwner == grid.Col() )
            localColOffset += thisBlockSize;
        if( rowOwner == grid.Row() )
            localRowOffset += thisBlockSize;
        rowOwner = nextRowOwner;
        colOwner = nextColOwner;
    }
    hMainWin.ProcessQueues();
    hSubWin.ProcessQueues();
}

template<typename Field>
void GatherTridiagonal
( const DistMatrix<Field,MC,MR,BLOCK>& H,
  const IR& winInd,
        DistMatrix<Field,STAR,STAR>& hMainWin,
        DistMatrix<Field,STAR,STAR>& hSubWin,
        DistMatrix<Field,STAR,STAR>& hSuperWin )
{
    EL_DEBUG_CSE
    const Int winSize = winInd.end - winInd.beg;
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();
    const auto& HLoc = H.LockedMatrix();
    EL_DEBUG_ONLY(
      if( H.BlockHeight() != H.BlockWidth() )
          LogicError("Assumed square distribution blocks");
      if( H.ColCut() != H.RowCut() )
          LogicError("Assumed symmetric cuts");
      if( blockSize < 2 )
          LogicError("Assumed blocks of size at least two");
    )
    Zeros( hMainWin, winSize, 1 );
    Zeros( hSubWin, winSize-1, 1 );
    Zeros( hSuperWin, winSize-1, 1 );

    // TODO(poulson): Optimize these redistributions
    hMainWin.Reserve( winSize );
    hSubWin.Reserve( winSize-1 );
    hSuperWin.Reserve( winSize-1 );

    Int localRowOffset = H.LocalRowOffset( winInd.beg );
    Int localColOffset = H.LocalColOffset( winInd.beg );
    int rowOwner = H.RowOwner( winInd.beg );
    int colOwner = H.ColOwner( winInd.beg );

    Int thisCut = Mod( H.ColCut()+winInd.beg, H.BlockHeight() );
    Int winIndex = 0;
    while( winIndex < winSize )
    {
        const Int thisBlockSize = Min( blockSize-thisCut, winSize-winIndex );
        const int nextRowOwner = Mod( rowOwner+1, grid.Height() );
        const int nextColOwner = Mod( colOwner+1, grid.Width() );

        //  ---------------------------
        // | m S x x x x | x x x x x x |
        // | s m S x x x | x x x x x x |
        // |   s m S x x | x x x x x x |
        // |     s m S x | x x x x x x |
        // |       s m S | x x x x x x |
        // |         s m | S x x x x x |
        // |-------------|-------------|
        // |           s | m S x x x x |
        // |             | s m S x x x |
        // |             |   s m S x x |
        // |             |     s m S x |
        // |             |       s m S |
        // |             |         s m |
        //  ---------------------------
        const bool haveLastOff = (winSize > winIndex+thisBlockSize);

        // Handle this main diagonal block
        if( colOwner == grid.Col() && rowOwner == grid.Row() )
        {
            for( Int offset=0; offset<thisBlockSize; ++offset )
            {
                hMainWin.QueueUpdate
                ( winIndex+offset, 0,
                  HLoc(localRowOffset+offset,localColOffset+offset) );
                if( offset < thisBlockSize-1 )
                    hSubWin.QueueUpdate
                    ( winIndex+offset, 0,
                      HLoc(localRowOffset+offset+1,localColOffset+offset) );
                if( offset < thisBlockSize-1 )
                    hSuperWin.QueueUpdate
                    ( winIndex+offset, 0,
                      HLoc(localRowOffset+offset,localColOffset+offset+1) );
            }
        }
        // Handle the last subdiagonal (if it exists)
        if( haveLastOff &&
            nextRowOwner == grid.Row() && colOwner == grid.Col() )
        {
            const Int subLocalRowOffset = localRowOffset +
              ( grid.Height() == 1 ? thisBlockSize : 0 );
            hSubWin.QueueUpdate
            ( winIndex+thisBlockSize-1, 0,
              HLoc(subLocalRowOffset,localColOffset+thisBlockSize-1) );
        }
        // Handle the last superdiagonal (if it exists)
        if( haveLastOff &&
            rowOwner == grid.Row() && nextColOwner == grid.Col() )
        {
            const Int subLocalColOffset = localColOffset +
              ( grid.Width() == 1 ? thisBlockSize : 0 );
            hSuperWin.QueueUpdate
            ( winIndex+thisBlockSize-1, 0,
              HLoc(localRowOffset+thisBlockSize-1,subLocalColOffset) );
        }

        winIndex += thisBlockSize;
        thisCut = 0;
        if( colOwner == grid.Col() )
            localColOffset += thisBlockSize;
        if( rowOwner == grid.Row() )
            localRowOffset += thisBlockSize;
        rowOwner = nextRowOwner;
        colOwner = nextColOwner;
    }
    hMainWin.ProcessQueues();
    hSubWin.ProcessQueues();
    hSuperWin.ProcessQueues();
}

} // namespace util
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_UTIL_GATHER_HPP
