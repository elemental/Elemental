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

template<typename F>
void GatherSubdiagonal
( const DistMatrix<F,MC,MR,BLOCK>& H,
  const IR& winInd,
        DistMatrix<Base<F>,STAR,STAR>& hSubWin )
{
    DEBUG_CSE
    const Int winSize = winInd.end - winInd.beg;
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();
    const auto& HLoc = H.LockedMatrix();
    DEBUG_ONLY(
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

        if( colOwner == grid.Col() )
        {
            const bool haveLastSub = (winSize > winIndex+thisBlockSize);
            if( rowOwner == grid.Row() )
            {
                for( Int offset=0; offset<thisBlockSize; ++offset )
                {
                    if( offset < thisBlockSize-1 ||
                        (haveLastSub && grid.Height()==1) )
                        hSubWin.QueueUpdate
                        ( winIndex+offset, 0,
                          RealPart(HLoc(localRowOffset+offset+1,
                                        localColOffset+offset)) );
                }
            }
            else if( haveLastSub && nextRowOwner == grid.Row() )
            {
                // Grab the entry in the top-right of this block
                hSubWin.QueueUpdate
                ( winIndex+thisBlockSize-1, 0,
                  RealPart(HLoc(localRowOffset,
                                localColOffset+thisBlockSize-1)) );
            }
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

template<typename F>
void GatherBidiagonal
( const DistMatrix<F,MC,MR,BLOCK>& H,
  const IR& winInd,
        DistMatrix<F,STAR,STAR>& hMainWin,
        DistMatrix<Base<F>,STAR,STAR>& hSubWin )
{
    DEBUG_CSE
    const Int winSize = winInd.end - winInd.beg;
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();
    const auto& HLoc = H.LockedMatrix();
    DEBUG_ONLY(
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

        if( colOwner == grid.Col() )
        {
            const bool haveLastSub = (winSize > winIndex+thisBlockSize);
            if( rowOwner == grid.Row() )
            {
                for( Int offset=0; offset<thisBlockSize; ++offset )
                {
                    hMainWin.QueueUpdate
                    ( winIndex+offset, 0,
                      HLoc(localRowOffset+offset,localColOffset+offset) );
                    if( offset < thisBlockSize-1 ||
                        (haveLastSub && grid.Height()==1) )
                        hSubWin.QueueUpdate
                        ( winIndex+offset, 0,
                          RealPart(HLoc(localRowOffset+offset+1,
                                        localColOffset+offset)) );
                }
            }
            else if( haveLastSub && nextRowOwner == grid.Row() )
            {
                // Grab the entry in the top-right of this block
                hSubWin.QueueUpdate
                ( winIndex+thisBlockSize-1, 0,
                  RealPart(HLoc(localRowOffset,
                                localColOffset+thisBlockSize-1)) );
            }
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

template<typename F>
void GatherTridiagonal
( const DistMatrix<F,MC,MR,BLOCK>& H,
  const IR& winInd,
        DistMatrix<F,STAR,STAR>& hMainWin,
        DistMatrix<Base<F>,STAR,STAR>& hSubWin,
        DistMatrix<F,STAR,STAR>& hSuperWin )
{
    DEBUG_CSE
    const Int winSize = winInd.end - winInd.beg;
    const Int blockSize = H.BlockHeight();
    const Grid& grid = H.Grid();
    const auto& HLoc = H.LockedMatrix();
    DEBUG_ONLY(
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

        if( colOwner == grid.Col() )
        {
            const bool haveLastOff = (winSize > winIndex+thisBlockSize);
            if( rowOwner == grid.Row() )
            {
                for( Int offset=0; offset<thisBlockSize; ++offset )
                {
                    hMainWin.QueueUpdate
                    ( winIndex+offset, 0,
                      HLoc(localRowOffset+offset,localColOffset+offset) );
                    if( offset < thisBlockSize-1 ||
                        (haveLastOff && grid.Height()==1) )
                        hSubWin.QueueUpdate
                        ( winIndex+offset, 0,
                          RealPart(HLoc(localRowOffset+offset+1,
                                        localColOffset+offset)) );
                    if( offset < thisBlockSize-1 ||
                        (haveLastOff && grid.Width()==1) )
                        hSuperWin.QueueUpdate
                        ( winIndex+offset, 0,
                          HLoc(localRowOffset+offset,localColOffset+offset+1) );
                }
            }
            else if( haveLastOff && nextRowOwner == grid.Row() )
            {
                // Grab the entry in the top-right of this block
                hSubWin.QueueUpdate
                ( winIndex+thisBlockSize-1, 0,
                  RealPart(HLoc(localRowOffset,
                                localColOffset+thisBlockSize-1)) );
            }
        }
        else if( rowOwner == grid.Row() && nextColOwner == grid.Col() )
        {
            const bool haveLastOff = (winSize > winIndex+thisBlockSize);
            if( haveLastOff )
            {
                // Grab the entry in the bottom-left of this block
                hSuperWin.QueueUpdate
                ( winIndex+thisBlockSize-1, 0,
                  HLoc(localRowOffset+thisBlockSize-1,localColOffset) );
            }
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
