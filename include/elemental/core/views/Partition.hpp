/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_VIEWS_PARTITION_HPP
#define ELEM_VIEWS_PARTITION_HPP

#include "./View.hpp"

namespace elem {

// To make our life easier. Undef'd at the bottom of the header
#define M  Matrix<T>
#define DM DistMatrix<T,U,V>

// Partition downwards from the top
// ================================

template<typename T>
inline void
PartitionDown
( M& A, M& AT,
        M& AB, Int heightAT=Blocksize() ) 
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    View( AT, A, 0,        0, heightAT, A.Width() );
    View( AB, A, heightAT, 0, heightAB, A.Width() );
}

template<typename T, Dist U, Dist V>
inline void
PartitionDown
( DM& A, DM& AT,
         DM& AB, Int heightAT=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    View( AT, A, 0,        0, heightAT, A.Width() );
    View( AB, A, heightAT, 0, heightAB, A.Width() );
}

template<typename T>
inline void
LockedPartitionDown
( const M& A, M& AT,
              M& AB, Int heightAT=Blocksize() ) 
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    LockedView( AT, A, 0,        0, heightAT, A.Width() );
    LockedView( AB, A, heightAT, 0, heightAB, A.Width() );
}

template<typename T, Dist U, Dist V>
inline void
LockedPartitionDown
( const DM& A, DM& AT,
               DM& AB, Int heightAT=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    LockedView( AT, A, 0,        0, heightAT, A.Width() );
    LockedView( AB, A, heightAT, 0, heightAB, A.Width() );
}

// Partition upwards from the bottom
// =================================

template<typename T>
inline void
PartitionUp
( M& A, M& AT,
        M& AB, Int heightAB=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUp"))
    PartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T, Dist U, Dist V>
inline void
PartitionUp
( DM& A, DM& AT,
         DM& AB, Int heightAB=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUp"))
    PartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T>
inline void
LockedPartitionUp
( const M& A, M& AT,
              M& AB, Int heightAB=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUp"))
    LockedPartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T, Dist U, Dist V>
inline void
LockedPartitionUp
( const DM& A, DM& AT,
               DM& AB, Int heightAB=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUp"))
    LockedPartitionDown( A, AT, AB, A.Height()-heightAB );
}

// Partition rightwards from the left
// ==================================

template<typename T>
inline void
PartitionRight( M& A, M& AL, M& AR, Int widthAL=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    View( AL, A, 0, 0,       A.Height(), widthAL );
    View( AR, A, 0, widthAL, A.Height(), widthAR );
}

template<typename T, Dist U, Dist V>
inline void
PartitionRight( DM& A, DM& AL, DM& AR, Int widthAL=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    View( AL, A, 0, 0,       A.Height(), widthAL );
    View( AR, A, 0, widthAL, A.Height(), widthAR );
}

template<typename T>
inline void
LockedPartitionRight( const M& A, M& AL, M& AR, Int widthAL=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    LockedView( AL, A, 0, 0,       A.Height(), widthAL );
    LockedView( AR, A, 0, widthAL, A.Height(), widthAR );
}

template<typename T, Dist U, Dist V>
inline void
LockedPartitionRight( const DM& A, DM& AL, DM& AR, Int widthAL=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    LockedView( AL, A, 0, 0,       A.Height(), widthAL );
    LockedView( AR, A, 0, widthAL, A.Height(), widthAR );
}

// Partition leftwards from the right
// ==================================

template<typename T>
inline void
PartitionLeft( M& A, M& AL, M& AR, Int widthAR=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionLeft"))
    PartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T, Dist U, Dist V>
inline void
PartitionLeft( DM& A, DM& AL, DM& AR, Int widthAR=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionLeft"))
    PartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T>
inline void
LockedPartitionLeft( const M& A, M& AL, M& AR, Int widthAR=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionLeft"))
    LockedPartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T, Dist U, Dist V>
inline void
LockedPartitionLeft( const DM& A, DM& AL, DM& AR, Int widthAR=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionLeft"))
    LockedPartitionRight( A, AL, AR, A.Width()-widthAR );
}

// Partition downward on a particular diagonal
// ===========================================

template<typename T>
inline void
PartitionDownOffsetDiagonal
( Int offset,
  M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownOffsetDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    diagDist = Max(Min(diagDist,diagLength),0);
    
    const Int mCut = ( offset<=0 ? -offset+diagDist : diagDist );
    const Int nCut = ( offset<=0 ? diagDist : offset+diagDist );
    View( ATL, A, 0,    0,    mCut,   nCut   );
    View( ATR, A, 0,    nCut, mCut,   n-nCut );
    View( ABL, A, mCut, 0,    m-mCut, nCut   );
    View( ABR, A, mCut, nCut, m-mCut, n-nCut );
}

template<typename T, Dist U, Dist V>
inline void
PartitionDownOffsetDiagonal
( Int offset,
  DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownOffsetDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    diagDist = Max(Min(diagDist,diagLength),0);

    const Int mCut = ( offset<=0 ? -offset+diagDist : diagDist );
    const Int nCut = ( offset<=0 ? diagDist : offset+diagDist );
    View( ATL, A, 0,    0,    mCut,   nCut   );
    View( ATR, A, 0,    nCut, mCut,   n-nCut );
    View( ABL, A, mCut, 0,    m-mCut, nCut   );
    View( ABR, A, mCut, nCut, m-mCut, n-nCut );
}

template<typename T>
inline void
LockedPartitionDownOffsetDiagonal
( Int offset,
  const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownOffsetDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    diagDist = Max(Min(diagDist,diagLength),0);
    
    const Int mCut = ( offset<=0 ? -offset+diagDist : diagDist );
    const Int nCut = ( offset<=0 ? diagDist : offset+diagDist );
    LockedView( ATL, A, 0,    0,    mCut,   nCut   );
    LockedView( ATR, A, 0,    nCut, mCut,   n-nCut );
    LockedView( ABL, A, mCut, 0,    m-mCut, nCut   );
    LockedView( ABR, A, mCut, nCut, m-mCut, n-nCut );
}

template<typename T, Dist U, Dist V>
inline void
LockedPartitionDownOffsetDiagonal
( Int offset,
  const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownOffsetDiagonal"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int diagLength = A.DiagonalLength(offset);
    diagDist = Max(Min(diagDist,diagLength),0);
    
    const Int mCut = ( offset<=0 ? -offset+diagDist : diagDist );
    const Int nCut = ( offset<=0 ? diagDist : offset+diagDist );
    LockedView( ATL, A, 0,    0,    mCut,   nCut   );
    LockedView( ATR, A, 0,    nCut, mCut,   n-nCut );
    LockedView( ABL, A, mCut, 0,    m-mCut, nCut   );
    LockedView( ABR, A, mCut, nCut, m-mCut, n-nCut );
}

// Partition upwards on a particular diagonal
// ==========================================

template<typename T>
inline void
PartitionUpOffsetDiagonal
( Int offset,
  M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpOffsetDiagonal"))
    PartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T, Dist U, Dist V>
inline void
PartitionUpOffsetDiagonal
( Int offset,
  DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpOffsetDiagonal"))
    PartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T>
inline void
LockedPartitionUpOffsetDiagonal
( Int offset,
  const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpOffsetDiagonal"))
    LockedPartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T, Dist U, Dist V>
inline void
LockedPartitionUpOffsetDiagonal
( Int offset,
  const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpOffsetDiagonal"))
    LockedPartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

// Partition downwards on the main diagonal
// ========================================

template<typename T>
inline void
PartitionDownDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownDiagonal"))
    PartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T, Dist U, Dist V>
inline void
PartitionDownDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownDiagonal"))
    PartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
inline void
LockedPartitionDownDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownDiagonal"))
    LockedPartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T, Dist U, Dist V>
inline void
LockedPartitionDownDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownDiagonal"))
    LockedPartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

// Partition upwards on the main diagonal
// ======================================

template<typename T>
inline void
PartitionUpDiagonal
( M& A, M& ATL, M& ATR,
        M& ABL, M& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpDiagonal"))
    PartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T, Dist U, Dist V>
inline void
PartitionUpDiagonal
( DM& A, DM& ATL, DM& ATR,
         DM& ABL, DM& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpDiagonal"))
    PartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
inline void
LockedPartitionUpDiagonal
( const M& A, M& ATL, M& ATR,
              M& ABL, M& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpDiagonal"))
    LockedPartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T, Dist U, Dist V>
inline void
LockedPartitionUpDiagonal
( const DM& A, DM& ATL, DM& ATR,
               DM& ABL, DM& ABR, Int diagDist=Blocksize() )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpDiagonal"))
    LockedPartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

#undef DM
#undef M

} // namespace elem

#endif // ifndef ELEM_VIEWS_PARTITION_HPP
