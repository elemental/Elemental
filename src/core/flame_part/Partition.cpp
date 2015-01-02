/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Partition downwards from the top
// ================================

template<typename T>
void PartitionDown
( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAT ) 
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    View( AT, A, 0,        0, heightAT, A.Width() );
    View( AB, A, heightAT, 0, heightAB, A.Width() );
}

template<typename T>
void PartitionDown
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, 
  Int heightAT )
{
    DEBUG_ONLY(
      CallStackEntry cse("PartitionDown");
      AssertSameGrids( A, AT, AB );
      AssertSameDists( A, AT, AB );
    )
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    View( AT, A, 0,        0, heightAT, A.Width() );
    View( AB, A, heightAT, 0, heightAB, A.Width() );
}

template<typename T>
void LockedPartitionDown
( const Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAT ) 
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDown"))
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    LockedView( AT, A, 0,        0, heightAT, A.Width() );
    LockedView( AB, A, heightAT, 0, heightAB, A.Width() );
}

template<typename T>
void LockedPartitionDown
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, 
  Int heightAT )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedPartitionDown");
      AssertSameGrids( A, AT, AB );
      AssertSameDists( A, AT, AB );
    )
    heightAT = Max(Min(heightAT,A.Height()),0);
    const Int heightAB = A.Height()-heightAT;
    LockedView( AT, A, 0,        0, heightAT, A.Width() );
    LockedView( AB, A, heightAT, 0, heightAB, A.Width() );
}

// Partition upwards from the bottom
// =================================

template<typename T>
void PartitionUp
( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAB )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUp"))
    PartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T>
void PartitionUp
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, 
  Int heightAB )
{
    DEBUG_ONLY(
      CallStackEntry cse("PartitionUp");
      AssertSameGrids( A, AT, AB );
      AssertSameDists( A, AT, AB );
    )
    PartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T>
void LockedPartitionUp
( const Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAB )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUp"))
    LockedPartitionDown( A, AT, AB, A.Height()-heightAB );
}

template<typename T>
void LockedPartitionUp
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, 
  Int heightAB )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedPartitionUp");
      AssertSameGrids( A, AT, AB );
      AssertSameDists( A, AT, AB );
    )
    LockedPartitionDown( A, AT, AB, A.Height()-heightAB );
}

// Partition rightwards from the left
// ==================================

template<typename T>
void PartitionRight
( Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAL )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    View( AL, A, 0, 0,       A.Height(), widthAL );
    View( AR, A, 0, widthAL, A.Height(), widthAR );
}

template<typename T>
void PartitionRight
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, 
  Int widthAL )
{
    DEBUG_ONLY(
      CallStackEntry cse("PartitionRight");
      AssertSameGrids( A, AL, AR );
      AssertSameDists( A, AL, AR );
    )
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    View( AL, A, 0, 0,       A.Height(), widthAL );
    View( AR, A, 0, widthAL, A.Height(), widthAR );
}

template<typename T>
void LockedPartitionRight
( const Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAL )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionRight"))
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    LockedView( AL, A, 0, 0,       A.Height(), widthAL );
    LockedView( AR, A, 0, widthAL, A.Height(), widthAR );
}

template<typename T>
void LockedPartitionRight
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, 
  Int widthAL )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedPartitionRight");
      AssertSameGrids( A, AL, AR );
      AssertSameDists( A, AL, AR );
    )
    widthAL = Max(Min(widthAL,A.Width()),0);
    const Int widthAR = A.Width()-widthAL;
    LockedView( AL, A, 0, 0,       A.Height(), widthAL );
    LockedView( AR, A, 0, widthAL, A.Height(), widthAR );
}

// Partition leftwards from the right
// ==================================

template<typename T>
void PartitionLeft
( Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAR )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionLeft"))
    PartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T>
void PartitionLeft
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, 
  Int widthAR )
{
    DEBUG_ONLY(
      CallStackEntry cse("PartitionLeft");
      AssertSameGrids( A, AL, AR );
      AssertSameDists( A, AL, AR );
    )
    PartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T>
void LockedPartitionLeft
( const Matrix<T>& A, 
        Matrix<T>& AL, Matrix<T>& AR, Int widthAR )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionLeft"))
    LockedPartitionRight( A, AL, AR, A.Width()-widthAR );
}

template<typename T>
void LockedPartitionLeft
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, 
  Int widthAR )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedPartitionLeft");
      AssertSameGrids( A, AL, AR );
      AssertSameDists( A, AL, AR );
    )
    LockedPartitionRight( A, AL, AR, A.Width()-widthAR );
}

// Partition downward on a particular diagonal
// ===========================================

template<typename T>
void PartitionDownOffsetDiagonal
( Int offset,
  Matrix<T>& A, 
  Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist )
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
void PartitionDownOffsetDiagonal
( Int offset,
  AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist )
{
    DEBUG_ONLY(
      CallStackEntry cse("PartitionDownOffsetDiagonal");
      AssertSameGrids( A, ATL, ATR, ABL, ABR );
      AssertSameDists( A, ATL, ATR, ABL, ABR );
    )
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
void LockedPartitionDownOffsetDiagonal
( Int offset,
  const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist )
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

template<typename T>
void LockedPartitionDownOffsetDiagonal
( Int offset,
  const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
        AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedPartitionDownOffsetDiagonal");
      AssertSameGrids( A, ATL, ATR, ABL, ABR );
      AssertSameDists( A, ATL, ATR, ABL, ABR );
    )
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
void PartitionUpOffsetDiagonal
( Int offset,
  Matrix<T>& A, 
  Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpOffsetDiagonal"))
    PartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T>
void PartitionUpOffsetDiagonal
( Int offset,
  AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist )
{
    DEBUG_ONLY(
      CallStackEntry cse("PartitionUpOffsetDiagonal");
      AssertSameGrids( A, ATL, ATR, ABL, ABR );
      AssertSameDists( A, ATL, ATR, ABL, ABR );
    )
    PartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T>
void LockedPartitionUpOffsetDiagonal
( Int offset,
  const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpOffsetDiagonal"))
    LockedPartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

template<typename T>
void LockedPartitionUpOffsetDiagonal
( Int offset,
  const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
        AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedPartitionUpOffsetDiagonal");
      AssertSameGrids( A, ATL, ATR, ABL, ABR );
      AssertSameDists( A, ATL, ATR, ABL, ABR );
    )
    LockedPartitionDownOffsetDiagonal
    ( offset, A, ATL, ATR, ABL, ABR, A.DiagonalLength(offset)-diagDist );
}

// Partition downwards on the main diagonal
// ========================================

template<typename T>
void PartitionDownDiagonal
( Matrix<T>& A, 
  Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionDownDiagonal"))
    PartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
void PartitionDownDiagonal
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist )
{
    DEBUG_ONLY(
      CallStackEntry cse("PartitionDownDiagonal");
      AssertSameGrids( A, ATL, ATR, ABL, ABR );
      AssertSameDists( A, ATL, ATR, ABL, ABR );
    )
    PartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
void LockedPartitionDownDiagonal
( const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionDownDiagonal"))
    LockedPartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
void LockedPartitionDownDiagonal
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
        AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedPartitionDownDiagonal");
      AssertSameGrids( A, ATL, ATR, ABL, ABR );
      AssertSameDists( A, ATL, ATR, ABL, ABR );
    )
    LockedPartitionDownOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

// Partition upwards on the main diagonal
// ======================================

template<typename T>
void PartitionUpDiagonal
( Matrix<T>& A, 
  Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("PartitionUpDiagonal"))
    PartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
void PartitionUpDiagonal
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist )
{
    DEBUG_ONLY(
      CallStackEntry cse("PartitionUpDiagonal");
      AssertSameGrids( A, ATL, ATR, ABL, ABR );
      AssertSameDists( A, ATL, ATR, ABL, ABR );
    )
    PartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
void LockedPartitionUpDiagonal
( const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist )
{
    DEBUG_ONLY(CallStackEntry cse("LockedPartitionUpDiagonal"))
    LockedPartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

template<typename T>
void LockedPartitionUpDiagonal
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
        AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedPartitionUpDiagonal");
      AssertSameGrids( A, ATL, ATR, ABL, ABR );
      AssertSameDists( A, ATL, ATR, ABL, ABR );
    )
    LockedPartitionUpOffsetDiagonal( 0, A, ATL, ATR, ABL, ABR, diagDist );
}

#define PROTO(T) \
  /* Downwards */ \
  template void PartitionDown \
  ( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAT ); \
  template void LockedPartitionDown \
  ( const Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAT ); \
  template void PartitionDown \
  ( AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, Int heightAT ); \
  template void LockedPartitionDown \
  ( const AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, Int heightAT ); \
  /* Upwards */ \
  template void PartitionUp \
  ( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAB ); \
  template void LockedPartitionUp \
  ( const Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAB ); \
  template void PartitionUp \
  ( AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, Int heightAB ); \
  template void LockedPartitionUp \
  ( const AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, Int heightAB ); \
  /* Leftward */ \
  template void PartitionLeft \
  ( Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAR ); \
  template void LockedPartitionLeft \
  ( const Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAR ); \
  template void PartitionLeft \
  ( AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, Int widthAR ); \
  template void LockedPartitionLeft \
  ( const AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, Int widthAR ); \
  /* Rightward */ \
  template void PartitionRight \
  ( Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAL ); \
  template void LockedPartitionRight \
  ( const Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAL ); \
  template void PartitionRight \
  ( AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, Int widthAL ); \
  template void LockedPartitionRight \
  ( const AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, Int widthAL ); \
  /* Down the main diagonal */ \
  template void PartitionDownDiagonal \
  ( Matrix<T>& A, \
    Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist ); \
  template void LockedPartitionDownDiagonal \
  ( const Matrix<T>& A, \
    Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist ); \
  template void PartitionDownDiagonal \
  ( AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, Int diagDist ); \
  template void LockedPartitionDownDiagonal \
  ( const AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, Int diagDist ); \
  /* Down an offset diagonal */ \
  template void PartitionDownOffsetDiagonal \
  ( Int offset, Matrix<T>& A, \
    Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist ); \
  template void LockedPartitionDownOffsetDiagonal \
  ( Int offset, const Matrix<T>& A, \
    Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist ); \
  template void PartitionDownOffsetDiagonal \
  ( Int offset, AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, Int diagDist ); \
  template void LockedPartitionDownOffsetDiagonal \
  ( Int offset, const AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, Int diagDist ); \
  /* Up the main diagonal */ \
  template void PartitionUpDiagonal \
  ( Matrix<T>& A, \
    Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist ); \
  template void LockedPartitionUpDiagonal \
  ( const Matrix<T>& A, \
    Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist ); \
  template void PartitionUpDiagonal \
  ( AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, Int diagDist ); \
  template void LockedPartitionUpDiagonal \
  ( const AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, Int diagDist ); \
  /* Up an offset diagonal */ \
  template void PartitionUpOffsetDiagonal \
  ( Int offset, Matrix<T>& A, \
    Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist ); \
  template void LockedPartitionUpOffsetDiagonal \
  ( Int offset, const Matrix<T>& A, \
    Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist ); \
  template void PartitionUpOffsetDiagonal \
  ( Int offset, AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, Int diagDist ); \
  template void LockedPartitionUpOffsetDiagonal \
  ( Int offset, const AbstractDistMatrix<T>& A, \
    AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, Int diagDist );

#include "El/macros/Instantiate.h"

} // namespace El
