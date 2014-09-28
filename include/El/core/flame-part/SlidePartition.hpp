/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FLAMEPART_SLIDEPARTITION_HPP
#define EL_FLAMEPART_SLIDEPARTITION_HPP

namespace El {

// To make our life easier. Undef'd at the bottom of the header
#define M   Matrix<T>
#define ADM AbstractDistMatrix<T>

// Slide a partition upward
// ========================

template<typename T>
inline void
SlidePartitionUp
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionUp"))
    View( AT, A0 );
    Merge2x1( AB, A1, A2 );
}

template<typename T>
inline void
SlidePartitionUp
( ADM& AT, ADM& A0,
           ADM& A1,
  ADM& AB, ADM& A2 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlidePartitionUp");
      AssertSameGrids( AT, AB, A0, A1, A2 );
      AssertSameDists( AT, AB, A0, A1, A2 );
    )
    View( AT, A0 );
    Merge2x1( AB, A1, A2 );
}

template<typename T>
inline void
SlideLockedPartitionUp
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionUp"))
    LockedView( AT, A0 );
    LockedMerge2x1( AB, A1, A2 );
}

template<typename T>
inline void
SlideLockedPartitionUp
( ADM& AT, const ADM& A0,
           const ADM& A1,
  ADM& AB, const ADM& A2 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlideLockedPartitionUp");
      AssertSameGrids( AT, AB, A0, A1, A2 );
      AssertSameDists( AT, AB, A0, A1, A2 );
    )
    LockedView( AT, A0 );
    LockedMerge2x1( AB, A1, A2 );
}

// Slide a partition downward
// ==========================

template<typename T>
inline void
SlidePartitionDown
( M& AT, M& A0,
         M& A1,
  M& AB, M& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionDown"))
    Merge2x1( AT, A0, A1 );
    View( AB, A2 );
}

template<typename T>
inline void
SlidePartitionDown
( ADM& AT, ADM& A0,
           ADM& A1,
  ADM& AB, ADM& A2 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlidePartitionDown");
      AssertSameGrids( AT, AB, A0, A1, A2 );
      AssertSameDists( AT, AB, A0, A1, A2 );
    )
    Merge2x1( AT, A0, A1 );
    View( AB, A2 );
}

template<typename T>
inline void
SlideLockedPartitionDown
( M& AT, const M& A0,
         const M& A1,
  M& AB, const M& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionDown"))
    LockedMerge2x1( AT, A0, A1 );
    LockedView( AB, A2 );
}

template<typename T>
inline void
SlideLockedPartitionDown
( ADM& AT, const ADM& A0,
           const ADM& A1,
  ADM& AB, const ADM& A2 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlideLockedPartitionDown");
      AssertSameGrids( AT, AB, A0, A1, A2 );
      AssertSameDists( AT, AB, A0, A1, A2 );
    )
    LockedMerge2x1( AT, A0, A1 );
    LockedView( AB, A2 );
}

// Slide a partition leftward
// ==========================

template<typename T>
inline void
SlidePartitionLeft
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionLeft"))
    View( AL, A0 );
    Merge1x2( AR, A1, A2 );
}

template<typename T>
inline void
SlidePartitionLeft
( ADM& AL, ADM& AR,
  ADM& A0, ADM& A1, ADM& A2 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlidePartitionLeft");
      AssertSameGrids( AL, AR, A0, A1, A2 );
      AssertSameDists( AL, AR, A0, A1, A2 );
    )
    View( AL, A0 );
    Merge1x2( AR, A1, A2 );
}

template<typename T>
inline void
SlideLockedPartitionLeft
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionLeft"))
    LockedView( AL, A0 );
    LockedMerge1x2( AR, A1, A2 );
}

template<typename T>
inline void
SlideLockedPartitionLeft
( ADM& AL, ADM& AR,
  const ADM& A0, const ADM& A1, const ADM& A2 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlideLockedPartitionLeft");
      AssertSameGrids( AL, AR, A0, A1, A2 );
      AssertSameDists( AL, AR, A0, A1, A2 );
    )
    LockedView( AL, A0 );
    LockedMerge1x2( AR, A1, A2 );
}

// Slide a partition rightward
// ===========================

template<typename T>
inline void
SlidePartitionRight
( M& AL, M& AR,
  M& A0, M& A1, M& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionRight"))
    Merge1x2( AL, A0, A1 );
    View( AR, A2 );
}

template<typename T>
inline void
SlidePartitionRight
( ADM& AL, ADM& AR,
  ADM& A0, ADM& A1, ADM& A2 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlidePartitionRight");
      AssertSameGrids( AL, AR, A0, A1, A2 );
      AssertSameDists( AL, AR, A0, A1, A2 );
    )
    Merge1x2( AL, A0, A1 );
    View( AR, A2 );
}

template<typename T>
inline void
SlideLockedPartitionRight
( M& AL, M& AR,
  const M& A0, const M& A1, const M& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionRight"))
    LockedMerge1x2( AL, A0, A1 );
    LockedView( AR, A2 );
}

template<typename T>
inline void
SlideLockedPartitionRight
( ADM& AL, ADM& AR,
  const ADM& A0, const ADM& A1, const ADM& A2 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlideLockedPartitionRight");
      AssertSameGrids( AL, AR, A0, A1, A2 );
      AssertSameDists( AL, AR, A0, A1, A2 );
    )
    LockedMerge1x2( AL, A0, A1 );
    LockedView( AR, A2 );
}

// Slide a partition upward on a diagonal
// ======================================

template<typename T>
inline void
SlidePartitionUpDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionUpDiagonal"))
    View( ATL, A00 );
    Merge1x2( ATR, A01, A02 );
    Merge2x1( ABL, A10, A20 );
    Merge2x2( ABR, A11, A12,
                   A21, A22 );
}

template<typename T>
inline void
SlidePartitionUpDiagonal
( ADM& ATL, ADM& ATR, ADM& A00, ADM& A01, ADM& A02,
                      ADM& A10, ADM& A11, ADM& A12,
  ADM& ABL, ADM& ABR, ADM& A20, ADM& A21, ADM& A22 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlidePartitionUpDiagonal");
      AssertSameGrids
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
      AssertSameDists
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
    )
    View( ATL, A00 );
    Merge1x2( ATR, A01, A02 );
    Merge2x1( ABL, A10, A20 );
    Merge2x2( ABR, A11, A12,
                   A21, A22 );
}

template<typename T>
inline void
SlideLockedPartitionUpDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionUpDiagonal"))
    LockedView( ATL, A00 );
    LockedMerge1x2( ATR, A01, A02 );
    LockedMerge2x1( ABL, A10, A20 );
    LockedMerge2x2( ABR, A11, A12,
                         A21, A22 );
}

template<typename T>
inline void
SlideLockedPartitionUpDiagonal
( ADM& ATL, ADM& ATR, const ADM& A00, const ADM& A01, const ADM& A02,
                      const ADM& A10, const ADM& A11, const ADM& A12,
  ADM& ABL, ADM& ABR, const ADM& A20, const ADM& A21, const ADM& A22 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlideLockedPartitionUpDiagonal");
      AssertSameGrids
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
      AssertSameDists
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
    )
    LockedView( ATL, A00 );
    LockedMerge1x2( ATR, A01, A02 );
    LockedMerge2x1( ABL, A10, A20 );
    LockedMerge2x2( ABR, A11, A12,
                         A21, A22 );
}

// Slide a partition downward on a diagonal
// ========================================

template<typename T>
inline void
SlidePartitionDownDiagonal
( M& ATL, M& ATR, M& A00, M& A01, M& A02,
                  M& A10, M& A11, M& A12,
  M& ABL, M& ABR, M& A20, M& A21, M& A22 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionDownDiagonal"))
    Merge2x2( ATL, A00, A01,
                   A10, A11 );
    Merge2x1( ATR, A02, A12 );
    Merge1x2( ABL, A20, A21 );
    View( ABR, A22 );
}

template<typename T>
inline void
SlidePartitionDownDiagonal
( ADM& ATL, ADM& ATR, ADM& A00, ADM& A01, ADM& A02,
                      ADM& A10, ADM& A11, ADM& A12,
  ADM& ABL, ADM& ABR, ADM& A20, ADM& A21, ADM& A22 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlidePartitionDownDiagonal");
      AssertSameGrids
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
      AssertSameDists
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
    )
    Merge2x2( ATL, A00, A01,
                   A10, A11 );
    Merge2x1( ATR, A02, A12 );
    Merge1x2( ABL, A20, A21 );
    View( ABR, A22 );
}

template<typename T>
inline void
SlideLockedPartitionDownDiagonal
( M& ATL, M& ATR, const M& A00, const M& A01, const M& A02,
                  const M& A10, const M& A11, const M& A12,
  M& ABL, M& ABR, const M& A20, const M& A21, const M& A22 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionDownDiagonal"))
    LockedMerge2x2( ATL, A00, A01,
                         A10, A11 );
    LockedMerge2x1( ATR, A02, A12 );
    LockedMerge1x2( ABL, A20, A21 );
    LockedView( ABR, A22 );
}

template<typename T>
inline void
SlideLockedPartitionDownDiagonal
( ADM& ATL, ADM& ATR, const ADM& A00, const ADM& A01, const ADM& A02,
                      const ADM& A10, const ADM& A11, const ADM& A12,
  ADM& ABL, ADM& ABR, const ADM& A20, const ADM& A21, const ADM& A22 )
{
    DEBUG_ONLY(
      CallStackEntry cse("SlideLockedPartitionDownDiagonal");
      AssertSameGrids
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
      AssertSameDists
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
    )
    LockedMerge2x2( ATL, A00, A01,
                         A10, A11 );
    LockedMerge2x1( ATR, A02, A12 );
    LockedMerge1x2( ABL, A20, A21 );
    LockedView( ABR, A22 );
}

#undef ADM
#undef M

} // namespace El

#endif // ifndef EL_FLAMEPART_SLIDEPARTITION_HPP
