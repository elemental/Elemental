/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Slide a partition upward
// ========================

template<typename T>
void SlidePartitionUp
( Matrix<T>& AT, Matrix<T>& A0,
                 Matrix<T>& A1,
  Matrix<T>& AB, Matrix<T>& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionUp"))
    View( AT, A0 );
    Merge2x1( AB, A1, A2 );
}

template<typename T>
void SlidePartitionUp
( AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0,
                             AbstractDistMatrix<T>& A1,
  AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2 )
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
void SlideLockedPartitionUp
( Matrix<T>& AT, const Matrix<T>& A0,
                 const Matrix<T>& A1,
  Matrix<T>& AB, const Matrix<T>& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionUp"))
    LockedView( AT, A0 );
    LockedMerge2x1( AB, A1, A2 );
}

template<typename T>
void SlideLockedPartitionUp
( AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& A0,
                             const AbstractDistMatrix<T>& A1,
  AbstractDistMatrix<T>& AB, const AbstractDistMatrix<T>& A2 )
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
void SlidePartitionDown
( Matrix<T>& AT, Matrix<T>& A0,
                 Matrix<T>& A1,
  Matrix<T>& AB, Matrix<T>& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionDown"))
    Merge2x1( AT, A0, A1 );
    View( AB, A2 );
}

template<typename T>
void SlidePartitionDown
( AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0,
                             AbstractDistMatrix<T>& A1,
  AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2 )
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
void SlideLockedPartitionDown
( Matrix<T>& AT, const Matrix<T>& A0,
                 const Matrix<T>& A1,
  Matrix<T>& AB, const Matrix<T>& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionDown"))
    LockedMerge2x1( AT, A0, A1 );
    LockedView( AB, A2 );
}

template<typename T>
void SlideLockedPartitionDown
( AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& A0,
                             const AbstractDistMatrix<T>& A1,
  AbstractDistMatrix<T>& AB, const AbstractDistMatrix<T>& A2 )
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
void SlidePartitionLeft
( Matrix<T>& AL, Matrix<T>& AR,
  Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionLeft"))
    View( AL, A0 );
    Merge1x2( AR, A1, A2 );
}

template<typename T>
void SlidePartitionLeft
( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR,
  AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, AbstractDistMatrix<T>& A2 )
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
void SlideLockedPartitionLeft
( Matrix<T>& AL, Matrix<T>& AR,
  const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionLeft"))
    LockedView( AL, A0 );
    LockedMerge1x2( AR, A1, A2 );
}

template<typename T>
void SlideLockedPartitionLeft
( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR,
  const AbstractDistMatrix<T>& A0, const AbstractDistMatrix<T>& A1, const AbstractDistMatrix<T>& A2 )
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
void SlidePartitionRight
( Matrix<T>& AL, Matrix<T>& AR,
  Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionRight"))
    Merge1x2( AL, A0, A1 );
    View( AR, A2 );
}

template<typename T>
void SlidePartitionRight
( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR,
  AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, AbstractDistMatrix<T>& A2 )
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
void SlideLockedPartitionRight
( Matrix<T>& AL, Matrix<T>& AR,
  const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionRight"))
    LockedMerge1x2( AL, A0, A1 );
    LockedView( AR, A2 );
}

template<typename T>
void SlideLockedPartitionRight
( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR,
  const AbstractDistMatrix<T>& A0, const AbstractDistMatrix<T>& A1, const AbstractDistMatrix<T>& A2 )
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
void SlidePartitionUpDiagonal
( Matrix<T>& ATL, Matrix<T>& ATR, 
  Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02,
  Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12,
  Matrix<T>& ABL, Matrix<T>& ABR, 
  Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionUpDiagonal"))
    View( ATL, A00 );
    Merge1x2( ATR, A01, A02 );
    Merge2x1( ABL, A10, A20 );
    Merge2x2( ABR, A11, A12,
                   A21, A22 );
}

template<typename T>
void SlidePartitionUpDiagonal
( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, 
  AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, AbstractDistMatrix<T>& A02,
  AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, AbstractDistMatrix<T>& A12,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, AbstractDistMatrix<T>& A20, 
  AbstractDistMatrix<T>& A21, AbstractDistMatrix<T>& A22 )
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
void SlideLockedPartitionUpDiagonal
( Matrix<T>& ATL, Matrix<T>& ATR, 
  const Matrix<T>& A00, const Matrix<T>& A01, const Matrix<T>& A02,
  const Matrix<T>& A10, const Matrix<T>& A11, const Matrix<T>& A12,
  Matrix<T>& ABL, Matrix<T>& ABR, 
  const Matrix<T>& A20, const Matrix<T>& A21, const Matrix<T>& A22 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionUpDiagonal"))
    LockedView( ATL, A00 );
    LockedMerge1x2( ATR, A01, A02 );
    LockedMerge2x1( ABL, A10, A20 );
    LockedMerge2x2( ABR, A11, A12,
                         A21, A22 );
}

template<typename T>
void SlideLockedPartitionUpDiagonal
( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, 
  const AbstractDistMatrix<T>& A00, const AbstractDistMatrix<T>& A01, const AbstractDistMatrix<T>& A02,
  const AbstractDistMatrix<T>& A10, const AbstractDistMatrix<T>& A11, const AbstractDistMatrix<T>& A12,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  const AbstractDistMatrix<T>& A20, const AbstractDistMatrix<T>& A21, const AbstractDistMatrix<T>& A22 )
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
void SlidePartitionDownDiagonal
( Matrix<T>& ATL, Matrix<T>& ATR, 
  Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02,
  Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12,
  Matrix<T>& ABL, Matrix<T>& ABR, 
  Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22 )
{
    DEBUG_ONLY(CallStackEntry cse("SlidePartitionDownDiagonal"))
    Merge2x2( ATL, A00, A01,
                   A10, A11 );
    Merge2x1( ATR, A02, A12 );
    Merge1x2( ABL, A20, A21 );
    View( ABR, A22 );
}

template<typename T>
void SlidePartitionDownDiagonal
( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, 
  AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, AbstractDistMatrix<T>& A02,
  AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, AbstractDistMatrix<T>& A12,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, AbstractDistMatrix<T>& A22 )
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
void SlideLockedPartitionDownDiagonal
( Matrix<T>& ATL, Matrix<T>& ATR, 
  const Matrix<T>& A00, const Matrix<T>& A01, const Matrix<T>& A02,
  const Matrix<T>& A10, const Matrix<T>& A11, const Matrix<T>& A12,
  Matrix<T>& ABL, Matrix<T>& ABR, 
  const Matrix<T>& A20, const Matrix<T>& A21, const Matrix<T>& A22 )
{
    DEBUG_ONLY(CallStackEntry cse("SlideLockedPartitionDownDiagonal"))
    LockedMerge2x2( ATL, A00, A01,
                         A10, A11 );
    LockedMerge2x1( ATR, A02, A12 );
    LockedMerge1x2( ABL, A20, A21 );
    LockedView( ABR, A22 );
}

template<typename T>
void SlideLockedPartitionDownDiagonal
( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, 
  const AbstractDistMatrix<T>& A00, const AbstractDistMatrix<T>& A01, const AbstractDistMatrix<T>& A02,
  const AbstractDistMatrix<T>& A10, const AbstractDistMatrix<T>& A11, const AbstractDistMatrix<T>& A12,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  const AbstractDistMatrix<T>& A20, const AbstractDistMatrix<T>& A21, const AbstractDistMatrix<T>& A22 )
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

#define PROTO(T) \
  /* Downward */ \
  template void SlidePartitionDown \
  ( Matrix<T>& AT, Matrix<T>& A0, \
                   Matrix<T>& A1, \
    Matrix<T>& AB, Matrix<T>& A2 ); \
  template void SlidePartitionDown \
  ( AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0, \
                               AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2 ); \
  template void SlideLockedPartitionDown \
  ( Matrix<T>& AT, const Matrix<T>& A0, \
                   const Matrix<T>& A1, \
    Matrix<T>& AB, const Matrix<T>& A2 ); \
  template void SlideLockedPartitionDown \
  ( AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& A0, \
                               const AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& AB, const AbstractDistMatrix<T>& A2 ); \
  /* Upward */ \
  template void SlidePartitionUp \
  ( Matrix<T>& AT, Matrix<T>& A0, \
                   Matrix<T>& A1, \
    Matrix<T>& AB, Matrix<T>& A2 ); \
  template void SlidePartitionUp \
  ( AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0, \
                               AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2 ); \
  template void SlideLockedPartitionUp \
  ( Matrix<T>& AT, const Matrix<T>& A0, \
                   const Matrix<T>& A1, \
    Matrix<T>& AB, const Matrix<T>& A2 ); \
  template void SlideLockedPartitionUp \
  ( AbstractDistMatrix<T>& AT, const AbstractDistMatrix<T>& A0, \
                               const AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& AB, const AbstractDistMatrix<T>& A2 ); \
  /* Right */ \
  template void SlidePartitionRight \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2 ); \
  template void SlidePartitionRight \
  ( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, \
    AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, AbstractDistMatrix<T>& A2 ); \
  template void SlideLockedPartitionRight \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2 ); \
  template void SlideLockedPartitionRight \
  ( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, \
    const AbstractDistMatrix<T>& A0, const AbstractDistMatrix<T>& A1, const AbstractDistMatrix<T>& A2 ); \
  /* Left */ \
  template void SlidePartitionLeft \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2 ); \
  template void SlidePartitionLeft \
  ( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, \
    AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, AbstractDistMatrix<T>& A2 ); \
  template void SlideLockedPartitionLeft \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2 ); \
  template void SlideLockedPartitionLeft \
  ( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, \
    const AbstractDistMatrix<T>& A0, const AbstractDistMatrix<T>& A1, const AbstractDistMatrix<T>& A2 ); \
  /* Down diagonal */ \
  template void SlidePartitionDownDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22 ); \
  template void SlidePartitionDownDiagonal \
  ( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, AbstractDistMatrix<T>& A02, \
    AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, AbstractDistMatrix<T>& A12, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, \
    AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, AbstractDistMatrix<T>& A22 ); \
  template void SlideLockedPartitionDownDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    const Matrix<T>& A00, const Matrix<T>& A01, const Matrix<T>& A02, \
    const Matrix<T>& A10, const Matrix<T>& A11, const Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    const Matrix<T>& A20, const Matrix<T>& A21, const Matrix<T>& A22 ); \
  template void SlideLockedPartitionDownDiagonal \
  ( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    const AbstractDistMatrix<T>& A00, const AbstractDistMatrix<T>& A01, const AbstractDistMatrix<T>& A02, \
    const AbstractDistMatrix<T>& A10, const AbstractDistMatrix<T>& A11, const AbstractDistMatrix<T>& A12, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, \
    const AbstractDistMatrix<T>& A20, const AbstractDistMatrix<T>& A21, const AbstractDistMatrix<T>& A22 ); \
  /* Up diagonal */ \
  template void SlidePartitionUpDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22 ); \
  template void SlidePartitionUpDiagonal \
  ( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, AbstractDistMatrix<T>& A02, \
    AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, AbstractDistMatrix<T>& A12, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, \
    AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, AbstractDistMatrix<T>& A22 ); \
  template void SlideLockedPartitionUpDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    const Matrix<T>& A00, const Matrix<T>& A01, const Matrix<T>& A02, \
    const Matrix<T>& A10, const Matrix<T>& A11, const Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    const Matrix<T>& A20, const Matrix<T>& A21, const Matrix<T>& A22 ); \
  template void SlideLockedPartitionUpDiagonal \
  ( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    const AbstractDistMatrix<T>& A00, const AbstractDistMatrix<T>& A01, const AbstractDistMatrix<T>& A02, \
    const AbstractDistMatrix<T>& A10, const AbstractDistMatrix<T>& A11, const AbstractDistMatrix<T>& A12, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, \
    const AbstractDistMatrix<T>& A20, const AbstractDistMatrix<T>& A21, const AbstractDistMatrix<T>& A22 );

#include "El/macros/Instantiate.h"

} // namespace El
