/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>

namespace El {

// Slide a partition upward
// ========================

template<typename T>
void SlidePartitionUp
( Matrix<T>& AT, Matrix<T>& A0,
                 Matrix<T>& A1,
  Matrix<T>& AB, Matrix<T>& A2 )
{
    DEBUG_CSE
    View( AT, A0 );
    Merge2x1( AB, A1, A2 );
}

template<typename T>
void SlidePartitionUp
( ElementalMatrix<T>& AT, ElementalMatrix<T>& A0,
                          ElementalMatrix<T>& A1,
  ElementalMatrix<T>& AB, ElementalMatrix<T>& A2 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    LockedView( AT, A0 );
    LockedMerge2x1( AB, A1, A2 );
}

template<typename T>
void SlideLockedPartitionUp
( ElementalMatrix<T>& AT, const ElementalMatrix<T>& A0,
                          const ElementalMatrix<T>& A1,
  ElementalMatrix<T>& AB, const ElementalMatrix<T>& A2 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    Merge2x1( AT, A0, A1 );
    View( AB, A2 );
}

template<typename T>
void SlidePartitionDown
( ElementalMatrix<T>& AT, ElementalMatrix<T>& A0,
                          ElementalMatrix<T>& A1,
  ElementalMatrix<T>& AB, ElementalMatrix<T>& A2 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    LockedMerge2x1( AT, A0, A1 );
    LockedView( AB, A2 );
}

template<typename T>
void SlideLockedPartitionDown
( ElementalMatrix<T>& AT, const ElementalMatrix<T>& A0,
                          const ElementalMatrix<T>& A1,
  ElementalMatrix<T>& AB, const ElementalMatrix<T>& A2 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    View( AL, A0 );
    Merge1x2( AR, A1, A2 );
}

template<typename T>
void SlidePartitionLeft
( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR,
  ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, ElementalMatrix<T>& A2 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    LockedView( AL, A0 );
    LockedMerge1x2( AR, A1, A2 );
}

template<typename T>
void SlideLockedPartitionLeft
(       ElementalMatrix<T>& AL,
        ElementalMatrix<T>& AR,
  const ElementalMatrix<T>& A0,
  const ElementalMatrix<T>& A1,
  const ElementalMatrix<T>& A2 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    Merge1x2( AL, A0, A1 );
    View( AR, A2 );
}

template<typename T>
void SlidePartitionRight
( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR,
  ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, ElementalMatrix<T>& A2 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    LockedMerge1x2( AL, A0, A1 );
    LockedView( AR, A2 );
}

template<typename T>
void SlideLockedPartitionRight
(       ElementalMatrix<T>& AL,
        ElementalMatrix<T>& AR,
  const ElementalMatrix<T>& A0,
  const ElementalMatrix<T>& A1,
  const ElementalMatrix<T>& A2 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    View( ATL, A00 );
    Merge1x2( ATR, A01, A02 );
    Merge2x1( ABL, A10, A20 );
    Merge2x2( ABR, A11, A12,
                   A21, A22 );
}

template<typename T>
void SlidePartitionUpDiagonal
( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, 
  ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, ElementalMatrix<T>& A02,
  ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, ElementalMatrix<T>& A12,
  ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, ElementalMatrix<T>& A20, 
  ElementalMatrix<T>& A21, ElementalMatrix<T>& A22 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    LockedView( ATL, A00 );
    LockedMerge1x2( ATR, A01, A02 );
    LockedMerge2x1( ABL, A10, A20 );
    LockedMerge2x2( ABR, A11, A12,
                         A21, A22 );
}

template<typename T>
void SlideLockedPartitionUpDiagonal
(       ElementalMatrix<T>& ATL,
        ElementalMatrix<T>& ATR, 
  const ElementalMatrix<T>& A00,
  const ElementalMatrix<T>& A01,
  const ElementalMatrix<T>& A02,
  const ElementalMatrix<T>& A10,
  const ElementalMatrix<T>& A11,
  const ElementalMatrix<T>& A12,
        ElementalMatrix<T>& ABL,
        ElementalMatrix<T>& ABR, 
  const ElementalMatrix<T>& A20,
  const ElementalMatrix<T>& A21,
  const ElementalMatrix<T>& A22 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    Merge2x2( ATL, A00, A01,
                   A10, A11 );
    Merge2x1( ATR, A02, A12 );
    Merge1x2( ABL, A20, A21 );
    View( ABR, A22 );
}

template<typename T>
void SlidePartitionDownDiagonal
( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, 
  ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, ElementalMatrix<T>& A02,
  ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, ElementalMatrix<T>& A12,
  ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
  ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, ElementalMatrix<T>& A22 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
    DEBUG_CSE
    LockedMerge2x2( ATL, A00, A01,
                         A10, A11 );
    LockedMerge2x1( ATR, A02, A12 );
    LockedMerge1x2( ABL, A20, A21 );
    LockedView( ABR, A22 );
}

template<typename T>
void SlideLockedPartitionDownDiagonal
(       ElementalMatrix<T>& ATL,
        ElementalMatrix<T>& ATR, 
  const ElementalMatrix<T>& A00,
  const ElementalMatrix<T>& A01,
  const ElementalMatrix<T>& A02,
  const ElementalMatrix<T>& A10,
  const ElementalMatrix<T>& A11,
  const ElementalMatrix<T>& A12,
        ElementalMatrix<T>& ABL,
        ElementalMatrix<T>& ABR, 
  const ElementalMatrix<T>& A20,
  const ElementalMatrix<T>& A21,
  const ElementalMatrix<T>& A22 )
{
    DEBUG_CSE
    DEBUG_ONLY(
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
  ( ElementalMatrix<T>& AT, ElementalMatrix<T>& A0, \
                               ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& AB, ElementalMatrix<T>& A2 ); \
  template void SlideLockedPartitionDown \
  ( Matrix<T>& AT, const Matrix<T>& A0, \
                   const Matrix<T>& A1, \
    Matrix<T>& AB, const Matrix<T>& A2 ); \
  template void SlideLockedPartitionDown \
  ( ElementalMatrix<T>& AT, const ElementalMatrix<T>& A0, \
                               const ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& AB, const ElementalMatrix<T>& A2 ); \
  /* Upward */ \
  template void SlidePartitionUp \
  ( Matrix<T>& AT, Matrix<T>& A0, \
                   Matrix<T>& A1, \
    Matrix<T>& AB, Matrix<T>& A2 ); \
  template void SlidePartitionUp \
  ( ElementalMatrix<T>& AT, ElementalMatrix<T>& A0, \
                               ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& AB, ElementalMatrix<T>& A2 ); \
  template void SlideLockedPartitionUp \
  ( Matrix<T>& AT, const Matrix<T>& A0, \
                   const Matrix<T>& A1, \
    Matrix<T>& AB, const Matrix<T>& A2 ); \
  template void SlideLockedPartitionUp \
  ( ElementalMatrix<T>& AT, const ElementalMatrix<T>& A0, \
                               const ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& AB, const ElementalMatrix<T>& A2 ); \
  /* Right */ \
  template void SlidePartitionRight \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2 ); \
  template void SlidePartitionRight \
  ( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, \
    ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, ElementalMatrix<T>& A2 ); \
  template void SlideLockedPartitionRight \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2 ); \
  template void SlideLockedPartitionRight \
  ( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, \
    const ElementalMatrix<T>& A0, const ElementalMatrix<T>& A1, const ElementalMatrix<T>& A2 ); \
  /* Left */ \
  template void SlidePartitionLeft \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2 ); \
  template void SlidePartitionLeft \
  ( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, \
    ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, ElementalMatrix<T>& A2 ); \
  template void SlideLockedPartitionLeft \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    const Matrix<T>& A0, const Matrix<T>& A1, const Matrix<T>& A2 ); \
  template void SlideLockedPartitionLeft \
  ( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, \
    const ElementalMatrix<T>& A0, const ElementalMatrix<T>& A1, const ElementalMatrix<T>& A2 ); \
  /* Down diagonal */ \
  template void SlidePartitionDownDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22 ); \
  template void SlidePartitionDownDiagonal \
  ( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, \
    ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, ElementalMatrix<T>& A02, \
    ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, ElementalMatrix<T>& A12, \
    ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, \
    ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, ElementalMatrix<T>& A22 ); \
  template void SlideLockedPartitionDownDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    const Matrix<T>& A00, const Matrix<T>& A01, const Matrix<T>& A02, \
    const Matrix<T>& A10, const Matrix<T>& A11, const Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    const Matrix<T>& A20, const Matrix<T>& A21, const Matrix<T>& A22 ); \
  template void SlideLockedPartitionDownDiagonal \
  ( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, \
    const ElementalMatrix<T>& A00, const ElementalMatrix<T>& A01, const ElementalMatrix<T>& A02, \
    const ElementalMatrix<T>& A10, const ElementalMatrix<T>& A11, const ElementalMatrix<T>& A12, \
    ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, \
    const ElementalMatrix<T>& A20, const ElementalMatrix<T>& A21, const ElementalMatrix<T>& A22 ); \
  /* Up diagonal */ \
  template void SlidePartitionUpDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22 ); \
  template void SlidePartitionUpDiagonal \
  ( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, \
    ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, ElementalMatrix<T>& A02, \
    ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, ElementalMatrix<T>& A12, \
    ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, \
    ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, ElementalMatrix<T>& A22 ); \
  template void SlideLockedPartitionUpDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    const Matrix<T>& A00, const Matrix<T>& A01, const Matrix<T>& A02, \
    const Matrix<T>& A10, const Matrix<T>& A11, const Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    const Matrix<T>& A20, const Matrix<T>& A21, const Matrix<T>& A22 ); \
  template void SlideLockedPartitionUpDiagonal \
  ( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, \
    const ElementalMatrix<T>& A00, const ElementalMatrix<T>& A01, const ElementalMatrix<T>& A02, \
    const ElementalMatrix<T>& A10, const ElementalMatrix<T>& A11, const ElementalMatrix<T>& A12, \
    ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, \
    const ElementalMatrix<T>& A20, const ElementalMatrix<T>& A21, const ElementalMatrix<T>& A22 );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
