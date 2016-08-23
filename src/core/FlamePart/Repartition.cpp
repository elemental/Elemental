/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "El/core/FlamePart.hpp"

namespace El {

// Repartition upwards from the bottom
// ===================================

template<typename T>
void RepartitionUp
( Matrix<T>& AT, Matrix<T>& A0,
                 Matrix<T>& A1,
  Matrix<T>& AB, Matrix<T>& A2, Int A1Height )
{
    DEBUG_CSE
    PartitionUp( AT, A0, A1, A1Height );
    View( A2, AB );
}

template<typename T>
void RepartitionUp
( ElementalMatrix<T>& AT, ElementalMatrix<T>& A0,
                          ElementalMatrix<T>& A1,
  ElementalMatrix<T>& AB, ElementalMatrix<T>& A2, Int A1Height )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( AT, AB, A0, A1, A2 );
      AssertSameDists( AT, AB, A0, A1, A2 );
    )
    PartitionUp( AT, A0, A1, A1Height );
    View( A2, AB );
}

template<typename T>
void LockedRepartitionUp
( const Matrix<T>& AT, Matrix<T>& A0,
                       Matrix<T>& A1,
  const Matrix<T>& AB, Matrix<T>& A2, Int A1Height )
{
    DEBUG_CSE
    LockedPartitionUp( AT, A0, A1, A1Height );
    LockedView( A2, AB );
}

template<typename T>
void LockedRepartitionUp
( const ElementalMatrix<T>& AT, ElementalMatrix<T>& A0,
                                ElementalMatrix<T>& A1,
  const ElementalMatrix<T>& AB, ElementalMatrix<T>& A2, Int A1Height )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( AT, AB, A0, A1, A2 );
      AssertSameDists( AT, AB, A0, A1, A2 );
    )
    LockedPartitionUp( AT, A0, A1, A1Height );
    LockedView( A2, AB );
}

// Repartition downwards from the top
// ==================================

template<typename T>
void RepartitionDown
( Matrix<T>& AT, Matrix<T>& A0,
                 Matrix<T>& A1,
  Matrix<T>& AB, Matrix<T>& A2, Int A1Height )
{
    DEBUG_CSE
    View( A0, AT );
    PartitionDown( AB, A1, A2, A1Height );
}

template<typename T>
void RepartitionDown
( ElementalMatrix<T>& AT, ElementalMatrix<T>& A0,
                          ElementalMatrix<T>& A1,
  ElementalMatrix<T>& AB, ElementalMatrix<T>& A2, Int A1Height )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( AT, AB, A0, A1, A2 );
      AssertSameDists( AT, AB, A0, A1, A2 );
    )
    View( A0, AT );
    PartitionDown( AB, A1, A2, A1Height );
}

template<typename T>
void LockedRepartitionDown
( const Matrix<T>& AT, Matrix<T>& A0,
                       Matrix<T>& A1,
  const Matrix<T>& AB, Matrix<T>& A2, Int A1Height )
{
    DEBUG_CSE
    LockedView( A0, AT );
    LockedPartitionDown( AB, A1, A2, A1Height );
}

template<typename T>
void LockedRepartitionDown
( const ElementalMatrix<T>& AT, ElementalMatrix<T>& A0,
                                ElementalMatrix<T>& A1,
  const ElementalMatrix<T>& AB, ElementalMatrix<T>& A2, Int A1Height )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( AT, AB, A0, A1, A2 );
      AssertSameDists( AT, AB, A0, A1, A2 );
    )
    LockedView( A0, AT );
    LockedPartitionDown( AB, A1, A2, A1Height );
}

// Repartition leftwards from the right
// ====================================

template<typename T>
void RepartitionLeft
( Matrix<T>& AL, Matrix<T>& AR,
  Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int A1Width )
{
    DEBUG_CSE
    PartitionLeft( AL, A0, A1, A1Width );
    View( A2, AR );
}

template<typename T>
void RepartitionLeft
( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR,
  ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, 
  ElementalMatrix<T>& A2, Int A1Width )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( AL, AR, A0, A1, A2 );
      AssertSameDists( AL, AR, A0, A1, A2 );
    )
    PartitionLeft( AL, A0, A1, A1Width );
    View( A2, AR );
}

template<typename T>
void LockedRepartitionLeft
( const Matrix<T>& AL, const Matrix<T>& AR,
  Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int A1Width )
{
    DEBUG_CSE
    LockedPartitionLeft( AL, A0, A1, A1Width );
    LockedView( A2, AR );
}

template<typename T>
void LockedRepartitionLeft
( const ElementalMatrix<T>& AL, const ElementalMatrix<T>& AR,
  ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, 
  ElementalMatrix<T>& A2, Int A1Width )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( AL, AR, A0, A1, A2 );
      AssertSameDists( AL, AR, A0, A1, A2 );
    )
    LockedPartitionLeft( AL, A0, A1, A1Width );
    LockedView( A2, AR );
}

// Repartition rightward from the left
// ===================================

template<typename T>
void RepartitionRight
( Matrix<T>& AL, Matrix<T>& AR,
  Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int A1Width )
{
    DEBUG_CSE
    View( A0, AL );
    PartitionRight( AR, A1, A2, A1Width );
}

template<typename T>
void RepartitionRight
( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR,
  ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, 
  ElementalMatrix<T>& A2, Int A1Width )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( AL, AR, A0, A1, A2 );
      AssertSameDists( AL, AR, A0, A1, A2 );
    )
    View( A0, AL );
    PartitionRight( AR, A1, A2, A1Width );
}

template<typename T>
void LockedRepartitionRight
( const Matrix<T>& AL, const Matrix<T>& AR,
  Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int A1Width )
{
    DEBUG_CSE
    LockedView( A0, AL );
    LockedPartitionRight( AR, A1, A2, A1Width );
}

template<typename T>
void LockedRepartitionRight
( const ElementalMatrix<T>& AL, const ElementalMatrix<T>& AR,
  ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, 
  ElementalMatrix<T>& A2, Int A1Width )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( AL, AR, A0, A1, A2 );
      AssertSameDists( AL, AR, A0, A1, A2 );
    )
    LockedView( A0, AL );
    LockedPartitionRight( AR, A1, A2, A1Width );
}

// Repartition upwards on a diagonal
// =================================

template<typename T>
void RepartitionUpDiagonal
( Matrix<T>& ATL, Matrix<T>& ATR, 
  Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02,
  Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12,
  Matrix<T>& ABL, Matrix<T>& ABR, 
  Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize )
{
    DEBUG_CSE
    PartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    PartitionUp( ATR, A02, A12, A11.Height() );
    PartitionLeft( ABL, A20, A21, A11.Width() );
    View( A22, ABR );
}

template<typename T>
void RepartitionUpDiagonal
( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, 
  ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, 
  ElementalMatrix<T>& A02,
  ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, 
  ElementalMatrix<T>& A12,
  ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
  ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, 
  ElementalMatrix<T>& A22, Int bsize )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
      AssertSameDists
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
    )
    PartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    PartitionUp( ATR, A02, A12, A11.Height() );
    PartitionLeft( ABL, A20, A21, A11.Width() );
    View( A22, ABR );
}

template<typename T>
void LockedRepartitionUpDiagonal
( const Matrix<T>& ATL, const Matrix<T>& ATR, 
  Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02,
  Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12,
  const Matrix<T>& ABL, const Matrix<T>& ABR, 
  Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize )
{
    DEBUG_CSE
    LockedPartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    LockedPartitionUp( ATR, A02, A12, A11.Height() );
    LockedPartitionLeft( ABL, A20, A21, A11.Width() );
    LockedView( A22, ABR );
}

template<typename T>
void LockedRepartitionUpDiagonal
( const ElementalMatrix<T>& ATL, const ElementalMatrix<T>& ATR, 
  ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, 
  ElementalMatrix<T>& A02,
  ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, 
  ElementalMatrix<T>& A12,
  const ElementalMatrix<T>& ABL, const ElementalMatrix<T>& ABR, 
  ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, 
  ElementalMatrix<T>& A22, Int bsize )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
      AssertSameDists
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
    )
    LockedPartitionUpOffsetDiagonal
    ( ATL.Width()-ATL.Height(),
      ATL, A00, A01,
           A10, A11, bsize );
    LockedPartitionUp( ATR, A02, A12, A11.Height() );
    LockedPartitionLeft( ABL, A20, A21, A11.Width() );
    LockedView( A22, ABR );
}

// Repartition downwards on a diagonal
// ===================================

template<typename T>
void RepartitionDownDiagonal
( Matrix<T>& ATL, Matrix<T>& ATR, 
  Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02,
  Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12,
  Matrix<T>& ABL, Matrix<T>& ABR, 
  Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize )
{
    DEBUG_CSE
    View( A00, ATL );
    PartitionDownDiagonal( ABR, A11, A12,
                                A21, A22, bsize );
    PartitionDown( ABL, A10, A20, A11.Height() );
    PartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T>
void RepartitionDownDiagonal
( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, 
  ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, 
  ElementalMatrix<T>& A02,
  ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, 
  ElementalMatrix<T>& A12,
  ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
  ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, 
  ElementalMatrix<T>& A22, Int bsize )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
      AssertSameDists
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
    )
    View( A00, ATL );
    PartitionDownDiagonal( ABR, A11, A12,
                                A21, A22, bsize );
    PartitionDown( ABL, A10, A20, A11.Height() );
    PartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T>
void LockedRepartitionDownDiagonal
( const Matrix<T>& ATL, const Matrix<T>& ATR, 
  Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02,
  Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12,
  const Matrix<T>& ABL, const Matrix<T>& ABR, 
  Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize )
{
    DEBUG_CSE
    LockedView( A00, ATL );
    LockedPartitionDownDiagonal( ABR, A11, A12,
                                      A21, A22, bsize );
    LockedPartitionDown( ABL, A10, A20, A11.Height() );
    LockedPartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T>
void LockedRepartitionDownDiagonal
( const ElementalMatrix<T>& ATL, const ElementalMatrix<T>& ATR, 
  ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, 
  ElementalMatrix<T>& A02,
  ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, 
  ElementalMatrix<T>& A12,
  const ElementalMatrix<T>& ABL, const ElementalMatrix<T>& ABR, 
  ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, 
  ElementalMatrix<T>& A22, Int bsize )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
      AssertSameDists
      ( ATL, ATR, ABL, ABR, A00, A01, A02, A10, A11, A12, A20, A21, A22 );
    )
    LockedView( A00, ATL );
    LockedPartitionDownDiagonal( ABR, A11, A12,
                                      A21, A22, bsize );
    LockedPartitionDown( ABL, A10, A20, A11.Height() );
    LockedPartitionRight( ATR, A01, A02, A11.Width() );
}

#define PROTO(T) \
  /* Downward */ \
  template void RepartitionDown \
  ( Matrix<T>& AT, Matrix<T>& A0, \
                   Matrix<T>& A1, \
    Matrix<T>& AB, Matrix<T>& A2, Int bsize ); \
  template void RepartitionDown \
  ( ElementalMatrix<T>& AT, ElementalMatrix<T>& A0, \
                            ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& AB, ElementalMatrix<T>& A2, Int bsize ); \
  template void LockedRepartitionDown \
  ( const Matrix<T>& AT, Matrix<T>& A0, \
                         Matrix<T>& A1, \
    const Matrix<T>& AB, Matrix<T>& A2, Int bsize ); \
  template void LockedRepartitionDown \
  ( const ElementalMatrix<T>& AT, ElementalMatrix<T>& A0, \
                                  ElementalMatrix<T>& A1, \
    const ElementalMatrix<T>& AB, ElementalMatrix<T>& A2, Int bsize ); \
  /* Upward */ \
  template void RepartitionUp \
  ( Matrix<T>& AT, Matrix<T>& A0, \
                   Matrix<T>& A1, \
    Matrix<T>& AB, Matrix<T>& A2, Int bsize ); \
  template void RepartitionUp \
  ( ElementalMatrix<T>& AT, ElementalMatrix<T>& A0, \
                            ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& AB, ElementalMatrix<T>& A2, Int bsize ); \
  template void LockedRepartitionUp \
  ( const Matrix<T>& AT, Matrix<T>& A0, \
                         Matrix<T>& A1, \
    const Matrix<T>& AB, Matrix<T>& A2, Int bsize ); \
  template void LockedRepartitionUp \
  ( const ElementalMatrix<T>& AT, ElementalMatrix<T>& A0, \
                                  ElementalMatrix<T>& A1, \
    const ElementalMatrix<T>& AB, ElementalMatrix<T>& A2, Int bsize ); \
  /* Rightward */ \
  template void RepartitionRight \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int bsize ); \
  template void RepartitionRight \
  ( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, \
    ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& A2, Int bsize ); \
  template void LockedRepartitionRight \
  ( const Matrix<T>& AL, const Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int bsize ); \
  template void LockedRepartitionRight \
  ( const ElementalMatrix<T>& AL, const ElementalMatrix<T>& AR, \
    ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& A2, Int bsize ); \
  /* Leftward */ \
  template void RepartitionLeft \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int bsize ); \
  template void RepartitionLeft \
  ( ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, \
    ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& A2, Int bsize ); \
  template void LockedRepartitionLeft \
  ( const Matrix<T>& AL, const Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int bsize ); \
  template void LockedRepartitionLeft \
  ( const ElementalMatrix<T>& AL, const ElementalMatrix<T>& AR, \
    ElementalMatrix<T>& A0, ElementalMatrix<T>& A1, \
    ElementalMatrix<T>& A2, Int bsize ); \
  /* Down a diagonal */ \
  template void RepartitionDownDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize ); \
  template void RepartitionDownDiagonal \
  ( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, \
    ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, \
    ElementalMatrix<T>& A02, \
    ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, \
    ElementalMatrix<T>& A12, \
    ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, \
    ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, \
    ElementalMatrix<T>& A22, Int bsize ); \
  template void LockedRepartitionDownDiagonal \
  ( const Matrix<T>& ATL, const Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    const Matrix<T>& ABL, const Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize ); \
  template void LockedRepartitionDownDiagonal \
  ( const ElementalMatrix<T>& ATL, const ElementalMatrix<T>& ATR, \
    ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, \
    ElementalMatrix<T>& A02, \
    ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, \
    ElementalMatrix<T>& A12, \
    const ElementalMatrix<T>& ABL, const ElementalMatrix<T>& ABR, \
    ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, \
    ElementalMatrix<T>& A22, Int bsize ); \
  /* Up a diagonal */ \
  template void RepartitionUpDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize ); \
  template void RepartitionUpDiagonal \
  ( ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR, \
    ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, \
    ElementalMatrix<T>& A02, \
    ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, \
    ElementalMatrix<T>& A12, \
    ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, \
    ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, \
    ElementalMatrix<T>& A22, Int bsize ); \
  template void LockedRepartitionUpDiagonal \
  ( const Matrix<T>& ATL, const Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    const Matrix<T>& ABL, const Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize ); \
  template void LockedRepartitionUpDiagonal \
  ( const ElementalMatrix<T>& ATL, const ElementalMatrix<T>& ATR, \
    ElementalMatrix<T>& A00, ElementalMatrix<T>& A01, \
    ElementalMatrix<T>& A02, \
    ElementalMatrix<T>& A10, ElementalMatrix<T>& A11, \
    ElementalMatrix<T>& A12, \
    const ElementalMatrix<T>& ABL, const ElementalMatrix<T>& ABR, \
    ElementalMatrix<T>& A20, ElementalMatrix<T>& A21, \
    ElementalMatrix<T>& A22, Int bsize );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
