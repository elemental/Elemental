/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

// Repartition upwards from the bottom
// ===================================

template<typename T>
void RepartitionUp
( Matrix<T>& AT, Matrix<T>& A0,
                 Matrix<T>& A1,
  Matrix<T>& AB, Matrix<T>& A2, Int A1Height )
{
    DEBUG_ONLY(CallStackEntry cse("RepartitionUp"))
    PartitionUp( AT, A0, A1, A1Height );
    View( A2, AB );
}

template<typename T>
void RepartitionUp
( AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0,
                             AbstractDistMatrix<T>& A1,
  AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2, Int A1Height )
{
    DEBUG_ONLY(
      CallStackEntry cse("RepartitionUp");
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
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionUp"))
    LockedPartitionUp( AT, A0, A1, A1Height );
    LockedView( A2, AB );
}

template<typename T>
void LockedRepartitionUp
( const AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0,
                                   AbstractDistMatrix<T>& A1,
  const AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2, Int A1Height )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedRepartitionUp");
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
    DEBUG_ONLY(CallStackEntry cse("RepartitionDown"))
    View( A0, AT );
    PartitionDown( AB, A1, A2, A1Height );
}

template<typename T>
void RepartitionDown
( AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0,
                             AbstractDistMatrix<T>& A1,
  AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2, Int A1Height )
{
    DEBUG_ONLY(
      CallStackEntry cse("RepartitionDown");
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
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionDown"))
    LockedView( A0, AT );
    LockedPartitionDown( AB, A1, A2, A1Height );
}

template<typename T>
void LockedRepartitionDown
( const AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0,
                                   AbstractDistMatrix<T>& A1,
  const AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2, Int A1Height )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedRepartitionDown");
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
    DEBUG_ONLY(CallStackEntry cse("RepartitionLeft"))
    PartitionLeft( AL, A0, A1, A1Width );
    View( A2, AR );
}

template<typename T>
void RepartitionLeft
( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR,
  AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, 
  AbstractDistMatrix<T>& A2, Int A1Width )
{
    DEBUG_ONLY(
      CallStackEntry cse("RepartitionLeft");
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
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionLeft"))
    LockedPartitionLeft( AL, A0, A1, A1Width );
    LockedView( A2, AR );
}

template<typename T>
void LockedRepartitionLeft
( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR,
  AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, 
  AbstractDistMatrix<T>& A2, Int A1Width )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedRepartitionLeft");
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
    DEBUG_ONLY(CallStackEntry cse("RepartitionRight"))
    View( A0, AL );
    PartitionRight( AR, A1, A2, A1Width );
}

template<typename T>
void RepartitionRight
( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR,
  AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, 
  AbstractDistMatrix<T>& A2, Int A1Width )
{
    DEBUG_ONLY(
      CallStackEntry cse("RepartitionRight");
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
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionRight"))
    LockedView( A0, AL );
    LockedPartitionRight( AR, A1, A2, A1Width );
}

template<typename T>
void LockedRepartitionRight
( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR,
  AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, 
  AbstractDistMatrix<T>& A2, Int A1Width )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedRepartitionRight");
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
    DEBUG_ONLY(CallStackEntry cse("RepartitionUpDiagonal"))
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
( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, 
  AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, 
  AbstractDistMatrix<T>& A02,
  AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, 
  AbstractDistMatrix<T>& A12,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, 
  AbstractDistMatrix<T>& A22, Int bsize )
{
    DEBUG_ONLY(
      CallStackEntry cse("RepartitionUpDiagonal");
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
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionUpDiagonal"))
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
( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR, 
  AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, 
  AbstractDistMatrix<T>& A02,
  AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, 
  AbstractDistMatrix<T>& A12,
  const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR, 
  AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, 
  AbstractDistMatrix<T>& A22, Int bsize )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedRepartitionUpDiagonal");
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
    DEBUG_ONLY(CallStackEntry cse("RepartitionDownDiagonal"))
    View( A00, ATL );
    PartitionDownDiagonal( ABR, A11, A12,
                                A21, A22, bsize );
    PartitionDown( ABL, A10, A20, A11.Height() );
    PartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T>
void RepartitionDownDiagonal
( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, 
  AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, 
  AbstractDistMatrix<T>& A02,
  AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, 
  AbstractDistMatrix<T>& A12,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, 
  AbstractDistMatrix<T>& A22, Int bsize )
{
    DEBUG_ONLY(
      CallStackEntry cse("RepartitionDownDiagonal");
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
    DEBUG_ONLY(CallStackEntry cse("LockedRepartitionDownDiagonal"))
    LockedView( A00, ATL );
    LockedPartitionDownDiagonal( ABR, A11, A12,
                                      A21, A22, bsize );
    LockedPartitionDown( ABL, A10, A20, A11.Height() );
    LockedPartitionRight( ATR, A01, A02, A11.Width() );
}

template<typename T>
void LockedRepartitionDownDiagonal
( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR, 
  AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, 
  AbstractDistMatrix<T>& A02,
  AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, 
  AbstractDistMatrix<T>& A12,
  const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR, 
  AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, 
  AbstractDistMatrix<T>& A22, Int bsize )
{
    DEBUG_ONLY(
      CallStackEntry cse("LockedRepartitionDownDiagonal");
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
  ( AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0, \
                               AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2, Int bsize ); \
  template void LockedRepartitionDown \
  ( const Matrix<T>& AT, Matrix<T>& A0, \
                         Matrix<T>& A1, \
    const Matrix<T>& AB, Matrix<T>& A2, Int bsize ); \
  template void LockedRepartitionDown \
  ( const AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0, \
                                     AbstractDistMatrix<T>& A1, \
    const AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2, Int bsize ); \
  /* Upward */ \
  template void RepartitionUp \
  ( Matrix<T>& AT, Matrix<T>& A0, \
                   Matrix<T>& A1, \
    Matrix<T>& AB, Matrix<T>& A2, Int bsize ); \
  template void RepartitionUp \
  ( AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0, \
                               AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2, Int bsize ); \
  template void LockedRepartitionUp \
  ( const Matrix<T>& AT, Matrix<T>& A0, \
                         Matrix<T>& A1, \
    const Matrix<T>& AB, Matrix<T>& A2, Int bsize ); \
  template void LockedRepartitionUp \
  ( const AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& A0, \
                                     AbstractDistMatrix<T>& A1, \
    const AbstractDistMatrix<T>& AB, AbstractDistMatrix<T>& A2, Int bsize ); \
  /* Rightward */ \
  template void RepartitionRight \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int bsize ); \
  template void RepartitionRight \
  ( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, \
    AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& A2, Int bsize ); \
  template void LockedRepartitionRight \
  ( const Matrix<T>& AL, const Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int bsize ); \
  template void LockedRepartitionRight \
  ( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR, \
    AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& A2, Int bsize ); \
  /* Leftward */ \
  template void RepartitionLeft \
  ( Matrix<T>& AL, Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int bsize ); \
  template void RepartitionLeft \
  ( AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, \
    AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& A2, Int bsize ); \
  template void LockedRepartitionLeft \
  ( const Matrix<T>& AL, const Matrix<T>& AR, \
    Matrix<T>& A0, Matrix<T>& A1, Matrix<T>& A2, Int bsize ); \
  template void LockedRepartitionLeft \
  ( const AbstractDistMatrix<T>& AL, const AbstractDistMatrix<T>& AR, \
    AbstractDistMatrix<T>& A0, AbstractDistMatrix<T>& A1, \
    AbstractDistMatrix<T>& A2, Int bsize ); \
  /* Down a diagonal */ \
  template void RepartitionDownDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize ); \
  template void RepartitionDownDiagonal \
  ( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, \
    AbstractDistMatrix<T>& A02, \
    AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, \
    AbstractDistMatrix<T>& A12, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, \
    AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, \
    AbstractDistMatrix<T>& A22, Int bsize ); \
  template void LockedRepartitionDownDiagonal \
  ( const Matrix<T>& ATL, const Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    const Matrix<T>& ABL, const Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize ); \
  template void LockedRepartitionDownDiagonal \
  ( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, \
    AbstractDistMatrix<T>& A02, \
    AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, \
    AbstractDistMatrix<T>& A12, \
    const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR, \
    AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, \
    AbstractDistMatrix<T>& A22, Int bsize ); \
  /* Up a diagonal */ \
  template void RepartitionUpDiagonal \
  ( Matrix<T>& ATL, Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    Matrix<T>& ABL, Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize ); \
  template void RepartitionUpDiagonal \
  ( AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, \
    AbstractDistMatrix<T>& A02, \
    AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, \
    AbstractDistMatrix<T>& A12, \
    AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, \
    AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, \
    AbstractDistMatrix<T>& A22, Int bsize ); \
  template void LockedRepartitionUpDiagonal \
  ( const Matrix<T>& ATL, const Matrix<T>& ATR, \
    Matrix<T>& A00, Matrix<T>& A01, Matrix<T>& A02, \
    Matrix<T>& A10, Matrix<T>& A11, Matrix<T>& A12, \
    const Matrix<T>& ABL, const Matrix<T>& ABR, \
    Matrix<T>& A20, Matrix<T>& A21, Matrix<T>& A22, Int bsize ); \
  template void LockedRepartitionUpDiagonal \
  ( const AbstractDistMatrix<T>& ATL, const AbstractDistMatrix<T>& ATR, \
    AbstractDistMatrix<T>& A00, AbstractDistMatrix<T>& A01, \
    AbstractDistMatrix<T>& A02, \
    AbstractDistMatrix<T>& A10, AbstractDistMatrix<T>& A11, \
    AbstractDistMatrix<T>& A12, \
    const AbstractDistMatrix<T>& ABL, const AbstractDistMatrix<T>& ABR, \
    AbstractDistMatrix<T>& A20, AbstractDistMatrix<T>& A21, \
    AbstractDistMatrix<T>& A22, Int bsize );

#include "El/macros/Instantiate.h"

} // namespace El
