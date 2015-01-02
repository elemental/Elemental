/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FLAMEPART_PARTITION_HPP
#define EL_FLAMEPART_PARTITION_HPP

namespace El {

// Partition downwards from the top
// ================================

template<typename T>
void PartitionDown
( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAT=Blocksize() );

template<typename T>
void PartitionDown
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, 
  Int heightAT=Blocksize() );

template<typename T>
void LockedPartitionDown
( const Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAT=Blocksize() );

template<typename T>
void LockedPartitionDown
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, 
  Int heightAT=Blocksize() );

// Partition upwards from the bottom
// =================================

template<typename T>
void PartitionUp
( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAB=Blocksize() );

template<typename T>
void PartitionUp
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, 
  Int heightAB=Blocksize() );

template<typename T>
void LockedPartitionUp
( const Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAB=Blocksize() );

template<typename T>
void LockedPartitionUp
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& AT, AbstractDistMatrix<T>& AB, 
  Int heightAB=Blocksize() );

// Partition rightwards from the left
// ==================================

template<typename T>
void PartitionRight
( Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAL=Blocksize() );

template<typename T>
void PartitionRight
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, 
  Int widthAL=Blocksize() );

template<typename T>
void LockedPartitionRight
( const Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAL=Blocksize() );

template<typename T>
void LockedPartitionRight
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, 
  Int widthAL=Blocksize() );

// Partition leftwards from the right
// ==================================

template<typename T>
void PartitionLeft
( Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAR=Blocksize() );

template<typename T>
void PartitionLeft
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, 
  Int widthAR=Blocksize() );

template<typename T>
void LockedPartitionLeft
( const Matrix<T>& A, 
        Matrix<T>& AL, Matrix<T>& AR, Int widthAR=Blocksize() );

template<typename T>
void LockedPartitionLeft
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& AL, AbstractDistMatrix<T>& AR, 
  Int widthAR=Blocksize() );

// Partition downward on a particular diagonal
// ===========================================

template<typename T>
void PartitionDownOffsetDiagonal
( Int offset,
  Matrix<T>& A, 
  Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void PartitionDownOffsetDiagonal
( Int offset,
  AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionDownOffsetDiagonal
( Int offset,
  const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionDownOffsetDiagonal
( Int offset,
  const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
        AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

// Partition upwards on a particular diagonal
// ==========================================

template<typename T>
void PartitionUpOffsetDiagonal
( Int offset,
  Matrix<T>& A, 
  Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void PartitionUpOffsetDiagonal
( Int offset,
  AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionUpOffsetDiagonal
( Int offset,
  const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionUpOffsetDiagonal
( Int offset,
  const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
        AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

// Partition downwards on the main diagonal
// ========================================

template<typename T>
void PartitionDownDiagonal
( Matrix<T>& A, 
  Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void PartitionDownDiagonal
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionDownDiagonal
( const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionDownDiagonal
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
        AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

// Partition upwards on the main diagonal
// ======================================

template<typename T>
void PartitionUpDiagonal
( Matrix<T>& A, 
  Matrix<T>& ATL, Matrix<T>& ATR,
  Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void PartitionUpDiagonal
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
  AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionUpDiagonal
( const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionUpDiagonal
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& ATL, AbstractDistMatrix<T>& ATR,
        AbstractDistMatrix<T>& ABL, AbstractDistMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

} // namespace El

#endif // ifndef EL_FLAMEPART_PARTITION_HPP
