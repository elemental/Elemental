/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
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
( ElementalMatrix<T>& A, 
  ElementalMatrix<T>& AT, ElementalMatrix<T>& AB, 
  Int heightAT=Blocksize() );

template<typename T>
void LockedPartitionDown
( const Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAT=Blocksize() );

template<typename T>
void LockedPartitionDown
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& AT, ElementalMatrix<T>& AB, 
  Int heightAT=Blocksize() );

// Partition upwards from the bottom
// =================================

template<typename T>
void PartitionUp
( Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAB=Blocksize() );

template<typename T>
void PartitionUp
( ElementalMatrix<T>& A, 
  ElementalMatrix<T>& AT, ElementalMatrix<T>& AB, 
  Int heightAB=Blocksize() );

template<typename T>
void LockedPartitionUp
( const Matrix<T>& A, Matrix<T>& AT, Matrix<T>& AB, Int heightAB=Blocksize() );

template<typename T>
void LockedPartitionUp
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& AT, ElementalMatrix<T>& AB, 
  Int heightAB=Blocksize() );

// Partition rightwards from the left
// ==================================

template<typename T>
void PartitionRight
( Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAL=Blocksize() );

template<typename T>
void PartitionRight
( ElementalMatrix<T>& A, 
  ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, 
  Int widthAL=Blocksize() );

template<typename T>
void LockedPartitionRight
( const Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAL=Blocksize() );

template<typename T>
void LockedPartitionRight
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, 
  Int widthAL=Blocksize() );

// Partition leftwards from the right
// ==================================

template<typename T>
void PartitionLeft
( Matrix<T>& A, Matrix<T>& AL, Matrix<T>& AR, Int widthAR=Blocksize() );

template<typename T>
void PartitionLeft
( ElementalMatrix<T>& A, 
  ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, 
  Int widthAR=Blocksize() );

template<typename T>
void LockedPartitionLeft
( const Matrix<T>& A, 
        Matrix<T>& AL, Matrix<T>& AR, Int widthAR=Blocksize() );

template<typename T>
void LockedPartitionLeft
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& AL, ElementalMatrix<T>& AR, 
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
  ElementalMatrix<T>& A, 
  ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR,
  ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
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
  const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR,
        ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
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
  ElementalMatrix<T>& A, 
  ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR,
  ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
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
  const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR,
        ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
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
( ElementalMatrix<T>& A, 
  ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR,
  ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionDownDiagonal
( const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionDownDiagonal
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR,
        ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
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
( ElementalMatrix<T>& A, 
  ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR,
  ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionUpDiagonal
( const Matrix<T>& A, 
        Matrix<T>& ATL, Matrix<T>& ATR,
        Matrix<T>& ABL, Matrix<T>& ABR, Int diagDist=Blocksize() );

template<typename T>
void LockedPartitionUpDiagonal
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& ATL, ElementalMatrix<T>& ATR,
        ElementalMatrix<T>& ABL, ElementalMatrix<T>& ABR, 
  Int diagDist=Blocksize() );

} // namespace El

#endif // ifndef EL_FLAMEPART_PARTITION_HPP
