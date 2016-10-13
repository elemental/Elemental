/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_VIEW_DECL_HPP
#define EL_VIEW_DECL_HPP

namespace El {

// View an entire matrix
// =====================

// (Sequential) matrix
// -------------------

template<typename T>
void View( Matrix<T>& A, Matrix<T>& B );
template<typename T>
void LockedView( Matrix<T>& A, const Matrix<T>& B );

template<typename T>
Matrix<T> View( Matrix<T>& B );
template<typename T>
Matrix<T> LockedView( const Matrix<T>& B );

// ElementalMatrix
// ---------------

template<typename T>
void View( ElementalMatrix<T>& A, ElementalMatrix<T>& B );
template<typename T>
void LockedView( ElementalMatrix<T>& A, const ElementalMatrix<T>& B );

// Return by value
// ^^^^^^^^^^^^^^^
template<typename T,Dist U,Dist V,DistWrap wrapType>
DistMatrix<T,U,V,wrapType> View( DistMatrix<T,U,V,wrapType>& B );
template<typename T,Dist U,Dist V,DistWrap wrapType>
DistMatrix<T,U,V,wrapType> LockedView( const DistMatrix<T,U,V,wrapType>& B );

// BlockMatrix
// -----------
template<typename T>
void View( BlockMatrix<T>& A, BlockMatrix<T>& B );
template<typename T>
void LockedView( BlockMatrix<T>& A, const BlockMatrix<T>& B );

// Mixed
// -----
template<typename T>
void View
( BlockMatrix<T>& A,
  ElementalMatrix<T>& B );
template<typename T>
void LockedView
(       BlockMatrix<T>& A,
  const ElementalMatrix<T>& B );

template<typename T>
void View
( ElementalMatrix<T>& A,
  BlockMatrix<T>& B );
template<typename T>
void LockedView
(       ElementalMatrix<T>& A,
  const BlockMatrix<T>& B );

// AbstractDistMatrix
// ------------------
template<typename T>
void View
( AbstractDistMatrix<T>& A,
  AbstractDistMatrix<T>& B );
template<typename T>
void LockedView
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B );

// View a contiguous submatrix
// ===========================

// (Sequential) Matrix
// -------------------

template<typename T>
void View
( Matrix<T>& A,
  Matrix<T>& B,
  Int i, Int j,
  Int height, Int width );
template<typename T>
void LockedView
(       Matrix<T>& A,
  const Matrix<T>& B,
  Int i, Int j,
  Int height, Int width );

template<typename T>
void View
( Matrix<T>& A,
  Matrix<T>& B,
  Range<Int> I, Range<Int> J );
template<typename T>
void LockedView
(       Matrix<T>& A,
  const Matrix<T>& B,
  Range<Int> I, Range<Int> J );

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T>
Matrix<T> View( Matrix<T>& B, Int i, Int j, Int height, Int width );
template<typename T>
Matrix<T> LockedView( const Matrix<T>& B, Int i, Int j, Int height, Int width );

template<typename T>
Matrix<T> View( Matrix<T>& B, Range<Int> I, Range<Int> J );
template<typename T>
Matrix<T> LockedView( const Matrix<T>& B, Range<Int> I, Range<Int> J );

// ElementalMatrix
// ---------------

template<typename T>
void View
( ElementalMatrix<T>& A,
  ElementalMatrix<T>& B,
  Int i, Int j, Int height, Int width );
template<typename T>
void LockedView
(       ElementalMatrix<T>& A,
  const ElementalMatrix<T>& B,
  Int i, Int j, Int height, Int width );

template<typename T>
void View
( ElementalMatrix<T>& A,
  ElementalMatrix<T>& B, 
  Range<Int> I, Range<Int> J );
template<typename T>
void LockedView
(       ElementalMatrix<T>& A,
  const ElementalMatrix<T>& B, 
  Range<Int> I, Range<Int> J );

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T,Dist U,Dist V,DistWrap wrapType>
DistMatrix<T,U,V,wrapType> View
( DistMatrix<T,U,V,wrapType>& B, Int i, Int j, Int height, Int width );
template<typename T,Dist U,Dist V,DistWrap wrapType>
DistMatrix<T,U,V,wrapType> LockedView
( const DistMatrix<T,U,V,wrapType>& B, Int i, Int j, Int height, Int width );

template<typename T,Dist U,Dist V,DistWrap wrapType>
DistMatrix<T,U,V,wrapType> View
( DistMatrix<T,U,V,wrapType>& B, Range<Int> I, Range<Int> J );
template<typename T,Dist U,Dist V,DistWrap wrapType>
DistMatrix<T,U,V,wrapType> LockedView
( const DistMatrix<T,U,V,wrapType>& B, Range<Int> I, Range<Int> J );

// BlockMatrix
// -----------

template<typename T>
void View
( BlockMatrix<T>& A,
  BlockMatrix<T>& B,
  Int i,
  Int j,
  Int height,
  Int width );
template<typename T>
void LockedView
(       BlockMatrix<T>& A,
  const BlockMatrix<T>& B,
  Int i, Int j, Int height, Int width );

template<typename T>
void View
( BlockMatrix<T>& A,
  BlockMatrix<T>& B, 
  Range<Int> I, Range<Int> J );
template<typename T>
void LockedView
(       BlockMatrix<T>& A,
  const BlockMatrix<T>& B, 
  Range<Int> I, Range<Int> J );

// AbstractDistMatrix
// ------------------
template<typename T>
void View
( AbstractDistMatrix<T>& A,
  AbstractDistMatrix<T>& B,
  Int i, Int j, Int height, Int width );
template<typename T>
void LockedView
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
  Int i, Int j, Int height, Int width );

template<typename T>
void View
( AbstractDistMatrix<T>& A,
  AbstractDistMatrix<T>& B, 
  Range<Int> I, Range<Int> J );
template<typename T>
void LockedView
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B, 
  Range<Int> I, Range<Int> J );

} // namespace El

#endif // ifndef EL_VIEW_DECL_HPP
