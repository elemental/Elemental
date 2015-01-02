/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FLAMEPART_MERGE_HPP
#define EL_FLAMEPART_MERGE_HPP

namespace El {

// Horizontally merge two contiguous matrices
// ==========================================
template<typename T>
void Merge1x2( Matrix<T>& A, Matrix<T>& BL, Matrix<T>& BR );
template<typename T>
void Merge1x2
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& BL, AbstractDistMatrix<T>& BR );

template<typename T>
void LockedMerge1x2( Matrix<T>& A, const Matrix<T>& BL, const Matrix<T>& BR );
template<typename T>
void LockedMerge1x2
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& BL, const AbstractDistMatrix<T>& BR );

// Return by value
// ^^^^^^^^^^^^^^^
template<typename T>
Matrix<T> Merge1x2( Matrix<T>& BL, Matrix<T>& BR );
template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> Merge1x2( DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR );

template<typename T>
Matrix<T> LockedMerge1x2( const Matrix<T>& BL, const Matrix<T>& BR );
template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> LockedMerge1x2
( const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR );

// Vertically merge two contiguous matrices
// ========================================
template<typename T>
void Merge2x1( Matrix<T>& A, Matrix<T>& BT, Matrix<T>& BB );
template<typename T>
void Merge2x1
( AbstractDistMatrix<T>& A, 
  AbstractDistMatrix<T>& BT, AbstractDistMatrix<T>& BB );

template<typename T>
void LockedMerge2x1( Matrix<T>& A, const Matrix<T>& BT, const Matrix<T>& BB );
template<typename T>
void LockedMerge2x1
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& BT,
  const AbstractDistMatrix<T>& BB );

// Return by value
// ^^^^^^^^^^^^^^^
template<typename T>
Matrix<T> Merge2x1( Matrix<T>& BT, Matrix<T>& BB );
template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> Merge2x1( DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB );

template<typename T>
Matrix<T> LockedMerge2x1( const Matrix<T>& BT, const Matrix<T>& BB );
template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> LockedMerge2x1
( const DistMatrix<T,U,V>& BT, const DistMatrix<T,U,V>& BB );

// Merge a contiguous 2x2 block of matrices
// ========================================
template<typename T>
void Merge2x2
( Matrix<T>& A,
  Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR );
template<typename T>
void Merge2x2
( AbstractDistMatrix<T>& A,
  AbstractDistMatrix<T>& BTL, AbstractDistMatrix<T>& BTR,
  AbstractDistMatrix<T>& BBL, AbstractDistMatrix<T>& BBR );

template<typename T>
void LockedMerge2x2
(       Matrix<T>& A,
  const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR );
template<typename T>
void LockedMerge2x2
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& BTL, const AbstractDistMatrix<T>& BTR,
  const AbstractDistMatrix<T>& BBL, const AbstractDistMatrix<T>& BBR );

// Return by value
// ^^^^^^^^^^^^^^^

template<typename T>
Matrix<T> Merge2x2
( Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR );
template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> Merge2x2
( DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR,
  DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR );

template<typename T>
Matrix<T> LockedMerge2x2
( const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR );
template<typename T,Dist U,Dist V>
DistMatrix<T,U,V> LockedMerge2x2
( const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR,
  const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR );

} // namespace El

#endif // ifndef EL_FLAMEPART_MERGE_HPP
