/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CORE_VIEW_DECL_HPP
#define ELEM_CORE_VIEW_DECL_HPP

namespace elem {

//
// Viewing a full matrix
//

template<typename T>
void View( Matrix<T>& A, Matrix<T>& B );
template<typename T>
Matrix<T> View( Matrix<T>& B );
template<typename T,Distribution U,Distribution V>
void View( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> View( DistMatrix<T,U,V>& B );

template<typename T>
void LockedView( Matrix<T>& A, const Matrix<T>& B );
template<typename T>
Matrix<T> LockedView( const Matrix<T>& B );
template<typename T,Distribution U,Distribution V>
void LockedView( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> LockedView( const DistMatrix<T,U,V>& B );

//
// Viewing a submatrix
//

template<typename T>
void View( Matrix<T>& A, Matrix<T>& B, Int i, Int j, Int height, Int width );
template<typename T>
Matrix<T> View( Matrix<T>& B, Int i, Int j, Int height, Int width );
template<typename T,Distribution U,Distribution V>
void View
( DistMatrix<T,U,V>& A, DistMatrix<T,U,V>& B,
  Int i, Int j, Int height, Int width );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> View
( DistMatrix<T,U,V>& B, Int i, Int j, Int height, Int width );

template<typename T>
void LockedView
( Matrix<T>& A, const Matrix<T>& B,
  Int i, Int j, Int height, Int width );
template<typename T>
Matrix<T> LockedView( const Matrix<T>& B, Int i, Int j, Int height, Int width );
template<typename T,Distribution U,Distribution V>
void LockedView
( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B,
  Int i, Int j, Int height, Int width );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> LockedView
( const DistMatrix<T,U,V>& B, Int i, Int j, Int height, Int width );

//
// View two horizontally connected matrices
//

template<typename T>
void View1x2
( Matrix<T>& A,
  Matrix<T>& BL, Matrix<T>& BR );
template<typename T>
Matrix<T> View1x2( Matrix<T>& BL, Matrix<T>& BR );
template<typename T,Distribution U,Distribution V>
void View1x2
( DistMatrix<T,U,V>& A,
  DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> View1x2( DistMatrix<T,U,V>& BL, DistMatrix<T,U,V>& BR );

template<typename T>
void LockedView1x2
(       Matrix<T>& A,
  const Matrix<T>& BL,
  const Matrix<T>& BR );
template<typename T>
Matrix<T> LockedView1x2( const Matrix<T>& BL, const Matrix<T>& BR );
template<typename T,Distribution U,Distribution V>
void LockedView1x2
(       DistMatrix<T,U,V>& A,
  const DistMatrix<T,U,V>& BL,
  const DistMatrix<T,U,V>& BR );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> LockedView1x2
( const DistMatrix<T,U,V>& BL, const DistMatrix<T,U,V>& BR );

//
// View two vertically connected matrices
//

template<typename T>
void View2x1
( Matrix<T>& A,
  Matrix<T>& BT,
  Matrix<T>& BB );
template<typename T>
Matrix<T> View2x1( Matrix<T>& BT, Matrix<T>& BB );
template<typename T,Distribution U,Distribution V>
void View2x1
( DistMatrix<T,U,V>& A,
  DistMatrix<T,U,V>& BT,
  DistMatrix<T,U,V>& BB );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> View2x1( DistMatrix<T,U,V>& BT, DistMatrix<T,U,V>& BB );

template<typename T>
void LockedView2x1
(       Matrix<T>& A,
  const Matrix<T>& BT,
  const Matrix<T>& BB );
template<typename T>
Matrix<T> LockedView2x1( const Matrix<T>& BT, const Matrix<T>& BB );
template<typename T,Distribution U,Distribution V>
void LockedView2x1
(       DistMatrix<T,U,V>& A,
  const DistMatrix<T,U,V>& BT,
  const DistMatrix<T,U,V>& BB );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> LockedView2x1
( const DistMatrix<T,U,V>& BT, const DistMatrix<T,U,V>& BB );

//
// View a two-by-two set of connected matrices
//

template<typename T>
void View2x2
( Matrix<T>& A,
  Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR );
template<typename T>
Matrix<T> View2x2
( Matrix<T>& BTL, Matrix<T>& BTR,
  Matrix<T>& BBL, Matrix<T>& BBR );
template<typename T,Distribution U,Distribution V>
void View2x2
( DistMatrix<T,U,V>& A,
  DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR,
  DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> View2x2
( DistMatrix<T,U,V>& BTL, DistMatrix<T,U,V>& BTR,
  DistMatrix<T,U,V>& BBL, DistMatrix<T,U,V>& BBR );

template<typename T>
void LockedView2x2
(       Matrix<T>& A,
  const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR );
template<typename T>
Matrix<T> LockedView2x2
( const Matrix<T>& BTL, const Matrix<T>& BTR,
  const Matrix<T>& BBL, const Matrix<T>& BBR );
template<typename T,Distribution U,Distribution V>
void LockedView2x2
(       DistMatrix<T,U,V>& A,
  const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR,
  const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR );
template<typename T,Distribution U,Distribution V>
DistMatrix<T,U,V> LockedView2x2
( const DistMatrix<T,U,V>& BTL, const DistMatrix<T,U,V>& BTR,
  const DistMatrix<T,U,V>& BBL, const DistMatrix<T,U,V>& BBR );

// Utilities for handling the extra information needed for [MD,* ] and [* ,MD]
template<typename T,Distribution U,Distribution V>
void HandleDiagPath
( DistMatrix<T,U,V>& A, const DistMatrix<T,U,V>& B );
template<typename T>
void HandleDiagPath
( DistMatrix<T,MD,STAR>& A, const DistMatrix<T,MD,STAR>& B );
template<typename T>
void HandleDiagPath
( DistMatrix<T,STAR,MD>& A, const DistMatrix<T,STAR,MD>& B );

} // namespace elem

#endif // ifndef ELEM_CORE_VIEW_DECL_HPP
