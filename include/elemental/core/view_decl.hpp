/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef CORE_VIEW_DECL_HPP
#define CORE_VIEW_DECL_HPP

namespace elem {

//
// Viewing a full matrix
//

template<typename T,typename Int>
void View
( Matrix<T,Int>& A, Matrix<T,Int>& B );
template<typename T,Distribution U,Distribution V,typename Int>
void View
( DistMatrix<T,U,V,Int>& A, DistMatrix<T,U,V,Int>& B );

template<typename T,typename Int>
void LockedView
( Matrix<T,Int>& A, const Matrix<T,Int>& B );
template<typename T,Distribution U,Distribution V,typename Int>
void LockedView
( DistMatrix<T,U,V,Int>& A, const DistMatrix<T,U,V,Int>& B );

//
// Viewing a submatrix
//

template<typename T,typename Int>
void View
( Matrix<T,Int>& A, Matrix<T,Int>& B,
  Int i, Int j, Int height, Int width );
template<typename T,Distribution U,Distribution V,typename Int>
void View
( DistMatrix<T,U,V,Int>& A, DistMatrix<T,U,V,Int>& B,
  Int i, Int j, Int height, Int width );

template<typename T,typename Int>
void LockedView
( Matrix<T,Int>& A, const Matrix<T,Int>& B,
  Int i, Int j, Int height, Int width );
template<typename T,Distribution U,Distribution V,typename Int>
void LockedView
( DistMatrix<T,U,V,Int>& A, const DistMatrix<T,U,V,Int>& B,
  Int i, Int j, Int height, Int width );

//
// View two horizontally connected matrices
//

template<typename T,typename Int>
void View1x2
( Matrix<T,Int>& A,
  Matrix<T,Int>& BL, Matrix<T,Int>& BR );
template<typename T,Distribution U,Distribution V,typename Int>
void View1x2
( DistMatrix<T,U,V,Int>& A,
  DistMatrix<T,U,V,Int>& BL, DistMatrix<T,U,V,Int>& BR );

template<typename T,typename Int>
void LockedView1x2
(       Matrix<T,Int>& A,
  const Matrix<T,Int>& BL,
  const Matrix<T,Int>& BR );
template<typename T,Distribution U,Distribution V,typename Int>
void LockedView1x2
(       DistMatrix<T,U,V,Int>& A,
  const DistMatrix<T,U,V,Int>& BL,
  const DistMatrix<T,U,V,Int>& BR );

//
// View two vertically connected matrices
//

template<typename T,typename Int>
void View2x1
( Matrix<T,Int>& A,
  Matrix<T,Int>& BT,
  Matrix<T,Int>& BB );
template<typename T,Distribution U,Distribution V,typename Int>
void View2x1
( DistMatrix<T,U,V,Int>& A,
  DistMatrix<T,U,V,Int>& BT,
  DistMatrix<T,U,V,Int>& BB );

template<typename T,typename Int>
void LockedView2x1
(       Matrix<T,Int>& A,
  const Matrix<T,Int>& BT,
  const Matrix<T,Int>& BB );
template<typename T,Distribution U,Distribution V,typename Int>
void LockedView2x1
(       DistMatrix<T,U,V,Int>& A,
  const DistMatrix<T,U,V,Int>& BT,
  const DistMatrix<T,U,V,Int>& BB );

//
// View a two-by-two set of connected matrices
//

template<typename T,typename Int>
void View2x2
( Matrix<T,Int>& A,
  Matrix<T,Int>& BTL, Matrix<T,Int>& BTR,
  Matrix<T,Int>& BBL, Matrix<T,Int>& BBR );
template<typename T,Distribution U,Distribution V,typename Int>
void View2x2
( DistMatrix<T,U,V,Int>& A,
  DistMatrix<T,U,V,Int>& BTL, DistMatrix<T,U,V,Int>& BTR,
  DistMatrix<T,U,V,Int>& BBL, DistMatrix<T,U,V,Int>& BBR );

template<typename T,typename Int>
void LockedView2x2
(       Matrix<T,Int>& A,
  const Matrix<T,Int>& BTL,
  const Matrix<T,Int>& BTR,
  const Matrix<T,Int>& BBL,
  const Matrix<T,Int>& BBR );
template<typename T,Distribution U,Distribution V,typename Int>
void LockedView2x2
(       DistMatrix<T,U,V,Int>& A,
  const DistMatrix<T,U,V,Int>& BTL,
  const DistMatrix<T,U,V,Int>& BTR,
  const DistMatrix<T,U,V,Int>& BBL,
  const DistMatrix<T,U,V,Int>& BBR );

// Utilities for handling the extra information needed for [MD,* ] and [* ,MD]
template<typename T,Distribution U,Distribution V,typename Int>
void HandleDiagPath
( DistMatrix<T,U,V,Int>& A, const DistMatrix<T,U,V,Int>& B );
template<typename T,typename Int>
void HandleDiagPath
( DistMatrix<T,MD,STAR,Int>& A, const DistMatrix<T,MD,STAR,Int>& B );
template<typename T,typename Int>
void HandleDiagPath
( DistMatrix<T,STAR,MD,Int>& A, const DistMatrix<T,STAR,MD,Int>& B );

} // namespace elem

#endif // ifndef CORE_VIEW_DECL_HPP
