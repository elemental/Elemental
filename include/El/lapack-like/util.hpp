/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_UTIL_HPP
#define EL_UTIL_HPP

namespace El {

// ApplyPackedReflectors
// =====================
template<typename F>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo,
  VerticalOrHorizontal dir, ForwardOrBackward order,
  Conjugation conjugation,
  Int offset, const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A );
template<typename F>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo,
  VerticalOrHorizontal dir, ForwardOrBackward order,
  Conjugation conjugation, Int offset,
  const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A );
template<typename F>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo,
  VerticalOrHorizontal dir, ForwardOrBackward order,
  Conjugation conjugation, Int offset,
  const DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t, DistMatrix<F>& A );

// ExpandPackedReflectors
// ======================
template<typename F>
void ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  Int offset, Matrix<F>& H, const Matrix<F>& t );
template<typename F>
void ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  Int offset, DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t );
template<typename F>
void ExpandPackedReflectors
( UpperOrLower uplo, VerticalOrHorizontal dir, Conjugation conjugation,
  Int offset, DistMatrix<F>& H, const DistMatrix<F,STAR,STAR>& t );

// HyperbolicReflector
// ===================
template<typename F>
F LeftHyperbolicReflector( F& chi, Matrix<F>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F LeftHyperbolicReflector( F& chi, DistMatrix<F,U,V>& x );
template<typename F>
F LeftHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F LeftHyperbolicReflector( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x );

template<typename F>
F RightHyperbolicReflector( F& chi, Matrix<F>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F RightHyperbolicReflector( F& chi, DistMatrix<F,U,V>& x );
template<typename F>
F RightHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F RightHyperbolicReflector( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x );

namespace hyp_reflector {

// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F Col( F& chi, DistMatrix<F,U,V>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F Col( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F Row( F& chi, DistMatrix<F,U,V>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F Row( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x );

} // namespace reflector

// Median
// ======
template<typename Real>
ValueInt<Real> Median( const Matrix<Real>& x );
template<typename Real,Dist U,Dist V>
ValueInt<Real> Median( const DistMatrix<Real,U,V>& x );

// Reflector
// =========
template<typename F>
F LeftReflector( F& chi, Matrix<F>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F LeftReflector( F& chi, DistMatrix<F,U,V>& x );
template<typename F>
F LeftReflector( Matrix<F>& chi, Matrix<F>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F LeftReflector( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x );

template<typename F>
F RightReflector( F& chi, Matrix<F>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F RightReflector( F& chi, DistMatrix<F,U,V>& x );
template<typename F>
F RightReflector( Matrix<F>& chi, Matrix<F>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F RightReflector( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x );

namespace reflector {

// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F Col( F& chi, DistMatrix<F,U,V>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F Col( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F Row( F& chi, DistMatrix<F,U,V>& x );
// NOTE: Only instantiated for (U,V)=(MC,MR)
template<typename F,Dist U,Dist V>
F Row( DistMatrix<F,U,V>& chi, DistMatrix<F,U,V>& x );

} // namespace reflector

// Sort
// ====
template<typename Real>
void Sort( Matrix<Real>& X, SortType sort=ASCENDING );
template<typename Real,Dist U,Dist V>
void Sort( DistMatrix<Real,U,V>& X, SortType sort=ASCENDING );

template<typename Real>
std::vector<ValueInt<Real>> TaggedSort
( const Matrix<Real>& x, SortType sort=ASCENDING );
template<typename Real,Dist U,Dist V>
std::vector<ValueInt<Real>> TaggedSort
( const DistMatrix<Real,U,V>& x, SortType sort=ASCENDING );

} // namespace El

#endif // ifndef EL_UTIL_HPP
