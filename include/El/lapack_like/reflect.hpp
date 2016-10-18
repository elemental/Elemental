/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_REFLECT_HPP
#define EL_REFLECT_HPP

namespace El {

// ApplyPackedReflectors
// =====================
template<typename F>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo,
  VerticalOrHorizontal dir, ForwardOrBackward order,
  Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A );
template<typename F>
void ApplyPackedReflectors
( LeftOrRight side, UpperOrLower uplo,
  VerticalOrHorizontal dir, ForwardOrBackward order,
  Conjugation conjugation, Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalars, 
        AbstractDistMatrix<F>& A );

// ExpandPackedReflectors
// ======================
template<typename F>
void ExpandPackedReflectors
( UpperOrLower uplo,
  VerticalOrHorizontal dir,
  Conjugation conjugation,
  Int offset,
        Matrix<F>& H,
  const Matrix<F>& householderScalars );
template<typename F>
void ExpandPackedReflectors
( UpperOrLower uplo,
  VerticalOrHorizontal dir,
  Conjugation conjugation,
  Int offset,
        AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalars );

// HyperbolicReflector
// ===================
template<typename F>
F LeftHyperbolicReflector( F& chi, Matrix<F>& x );
template<typename F>
F LeftHyperbolicReflector( F& chi, AbstractDistMatrix<F>& x );
template<typename F>
F LeftHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x );
template<typename F>
F LeftHyperbolicReflector
( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );

template<typename F>
F RightHyperbolicReflector( F& chi, Matrix<F>& x );
template<typename F>
F RightHyperbolicReflector( F& chi, AbstractDistMatrix<F>& x );
template<typename F>
F RightHyperbolicReflector( Matrix<F>& chi, Matrix<F>& x );
template<typename F>
F RightHyperbolicReflector
( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );

namespace hyp_reflector {

template<typename F>
F Col( F& chi, AbstractDistMatrix<F>& x );
template<typename F>
F Col( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );
template<typename F>
F Row( F& chi, AbstractDistMatrix<F>& x );
template<typename F>
F Row( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );

} // namespace reflector

// Reflector
// =========
// Since LAPACK chooses to use the identity matrix, rather than a single
// coordinate negation, in cases where the mass is already entirely in the
// first entry, and the identity matrix cannot be represented as a Householder
// reflector, Elemental does not ever directly call LAPACK's Householder
// routines. Otherwise, the logic of routines such as ApplyPackedReflectors
// would need to be (unnecessarily) complicated.
//
// Furthermore, LAPACK defines H = I - tau [1; v] [1; v]' such that
// adjoint(H) [chi; x] = [beta; 0], but Elemental instead defines
// H = I - tau [1; v] [1; v]' such that H [chi; x] = [beta; 0].
//
template<typename F>
F LeftReflector( F& chi, Matrix<F>& x );
template<typename F>
F LeftReflector( F& chi, AbstractDistMatrix<F>& x );
template<typename F>
F LeftReflector( Matrix<F>& chi, Matrix<F>& x );
template<typename F>
F LeftReflector( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );

template<typename F>
F RightReflector( F& chi, Matrix<F>& x );
template<typename F>
F RightReflector( F& chi, AbstractDistMatrix<F>& x );
template<typename F>
F RightReflector( Matrix<F>& chi, Matrix<F>& x );
template<typename F>
F RightReflector( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );

namespace reflector {

template<typename F>
F Col( F& chi, AbstractDistMatrix<F>& x );
template<typename F>
F Col( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );
template<typename F>
F Row( F& chi, AbstractDistMatrix<F>& x );
template<typename F>
F Row( AbstractDistMatrix<F>& chi, AbstractDistMatrix<F>& x );

} // namespace reflector

} // namespace El

#endif // ifndef EL_REFLECT_HPP
