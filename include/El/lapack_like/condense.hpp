/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CONDENSE_HPP
#define EL_CONDENSE_HPP

#include <El/blas_like/level2.hpp>

namespace El {

// Bidiag
// ======

// Return the packed reduction to bidiagonal form
// ----------------------------------------------
template<typename Field>
void Bidiag
( Matrix<Field>& A,
  Matrix<Field>& householderScalarsP,
  Matrix<Field>& householderScalarsQ );

template<typename Field>
void Bidiag
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalarsP,
  AbstractDistMatrix<Field>& householderScalarsQ );
template<typename Field>
void Bidiag
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalarsP,
  AbstractDistMatrix<Field>& householderScalarsQ );

namespace bidiag {

// Overwrite A with B = Q^H A P and additionally return P and Q
template<typename Field>
void Explicit
( Matrix<Field>& A,
  Matrix<Field>& P,
  Matrix<Field>& Q );
template<typename Field>
void Explicit
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& P,
  AbstractDistMatrix<Field>& Q );
template<typename Field>
void Explicit
( DistMatrix<Field>& A,
  DistMatrix<Field>& P,
  DistMatrix<Field>& Q );

// Only return the condensed bidiagonal matrix
template<typename Field>
void ExplicitCondensed( Matrix<Field>& A );
template<typename Field>
void ExplicitCondensed( AbstractDistMatrix<Field>& A );

template<typename Field>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
        Matrix<Field>& B );
template<typename Field>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
        AbstractDistMatrix<Field>& B );

template<typename Field>
void ApplyP
( LeftOrRight side, Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
        Matrix<Field>& B );
template<typename Field>
void ApplyP
( LeftOrRight side, Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
        AbstractDistMatrix<Field>& B );

} // namespace bidiag

// HermitianTridiag
// ================

namespace HermitianTridiagApproachNS {
enum HermitianTridiagApproach
{
    HERMITIAN_TRIDIAG_NORMAL, // Keep the current grid
    HERMITIAN_TRIDIAG_SQUARE, // Drop to a square process grid
    HERMITIAN_TRIDIAG_DEFAULT // Square grid algorithm only if already square
};
}
using namespace HermitianTridiagApproachNS;

template<typename Field>
struct HermitianTridiagCtrl
{
    HermitianTridiagApproach approach=HERMITIAN_TRIDIAG_SQUARE;
    GridOrder order=ROW_MAJOR;
    SymvCtrl<Field> symvCtrl;
};

template<typename Field>
void HermitianTridiag
( UpperOrLower uplo, Matrix<Field>& A, Matrix<Field>& householderScalars );
template<typename Field>
void HermitianTridiag
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalars,
  const HermitianTridiagCtrl<Field>& ctrl=HermitianTridiagCtrl<Field>() );

namespace herm_tridiag {

template<typename Field>
void ExplicitCondensed( UpperOrLower uplo, Matrix<Field>& A );
template<typename Field>
void ExplicitCondensed
( UpperOrLower uplo, AbstractDistMatrix<Field>& A,
  const HermitianTridiagCtrl<Field>& ctrl=HermitianTridiagCtrl<Field>() );

template<typename Field>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
        Matrix<Field>& B );
template<typename Field>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
        AbstractDistMatrix<Field>& B );

} // namespace herm_tridiag

// Hessenberg
// ==========
template<typename Field>
void Hessenberg
( UpperOrLower uplo, Matrix<Field>& A, Matrix<Field>& householderScalars );
template<typename Field>
void Hessenberg
( UpperOrLower uplo, AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& householderScalars );

namespace hessenberg {

template<typename Field>
void ExplicitCondensed( UpperOrLower uplo, Matrix<Field>& A );
template<typename Field>
void ExplicitCondensed( UpperOrLower uplo, AbstractDistMatrix<Field>& A );

template<typename Field>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
        Matrix<Field>& B );
template<typename Field>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
        AbstractDistMatrix<Field>& B );

template<typename Field>
void FormQ
( UpperOrLower uplo,
  const Matrix<Field>& A,
  const Matrix<Field>& householderScalars,
        Matrix<Field>& Q );
template<typename Field>
void FormQ
( UpperOrLower uplo,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& householderScalars,
        AbstractDistMatrix<Field>& Q );

} // namespace hessenberg

} // namespace El

#endif // ifndef EL_CONDENSE_HPP
