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
template<typename F>
void Bidiag
( Matrix<F>& A,
  Matrix<F>& phaseP,
  Matrix<F>& phaseQ );

template<typename F>
void Bidiag
( ElementalMatrix<F>& A, 
  ElementalMatrix<F>& phaseP,
  ElementalMatrix<F>& phaseQ );
template<typename F>
void Bidiag
( ElementalMatrix<F>& A, 
  ElementalMatrix<F>& phaseP,
  ElementalMatrix<F>& phaseQ );

namespace bidiag {

// Overwrite A with B = Q^H A P and additionally return P and Q
template<typename F>
void Explicit
( Matrix<F>& A, 
  Matrix<F>& P,
  Matrix<F>& Q );
template<typename F>
void Explicit
( ElementalMatrix<F>& A, 
  ElementalMatrix<F>& P,
  ElementalMatrix<F>& Q );
template<typename F>
void Explicit
( DistMatrix<F>& A, 
  DistMatrix<F>& P,
  DistMatrix<F>& Q );

// Only return the condensed bidiagonal matrix
template<typename F>
void ExplicitCondensed( Matrix<F>& A ); 
template<typename F>
void ExplicitCondensed( ElementalMatrix<F>& A );

template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A,
  const Matrix<F>& phase,
        Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, Orientation orientation,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& phase,
        ElementalMatrix<F>& B );

template<typename F>
void ApplyP
( LeftOrRight side, Orientation orientation,
  const Matrix<F>& A,
  const Matrix<F>& phase,
        Matrix<F>& B );
template<typename F>
void ApplyP
( LeftOrRight side, Orientation orientation,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& phase,
        ElementalMatrix<F>& B );

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

template<typename F>
struct HermitianTridiagCtrl 
{
    HermitianTridiagApproach approach=HERMITIAN_TRIDIAG_SQUARE;
    GridOrder order=ROW_MAJOR;
    SymvCtrl<F> symvCtrl;
};

template<typename F>
void HermitianTridiag( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& phase );
template<typename F>
void HermitianTridiag
( UpperOrLower uplo,
  ElementalMatrix<F>& A,
  ElementalMatrix<F>& phase,
  const HermitianTridiagCtrl<F>& ctrl=HermitianTridiagCtrl<F>() );

namespace herm_tridiag {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void ExplicitCondensed
( UpperOrLower uplo, ElementalMatrix<F>& A,
  const HermitianTridiagCtrl<F>& ctrl=HermitianTridiagCtrl<F>() );

template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<F>& A,
  const Matrix<F>& phase,
        Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& phase, 
        ElementalMatrix<F>& B );

} // namespace herm_tridiag

// Hessenberg
// ==========
template<typename F>
void Hessenberg( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& phase );
template<typename F>
void Hessenberg
( UpperOrLower uplo, ElementalMatrix<F>& A, ElementalMatrix<F>& phase );

namespace hessenberg {

template<typename F>
void ExplicitCondensed( UpperOrLower uplo, Matrix<F>& A );
template<typename F>
void ExplicitCondensed( UpperOrLower uplo, ElementalMatrix<F>& A );

template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const Matrix<F>& A,
  const Matrix<F>& phase,
        Matrix<F>& B );
template<typename F>
void ApplyQ
( LeftOrRight side, UpperOrLower uplo, Orientation orientation,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& phase, 
        ElementalMatrix<F>& B );

template<typename F>
void FormQ
( UpperOrLower uplo,
  const Matrix<F>& A,
  const Matrix<F>& phase,
        Matrix<F>& Q );
template<typename F>
void FormQ
( UpperOrLower uplo,
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& phase, 
        ElementalMatrix<F>& Q );

} // namespace hessenberg

} // namespace El

#endif // ifndef EL_CONDENSE_HPP
