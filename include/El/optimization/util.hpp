/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_UTIL_HPP
#define EL_OPTIMIZATION_UTIL_HPP

#include "./util/cone.hpp"
#include "./util/pos_orth.hpp"
#include "./util/soc.hpp"

namespace El {

// Coherence
// =========
template<typename F>
Base<F> Coherence( const Matrix<F>& A );
template<typename F>
Base<F> Coherence( const ElementalMatrix<F>& A );

// Covariance
// ==========
template<typename F>
void Covariance( const Matrix<F>& D, Matrix<F>& S );
template<typename F>
void Covariance( const ElementalMatrix<F>& D, ElementalMatrix<F>& S );

// Log barrier
// ===========
template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, const ElementalMatrix<F>& A );

template<typename F>
Base<F> LogBarrier
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false );
template<typename F>
Base<F> LogBarrier
( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite=false );

// Log-det divergence
// ==================
template<typename F>
Base<F> LogDetDiv
( UpperOrLower uplo,
  const Matrix<F>& A,
  const Matrix<F>& B );
template<typename F>
Base<F> LogDetDiv
( UpperOrLower uplo, 
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& B );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_UTIL_HPP
