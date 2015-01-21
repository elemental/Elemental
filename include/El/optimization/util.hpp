/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_OPTIMIZATION_UTIL_HPP
#define EL_OPTIMIZATION_UTIL_HPP

namespace El {

// Coherence
// =========
template<typename F>
Base<F> Coherence( const Matrix<F>& A );
template<typename F>
Base<F> Coherence( const AbstractDistMatrix<F>& A );

// Covariance
// ==========
template<typename F>
void Covariance( const Matrix<F>& D, Matrix<F>& S );
template<typename F>
void Covariance( const AbstractDistMatrix<F>& D, AbstractDistMatrix<F>& S );

// Log barrier
// ===========
template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, const Matrix<F>& A );
template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

template<typename F>
Base<F> LogBarrier
( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false );
template<typename F>
Base<F> LogBarrier
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool canOverwrite=false );

// Log-det divergence
// ==================
template<typename F>
Base<F> LogDetDiv
( UpperOrLower uplo, const Matrix<F>& A, const Matrix<F>& B );
template<typename F>
Base<F> LogDetDiv
( UpperOrLower uplo, 
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B );

// Maximum step within the positive cone
// =====================================
template<typename Real>
Real MaxStepInPositiveCone
( const Matrix<Real>& s, 
  const Matrix<Real>& ds, Real upperBound );
template<typename Real>
Real MaxStepInPositiveCone
( const AbstractDistMatrix<Real>& s, 
  const AbstractDistMatrix<Real>& ds, Real upperBound );
template<typename Real>
Real MaxStepInPositiveCone
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& ds, Real upperBound );

// Number of non-positive entries
// ==============================
template<typename Real>
Int NumNonPositive( const Matrix<Real>& A );
template<typename Real>
Int NumNonPositive( const SparseMatrix<Real>& A );
template<typename Real>
Int NumNonPositive( const AbstractDistMatrix<Real>& A );
template<typename Real>
Int NumNonPositive( const DistSparseMatrix<Real>& A );
template<typename Real>
Int NumNonPositive( const DistMultiVec<Real>& A );

// Number of non-SOC members
// =========================
template<typename Real>
Int NumNonSecondOrder
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
Int NumNonSecondOrder
( const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
Int NumNonSecondOrder
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_UTIL_HPP
