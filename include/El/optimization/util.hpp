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
  const Matrix<Real>& ds, 
  Real upperBound=std::numeric_limits<Real>::max() );
template<typename Real>
Real MaxStepInPositiveCone
( const AbstractDistMatrix<Real>& s, 
  const AbstractDistMatrix<Real>& ds, 
  Real upperBound=std::numeric_limits<Real>::max() );
template<typename Real>
Real MaxStepInPositiveCone
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& ds, 
  Real upperBound=std::numeric_limits<Real>::max() );

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

// SOC Identity
// ============
template<typename Real>
void SOCIdentity
(       Matrix<Real>& e, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCIdentity
(       AbstractDistMatrix<Real>& e, 
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds );
template<typename Real>
void SOCIdentity
(       DistMultiVec<Real>& e, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds );

// Dot products of sequences of second-order cones
// ===============================================
template<typename Real>
void SOCDots
( const Matrix<Real>& x, 
  const Matrix<Real>& y, 
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCDots
( const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Real>& y, 
        AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCDots
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// SOC Broadcast
// =============
// Replicate the entry in the root position in each SOC over the entire cone
template<typename Real>
void SOCBroadcast
(       Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCBroadcast
(       AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCBroadcast
(       DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// SOC Reflect
// ===========
template<typename Real>
void SOCReflect
(       Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCReflect
(       AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds );
template<typename Real>
void SOCReflect
(       DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds );

// SOC Determinants
// ================
template<typename Real>
void SOCDets
( const Matrix<Real>& x,
        Matrix<Real>& d,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCDets
( const AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCDets
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& d,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// Number of non-SOC members
// =========================
// Return the number of negative determinants
template<typename Real>
Int NumNonSOC
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
Int NumNonSOC
( const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
Int NumNonSOC
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// Apply an SOC vector as a linear operator
// ========================================
template<typename Real>
void SOCApply
( const Matrix<Real>& x, 
  const Matrix<Real>& y, 
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCApply
( const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Real>& y, 
        AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCApply
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// Apply the quadratic representation of a product of SOCs to a vector
// ===================================================================
template<typename Real>
void SOCApplyQuadratic
( const Matrix<Real>& x, 
  const Matrix<Real>& y, 
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCApplyQuadratic
( const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Real>& y, 
        AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCApplyQuadratic
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// Compute the inverse in the product SOC Jordan algebra
// =====================================================
template<typename Real>
void SOCInverse
( const Matrix<Real>& x,
        Matrix<Real>& xInv,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCInverse
( const AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& xInv,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCInverse
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& xInv,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// Compute the square-root in the product SOC Jordan algebra
// =========================================================
template<typename Real>
void SOCSquareRoot
( const Matrix<Real>& x,
        Matrix<Real>& xRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCSquareRoot
( const AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& xRoot,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCSquareRoot
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& xRoot,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// Compute an SOC Nesterov-Todd point
// ==================================
// The Nesterov-Todd point, w, is a member of the SOC whose quadratic 
// representation maps s to z, where s and z are both members of the SOC.
template<typename Real>
void SOCNesterovTodd
( const Matrix<Real>& s, 
  const Matrix<Real>& z, 
        Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCNesterovTodd
( const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Real>& z, 
        AbstractDistMatrix<Real>& w,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCNesterovTodd
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z, 
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// Maximum step in a product of second-order cones
// ===============================================
template<typename Real>
Real MaxStepInSOC
( const Matrix<Real>& x, 
  const Matrix<Real>& y,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  Real upperBound=std::numeric_limits<Real>::max() );
template<typename Real>
Real MaxStepInSOC
( const AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds,
  Real upperBound=std::numeric_limits<Real>::max(),
  Int cutoff=1000 );
template<typename Real>
Real MaxStepInSOC
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Real upperBound=std::numeric_limits<Real>::max(),
  Int cutoff=1000 );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_UTIL_HPP
