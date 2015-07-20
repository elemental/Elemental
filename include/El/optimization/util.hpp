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

// Compute a positive-orthant Nesterov-Todd point
// ==============================================
// The Nesterov-Todd point, w, is a member of the positive orthant whose 
// quadratic representation maps z to s.
template<typename Real>
void PositiveNesterovTodd
( const Matrix<Real>& s, 
  const Matrix<Real>& z, 
        Matrix<Real>& w );
template<typename Real>
void PositiveNesterovTodd
( const AbstractDistMatrix<Real>& s, 
  const AbstractDistMatrix<Real>& z, 
        AbstractDistMatrix<Real>& w );
template<typename Real>
void PositiveNesterovTodd
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z, 
        DistMultiVec<Real>& w );

// Cone Broadcast
// ==============
// Replicate the entry in the root position in each cone over the entire cone
template<typename Real>
void ConeBroadcast
(       Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void ConeBroadcast
(       AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void ConeBroadcast
(       DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// Cone AllReduce
// ==============
// Fill each subcone with the reduction over each cone
template<typename Real>
void ConeAllReduce
(       Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  mpi::Op op=mpi::SUM );
template<typename Real>
void ConeAllReduce
(       AbstractDistMatrix<Real>& x, 
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, 
  mpi::Op op=mpi::SUM, Int cutoff=1000 );
template<typename Real>
void ConeAllReduce
(       DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, 
  mpi::Op op=mpi::SUM, Int cutoff=1000 );

// A specialization of Ruiz scaling which respects a product of cones
// ==================================================================
template<typename F>
void ConeRuizEquil
(       Matrix<F>& A,
        Matrix<F>& B,
        Matrix<Base<F>>& dRowA,
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename F>
void ConeRuizEquil
(       AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& B,
        AbstractDistMatrix<Base<F>>& dRowA,
        AbstractDistMatrix<Base<F>>& dRowB,
        AbstractDistMatrix<Base<F>>& dCol,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

template<typename F>
void ConeRuizEquil
(       SparseMatrix<F>& A,
        SparseMatrix<F>& B,
        Matrix<Base<F>>& dRowA,
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename F>
void ConeRuizEquil
(       DistSparseMatrix<F>& A,
        DistSparseMatrix<F>& B,
        DistMultiVec<Base<F>>& dRowA,
        DistMultiVec<Base<F>>& dRowB,
        DistMultiVec<Base<F>>& dCol,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

// A specialization of GeomEquil which respects a product of cones
// ===============================================================
template<typename F>
void ConeGeomEquil
(       Matrix<F>& A,
        Matrix<F>& B,
        Matrix<Base<F>>& dRowA,
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename F>
void ConeGeomEquil
(       AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<F>& BPre,
        AbstractDistMatrix<Base<F>>& dRowAPre,
        AbstractDistMatrix<Base<F>>& dRowBPre,
        AbstractDistMatrix<Base<F>>& dColPre,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

template<typename F>
void ConeGeomEquil
(       SparseMatrix<F>& A,
        SparseMatrix<F>& B,
        Matrix<Base<F>>& dRowA,
        Matrix<Base<F>>& dRowB,
        Matrix<Base<F>>& dCol,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool progress=false );

template<typename F>
void ConeGeomEquil
(       DistSparseMatrix<F>& A,
        DistSparseMatrix<F>& B,
        DistMultiVec<Base<F>>& dRowA,
        DistMultiVec<Base<F>>& dRowB,
        DistMultiVec<Base<F>>& dCol,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000, bool progress=false );

// SOC Degree
// ==========
Int SOCDegree( const Matrix<Int>& firstInds );
Int SOCDegree( const AbstractDistMatrix<Int>& firstInds );
Int SOCDegree( const DistMultiVec<Int>& firstInds );

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

// SOC shift
// =========
// Add a multiple of the identity of the product cone
template<typename Real>
void SOCShift
(       Matrix<Real>& x, Real shift,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCShift
(       AbstractDistMatrix<Real>& x, Real shift,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds );
template<typename Real>
void SOCShift
(       DistMultiVec<Real>& x, Real shift,
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

// SOC lower norms
// ===============
template<typename Real>
void SOCLowerNorms
( const Matrix<Real>& x,
        Matrix<Real>& lowerNorms,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCLowerNorms
( const AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& lowerNorms,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCLowerNorms
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& lowerNorms,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// SOC min eigenvalues
// ===================
template<typename Real>
void SOCMinEig
( const Matrix<Real>& x,
        Matrix<Real>& minEigs,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCMinEig
( const AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& minEigs,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCMinEig
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& minEigs,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

template<typename Real>
Real SOCMinEig
( const Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
Real SOCMinEig
( const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
Real SOCMinEig
( const DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// SOC max eigenvalues
// ===================
template<typename Real>
void SOCMaxEig
( const Matrix<Real>& x,
        Matrix<Real>& maxEigs,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCMaxEig
( const AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& maxEigs,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCMaxEig
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& maxEigs,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

template<typename Real>
Real SOCMaxEig
( const Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real>
Real SOCMaxEig
( const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
Real SOCMaxEig
( const DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds, Int cutoff=1000 );

// TODO: SOC eigenvectors?

// Force into SOC
// ==============
template<typename Real>
void ForceIntoSOC
(       Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  Real minDist=0 );
template<typename Real>
void ForceIntoSOC
(       AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Real minDist=0, Int cutoff=1000 );
template<typename Real>
void ForceIntoSOC
(       DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Real minDist=0, Int cutoff=1000 );

// Force pair into SOC
// ===================
template<typename Real>
void ForcePairIntoSOC
(       Matrix<Real>& s,
        Matrix<Real>& z,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  Real wMaxNormLimit );
template<typename Real>
void ForcePairIntoSOC
(       AbstractDistMatrix<Real>& s,
        AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& w,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  Real wMaxNormLimit, Int cutoff=1000 );
template<typename Real>
void ForcePairIntoSOC
(       DistMultiVec<Real>& s,
        DistMultiVec<Real>& z,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Real wMaxNormLimit, Int cutoff=1000 );

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

// Overwrite y with x o y
// ----------------------
template<typename Real>
void SOCApply
( const Matrix<Real>& x, 
        Matrix<Real>& y,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCApply
( const AbstractDistMatrix<Real>& x, 
        AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCApply
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
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

// Overwrite y with Q_x y
// ----------------------
template<typename Real>
void SOCApplyQuadratic
( const Matrix<Real>& x, 
        Matrix<Real>& y,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real>
void SOCApplyQuadratic
( const AbstractDistMatrix<Real>& x, 
        AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Int>& orders, 
  const AbstractDistMatrix<Int>& firstInds, Int cutoff=1000 );
template<typename Real>
void SOCApplyQuadratic
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
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
// representation maps z to s, where s and z are both members of the SOC.
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
