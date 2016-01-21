/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_UTIL_SOC_HPP
#define EL_OPTIMIZATION_UTIL_SOC_HPP

namespace El {
namespace soc {

// TODO: SOC eigenvectors?

// Apply an SOC vector as a linear operator
// ========================================
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Apply
( const Matrix<Real>& x, 
  const Matrix<Real>& y, 
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Apply
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Real>& y, 
        ElementalMatrix<Real>& z,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
        Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Apply
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
        Int cutoff=1000 );

// Overwrite y with x o y
// ----------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Apply
( const Matrix<Real>& x, 
        Matrix<Real>& y,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Apply
( const ElementalMatrix<Real>& x, 
        ElementalMatrix<Real>& y,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Apply
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Apply the quadratic representation of a product of SOCs to a vector
// ===================================================================
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ApplyQuadratic
( const Matrix<Real>& x, 
  const Matrix<Real>& y, 
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ApplyQuadratic
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Real>& y, 
        ElementalMatrix<Real>& z,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ApplyQuadratic
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Overwrite y with Q_x y
// ----------------------
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ApplyQuadratic
( const Matrix<Real>& x, 
        Matrix<Real>& y,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ApplyQuadratic
( const ElementalMatrix<Real>& x, 
        ElementalMatrix<Real>& y,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void ApplyQuadratic
( const DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Degree
// ======
Int Degree( const Matrix<Int>& firstInds );
Int Degree( const ElementalMatrix<Int>& firstInds );
Int Degree( const DistMultiVec<Int>& firstInds );

// Determinants
// ============
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Dets
( const Matrix<Real>& x,
        Matrix<Real>& d,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Dets
( const ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& d,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Dets
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& d,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Dot products of sequences of second-order cones
// ===============================================
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Dots
( const Matrix<Real>& x, 
  const Matrix<Real>& y, 
        Matrix<Real>& z,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Dots
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Real>& y, 
        ElementalMatrix<Real>& z,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Dots
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Embedding maps
// ==============
void EmbeddingMaps
( const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Int>& sparseOrders,
        Matrix<Int>& sparseFirstInds,
        Matrix<Int>& origToSparseOrders,
        Matrix<Int>& origToSparseFirstInds,
        Matrix<Int>& sparseToOrigOrders,
        Matrix<Int>& sparseToOrigFirstInds,
  Int cutoffSparse );
void EmbeddingMaps
( const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Int>& sparseOrders,
        DistMultiVec<Int>& sparseFirstInds,
        DistMultiVec<Int>& origToSparseOrders,
        DistMultiVec<Int>& origToSparseFirstInds,
        DistMultiVec<Int>& sparseToOrigOrders,
        DistMultiVec<Int>& sparseToOrigFirstInds,
  Int cutoffSparse );

// Identity
// ========
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Identity
(       Matrix<Real>& e, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Identity
(       ElementalMatrix<Real>& e, 
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Identity
(       DistMultiVec<Real>& e, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds );

// Compute the inverse in the product SOC Jordan algebra
// =====================================================
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Inverse
( const Matrix<Real>& x,
        Matrix<Real>& xInv,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Inverse
( const ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& xInv,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Inverse
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& xInv,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Lower norms
// ===========
template<typename Real,typename=EnableIf<IsReal<Real>>>
void LowerNorms
( const Matrix<Real>& x,
        Matrix<Real>& lowerNorms,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void LowerNorms
( const ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& lowerNorms,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void LowerNorms
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& lowerNorms,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Max eigenvalues
// ===============
template<typename Real,typename=EnableIf<IsReal<Real>>>
void MaxEig
( const Matrix<Real>& x,
        Matrix<Real>& maxEigs,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void MaxEig
( const ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& maxEigs,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void MaxEig
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& maxEigs,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MaxEig
( const Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MaxEig
( const ElementalMatrix<Real>& x,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MaxEig
( const DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Maximum step in a product of second-order cones
// ===============================================
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MaxStep
( const Matrix<Real>& x, 
  const Matrix<Real>& y,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds,
  Real upperBound=std::numeric_limits<Real>::max() );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MaxStep
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Real>& y,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
  Real upperBound=std::numeric_limits<Real>::max(),
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MaxStep
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& y,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Real upperBound=std::numeric_limits<Real>::max(),
  Int cutoff=1000 );

// Min eigenvalues
// ===============
template<typename Real,typename=EnableIf<IsReal<Real>>>
void MinEig
( const Matrix<Real>& x,
        Matrix<Real>& minEigs,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void MinEig
( const ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& minEigs,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void MinEig
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& minEigs,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MinEig
( const Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MinEig
( const ElementalMatrix<Real>& x,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real MinEig
( const DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Compute an SOC Nesterov-Todd point
// ==================================
// The Nesterov-Todd point, w, is a member of the SOC whose quadratic 
// representation maps z to s, where s and z are both members of the SOC.
template<typename Real,typename=EnableIf<IsReal<Real>>>
void NesterovTodd
( const Matrix<Real>& s, 
  const Matrix<Real>& z, 
        Matrix<Real>& w,
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void NesterovTodd
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Real>& z, 
        ElementalMatrix<Real>& w,
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void NesterovTodd
( const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z, 
        DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Number of non-SOC members
// =========================
// Return the number of negative determinants
template<typename Real,typename=EnableIf<IsReal<Real>>>
Int NumOutside
( const Matrix<Real>& x, 
  const Matrix<Int>& orders, 
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Int NumOutside
( const ElementalMatrix<Real>& x, 
  const ElementalMatrix<Int>& orders, 
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
Int NumOutside
( const DistMultiVec<Real>& x, 
  const DistMultiVec<Int>& orders, 
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

// Push into SOC
// ==============
template<typename Real,typename=EnableIf<IsReal<Real>>>
void PushInto
(       Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  Real minDist=0 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void PushInto
(       ElementalMatrix<Real>& x,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Real minDist=0,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void PushInto
(       DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Real minDist=0,
  Int cutoff=1000 );

// Push pair into SOC
// ==================
template<typename Real,typename=EnableIf<IsReal<Real>>>
void PushPairInto
(       Matrix<Real>& s,
        Matrix<Real>& z,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  Real wMaxNormLimit );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void PushPairInto
(       ElementalMatrix<Real>& s,
        ElementalMatrix<Real>& z,
  const ElementalMatrix<Real>& w,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Real wMaxNormLimit,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void PushPairInto
(       DistMultiVec<Real>& s,
        DistMultiVec<Real>& z,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Real wMaxNormLimit,
  Int cutoff=1000 );

// Reflect
// =======
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Reflect
(       Matrix<Real>& x,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Reflect
(       ElementalMatrix<Real>& x,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Reflect
(       DistMultiVec<Real>& x,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds );

// Shift
// =====
// Add a multiple of the identity of the product cone
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Shift
(       Matrix<Real>& x,
        Real shift,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Shift
(       ElementalMatrix<Real>& x,
        Real shift,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void Shift
(       DistMultiVec<Real>& x,
        Real shift,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds );

// Compute the square-root in the product SOC Jordan algebra
// =========================================================
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SquareRoot
( const Matrix<Real>& x,
        Matrix<Real>& xRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SquareRoot
( const ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& xRoot,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
  Int cutoff=1000 );
template<typename Real,typename=EnableIf<IsReal<Real>>>
void SquareRoot
( const DistMultiVec<Real>& x,
        DistMultiVec<Real>& xRoot,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  Int cutoff=1000 );

} // namespace soc
} // namespace El

#endif // ifndef EL_OPTIMIZATION_UTIL_SOC_HPP
