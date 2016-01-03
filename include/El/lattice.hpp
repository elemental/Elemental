/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LATTICE_HPP
#define EL_LATTICE_HPP

namespace El {

// Deep insertion
// ==============
template<typename F>
void DeepColSwap( Matrix<F>& B, Int i, Int k );
template<typename F>
void DeepRowSwap( Matrix<F>& B, Int i, Int k );

// Lenstra-Lenstra-Lovasz (LLL) lattice reduction
// ==============================================
// A reduced basis, say D, is an LLL(delta) reduction of an m x n matrix B if
//
//    B U = D = Q R,
//
// where U is unimodular (integer-valued with absolute determinant of 1)
// and Q R is a floating-point QR factorization of D that satisfies the three
//  properties:
//
//   1. R has non-negative diagonal
//
//   2. R is (eta) size-reduced:
//
//        | R(i,j) / R(i,i) | < phi(F) eta,  for all i < j, and
//
//      where phi(F) is 1 for a real field F or sqrt(2) for a complex
//      field F, and
//
//   3. R is (delta) Lovasz reduced:
//
//        delta R(i,i)^2 <= R(i+1,i+1)^2 + |R(i,i+1)|^2,  for all i.
//
// Please see
//
//   Henri Cohen, "A course in computational algebraic number theory"
// 
// for more information on the "MLLL" variant of LLL used by Elemental to 
// handle linearly dependent vectors (the algorithm was originally suggested by
// Mike Pohst).
//

template<typename Real>
struct LLLInfo
{
    Real delta;
    Real eta; 
    Int rank;
    Int nullity; 
    Int numSwaps;
    Real logVol;
};

// Return the Gaussian estimate of the minimum-length vector
// 
//   GH(L) = (1/sqrt(pi)) Gamma(n/2+1)^{1/n} |det(L)|^{1/n}.
//
// where n is the rank of the lattice L.
template<typename Real,typename=EnableIf<IsReal<Real>>>
Real LatticeGaussianHeuristic( Int n, Real logVol )
{
    return Exp((LogGamma(Real(n)/Real(2)+Real(1))+logVol)/Real(n))/
           Sqrt(Pi<Real>());
}

enum LLLVariant {
  // A 'weak' LLL reduction only ensures that | R(i,i+1) / R(i,i) | is
  // bounded above by eta (or, for complex data, by sqrt(2) eta), but it often
  // produces much lower-quality basis vectors
  LLL_WEAK,
  LLL_NORMAL,
  // LLL with 'deep insertion' is no longer guaranteed to be polynomial time
  // but produces significantly higher quality bases than normal LLL.
  // See Schnorr and Euchner's "Lattice Basis Reduction: Improved Practical
  // Algorithms and Solving Subset Sum Problems".
  LLL_DEEP,
  // Going one step further, one can perform additional size reduction before
  // checking each deep insertion condition. See Schnorr's article
  // "Progress on LLL and Lattice Reduction" in the book "The LLL Algorithm",
  // edited by Nguyen and Vallee.
  LLL_DEEP_REDUCE
};

template<typename Real>
struct LLLCtrl
{
    Real delta=Real(3)/Real(4);
    Real eta=Real(1)/Real(2) + Pow(limits::Epsilon<Real>(),Real(0.9));

    LLLVariant variant=LLL_NORMAL;

    // Preprocessing with a "rank-obscuring" column-pivoted QR factorization
    // (in the manner suggested by Wubben et al.) can greatly decrease
    // the number of swaps within LLL in some circumstances
    bool presort=false;
    bool smallestFirst=true;

    // If the size-reduced column has a two-norm that is less than or
    // equal to `reorthogTol` times the  original two-norm, then reorthog.
    Real reorthogTol=0;

    // The number of times to execute the orthogonalization
    Int numOrthog=1;

    // If a size-reduced column has a two-norm less than or equal to 'zeroTol',
    // then it is interpreted as a zero vector (and forced to zero)
    Real zeroTol=Pow(limits::Epsilon<Real>(),Real(0.9));

    bool progress=false;
    bool time=false;

    // If 'jumpstart' is true, start LLL under the assumption that the first
    // 'startCol' columns are already processed
    bool jumpstart=false;
    Int startCol=0;
};

// TODO: Maintain B in BigInt form

template<typename F>
LLLInfo<Base<F>> LLL
( Matrix<F>& B,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo<Base<F>> LLL
( Matrix<F>& B,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo<Base<F>> LLL
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& R,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo<Base<F>> LLLWithQ
( Matrix<F>& B,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

template<typename F>
LLLInfo<Base<F>> LLLWithQ
( Matrix<F>& B,
  Matrix<F>& U,
  Matrix<F>& UInv,
  Matrix<F>& QR,
  Matrix<F>& t,
  Matrix<Base<F>>& d,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Perform a tree reduction of subsets of the original basis in order to 
// expose parallelism and perform as much work as possible in double-precision
// (which is often possible even for the SVP Challenge).
// This will not be substantially faster than the above LLL until Elemental
// supports different MPFR precisions simultaneously
template<typename F>
LLLInfo<Base<F>> RecursiveLLL
( Matrix<F>& B,
  Int cutoff=10,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );
template<typename F>
LLLInfo<Base<F>> RecursiveLLL
( Matrix<F>& B,
  Matrix<F>& R,
  Int cutoff=10,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Overwrite B, fill M with its (quasi-reduced) image of B, and fill K with the
// LLL-reduced basis for the kernel of B.
//
// This is essentially Algorithm 2.7.1 from Cohen's
// "A course in computational algebraic number theory". The main difference
// is that we avoid solving the normal equations and call a least squares
// solver.
// 
template<typename F>
void LatticeImageAndKernel
( Matrix<F>& B,
  Matrix<F>& M,
  Matrix<F>& K,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Overwrite B and fill K with the LLL-reduced basis for the kernel of B.
// This will eventually mirror Algorithm 2.7.2 from Cohen's
// "A course in computational algebraic number theory".
template<typename F>
void LatticeKernel
( Matrix<F>& B,
  Matrix<F>& K,
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Search for Z-dependence
// =======================
// Search for Z-dependence of a vector of real or complex numbers, z, via
// the quadratic form
//
//   Q(a) = || a ||_2^2 + N | z^T a |^2,
//
// which is generated by the basis matrix
//
//   
//   B = [I; sqrt(N) z^T],
//
// as Q(a) = a^T B^T B a = || B a ||_2^2. Cohen has advice for the choice of
// the (large) parameter N within subsection 2.7.2 within his book. 
//
// The number of (nearly) exact Z-dependences detected is returned.
//
template<typename F>
Int ZDependenceSearch
( const Matrix<F>& z,
        Base<F> NSqrt,
        Matrix<F>& B,
        Matrix<F>& U, 
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Search for an algebraic relation
// ================================
// Search for the (Gaussian) integer coefficients of a polynomial of alpha
// that is (nearly) zero.
template<typename F>
Int AlgebraicRelationSearch
( F alpha,
  Int n,
  Base<F> NSqrt,
  Matrix<F>& B,
  Matrix<F>& U, 
  const LLLCtrl<Base<F>>& ctrl=LLLCtrl<Base<F>>() );

// Schnorr-Euchner enumeration
// ===========================

namespace svp {

// If successful, fills 'v' with the integer coordinates of the columns of 
// an m x n matrix B (represented by its n x n upper-triangular Gaussian Normal
// Form; the 'R' from the QR factorization) which had a norm profile
// underneath the vector 'u' of upper bounds (|| (B v)(0:j) ||_2 < u(j)).
// Notice that the inequalities are strict.
//
// If not successful, the return value is a value greater than u(n-1) and 
// the contents of 'v' should be ignored.
//
// NOTE: There is not currently a complex implementation, though algorithms
//       exist.
template<typename F>
Base<F> BoundedEnumeration
( const Matrix<F>& R,
  const Matrix<Base<F>>& u,
        Matrix<F>& v );

} // namespace svp

// Given a reduced lattice B and its Gaussian Normal Form, R, either find a
// member of the lattice (given by B v, with v the output) with norm less than
// the upper bound (and return its norm), or return a value greater than the
// upper bound.
template<typename F>
Base<F> ShortVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  bool probabalistic=false );

// Given a reduced lattice B and its Gaussian Normal Form, R, find the shortest
// member of the lattice (with the shortest vector given by B v).
//
// The return value is the norm of the (approximately) shortest vector.
template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Matrix<F>& v,
  bool probabalistic=false );
// If an upper-bound on the shortest vector which is better than || b_0 ||_2 is
// available
template<typename F>
Base<F> ShortestVectorEnumeration
( const Matrix<F>& B,
  const Matrix<F>& R,
        Base<F> normUpperBound,
        Matrix<F>& v,
  bool probabalistic=false );

} // namespace El

#endif // ifndef EL_LATTICE_HPP
