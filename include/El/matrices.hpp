/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_MATRICES_HPP
#define EL_MATRICES_HPP

namespace El {

// Deterministic
// #############

// Classical
// =========

// Cauchy
// ------
template<typename F1,typename F2>
void Cauchy
( Matrix<F1>& A, const vector<F2>& x, const vector<F2>& y );
template<typename F1,typename F2>
void Cauchy
( AbstractDistMatrix<F1>& A,
  const vector<F2>& x, const vector<F2>& y );

// Cauchy-like
// -----------
template<typename F1,typename F2>
void CauchyLike
( Matrix<F1>& A,
  const vector<F2>& r, const vector<F2>& s,
  const vector<F2>& x, const vector<F2>& y );
template<typename F1,typename F2>
void CauchyLike
( AbstractDistMatrix<F1>& A,
  const vector<F2>& r, const vector<F2>& s,
  const vector<F2>& x, const vector<F2>& y );

// Circulant
// ---------
template<typename T>
void Circulant( Matrix<T>& A, const Matrix<T>& a );
template<typename T>
void Circulant( Matrix<T>& A, const vector<T>& a );

template<typename T>
void Circulant( AbstractDistMatrix<T>& A, const Matrix<T>& a );
template<typename T>
void Circulant( AbstractDistMatrix<T>& A, const vector<T>& a );

// Diagonal
// --------
template<typename S,typename T>
void Diagonal( Matrix<S>& D, const vector<T>& d );
template<typename S,typename T>
void Diagonal( Matrix<S>& D, const Matrix<T>& d );
template<typename S,typename T>
void Diagonal( AbstractDistMatrix<S>& D, const vector<T>& d );
template<typename S,typename T>
void Diagonal( AbstractDistMatrix<S>& D, const Matrix<T>& d );
template<typename S,typename T>
void Diagonal( AbstractDistMatrix<S>& D, const AbstractDistMatrix<T>& d );

template<typename S,typename T>
void Diagonal( SparseMatrix<S>& D, const Matrix<T>& d );
template<typename S,typename T>
void Diagonal( DistSparseMatrix<S>& D, const DistMultiVec<T>& d );

// Egorov
// ------
template<typename Real>
void Egorov
( Matrix<Complex<Real>>& A, function<Real(Int,Int)> phase, Int n );
template<typename Real>
void Egorov
( AbstractDistMatrix<Complex<Real>>& A,
  function<Real(Int,Int)> phase, Int n );

// Fiedler
// -------
template<typename Field>
void Fiedler( Matrix<Field>& A, const vector<Field>& c );
template<typename Field>
void Fiedler( AbstractDistMatrix<Field>& A, const vector<Field>& c );

// Fourier
// -------
template<typename Real>
void Fourier( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void Fourier( AbstractDistMatrix<Complex<Real>>& A, Int n );

// Greatest Common Denominator matrix
// ----------------------------------
template<typename T>
void GCDMatrix( Matrix<T>& G, Int m, Int n );
template<typename T>
void GCDMatrix( AbstractDistMatrix<T>& G, Int m, Int n );

// Hankel
// ------
template<typename T>
void Hankel( Matrix<T>& A, Int m, Int n, const vector<T>& a );
template<typename T>
void Hankel( AbstractDistMatrix<T>& A, Int m, Int n, const vector<T>& a );

// Hilbert
// -------
template<typename Field>
void Hilbert( Matrix<Field>& A, Int n );
template<typename Field>
void Hilbert( AbstractDistMatrix<Field>& A, Int n );

// Identity
// --------
template<typename T>
void MakeIdentity( Matrix<T>& I );
template<typename T>
void MakeIdentity( AbstractDistMatrix<T>& I );

template<typename T>
void Identity( Matrix<T>& I, Int m, Int n );
template<typename T>
void Identity( AbstractDistMatrix<T>& I, Int m, Int n );
template<typename T>
void Identity( SparseMatrix<T>& I, Int m, Int n );
template<typename T>
void Identity( DistSparseMatrix<T>& I, Int m, Int n );

// Jordan
// ------
template<typename T>
void Jordan( Matrix<T>& J, Int n, T lambda );
template<typename T>
void Jordan( AbstractDistMatrix<T>& J, Int n, T lambda );

// Ones
// ----
template<typename T>
void Ones( Matrix<T>& A, Int m, Int n );
template<typename T>
void Ones( AbstractDistMatrix<T>& A, Int m, Int n );
template<typename T>
void Ones( DistMultiVec<T>& A, Int m, Int n );
template<typename T>
void Ones( SparseMatrix<T>& A, Int m, Int n );
template<typename T>
void Ones( DistSparseMatrix<T>& A, Int m, Int n );

// Toeplitz
// --------
template<typename S,typename T>
void Toeplitz
( Matrix<S>& A, Int m, Int n, const vector<T>& a );
template<typename S,typename T>
void Toeplitz
( AbstractDistMatrix<S>& A, Int m, Int n, const vector<T>& a );

// Walsh
// -----
template<typename T>
void Walsh( Matrix<T>& A, Int k, bool binary=false );
template<typename T>
void Walsh( AbstractDistMatrix<T>& A, Int k, bool binary=false );

// Zeros
// -----
template<typename T>
void Zeros( Matrix<T>& A, Int m, Int n );
template<typename T>
void Zeros( AbstractDistMatrix<T>& A, Int m, Int n );
template<typename T>
void Zeros( SparseMatrix<T>& A, Int m, Int n );
template<typename T>
void Zeros( DistSparseMatrix<T>& A, Int m, Int n );
template<typename T>
void Zeros( DistMultiVec<T>& A, Int m, Int n );

// Integral equations
// ==================

// Fox-Li (aka the Landau matrix)
// ------------------------------
template<typename Real>
void FoxLi( Matrix<Complex<Real>>& A, Int n, Real omega );
template<typename Real>
void FoxLi( AbstractDistMatrix<Complex<Real>>& A, Int n, Real omega );

// Partial differential equations
// ==============================

// Helmholtz
// ---------
template<typename Field>
void Helmholtz( Matrix<Field>& H, Int nx, Field shift );
template<typename Field>
void Helmholtz( AbstractDistMatrix<Field>& H, Int nx, Field shift );
template<typename Field>
void Helmholtz( SparseMatrix<Field>& H, Int nx, Field shift );
template<typename Field>
void Helmholtz( DistSparseMatrix<Field>& H, Int nx, Field shift );

template<typename Field>
void Helmholtz( Matrix<Field>& H, Int nx, Int ny, Field shift );
template<typename Field>
void Helmholtz( AbstractDistMatrix<Field>& H, Int nx, Int ny, Field shift );
template<typename Field>
void Helmholtz( SparseMatrix<Field>& H, Int nx, Int ny, Field shift );
template<typename Field>
void Helmholtz( DistSparseMatrix<Field>& H, Int nx, Int ny, Field shift );

template<typename Field>
void Helmholtz
( Matrix<Field>& H, Int nx, Int ny, Int nz, Field shift );
template<typename Field>
void Helmholtz
( AbstractDistMatrix<Field>& H, Int nx, Int ny, Int nz, Field shift );
template<typename Field>
void Helmholtz
( SparseMatrix<Field>& H, Int nx, Int ny, Int nz, Field shift );
template<typename Field>
void Helmholtz
( DistSparseMatrix<Field>& H, Int nx, Int ny, Int nz, Field shift );

// Helmholtz PML
// -------------
template<typename Real>
void HelmholtzPML
( Matrix<Complex<Real>>& H, Int nx,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( AbstractDistMatrix<Complex<Real>>& H, Int nx,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( SparseMatrix<Complex<Real>>& H, Int nx,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( DistSparseMatrix<Complex<Real>>& H, Int nx,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );

template<typename Real>
void HelmholtzPML
( Matrix<Complex<Real>>& H, Int nx, Int ny,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( AbstractDistMatrix<Complex<Real>>& H, Int nx, Int ny,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( SparseMatrix<Complex<Real>>& H, Int nx, Int ny,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( DistSparseMatrix<Complex<Real>>& H, Int nx, Int ny,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );

template<typename Real>
void HelmholtzPML
( Matrix<Complex<Real>>& H, Int nx, Int ny, Int nz,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( AbstractDistMatrix<Complex<Real>>& H, Int nx, Int ny, Int nz,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( SparseMatrix<Complex<Real>>& H, Int nx, Int ny, Int nz,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( DistSparseMatrix<Complex<Real>>& H, Int nx, Int ny, Int nz,
  Complex<Real> omega, Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );

// Laplacian
// ---------
template<typename Field>
void Laplacian( Matrix<Field>& L, Int nx );
template<typename Field>
void Laplacian( AbstractDistMatrix<Field>& L, Int nx );
template<typename Field>
void Laplacian( SparseMatrix<Field>& L, Int nx );
template<typename Field>
void Laplacian( DistSparseMatrix<Field>& L, Int nx );

template<typename Field>
void Laplacian( Matrix<Field>& L, Int nx, Int ny );
template<typename Field>
void Laplacian( AbstractDistMatrix<Field>& L, Int nx, Int ny );
template<typename Field>
void Laplacian( SparseMatrix<Field>& L, Int nx, Int ny );
template<typename Field>
void Laplacian( DistSparseMatrix<Field>& L, Int nx, Int ny );

template<typename Field>
void Laplacian( Matrix<Field>& L, Int nx, Int ny, Int nz );
template<typename Field>
void Laplacian( AbstractDistMatrix<Field>& L, Int nx, Int ny, Int nz );
template<typename Field>
void Laplacian( SparseMatrix<Field>& L, Int nx, Int ny, Int nz );
template<typename Field>
void Laplacian( DistSparseMatrix<Field>& L, Int nx, Int ny, Int nz );

// Miscellaneous (to be categorized)
// =================================

// Demmel
// ------
template<typename Field>
void Demmel( Matrix<Field>& A, Int n );
template<typename Field>
void Demmel( AbstractDistMatrix<Field>& A, Int n );

// Druinsky-Toledo matrices
// ------------------------
// An example of Bunch-Kaufman A producing large element growth in
// floating-point arithmetic. Please see Theorem 5 from:
//     http://www.alexdruinsky.com/pdfs/bkbound-revised.pdf
template<typename Field>
void DruinskyToledo( Matrix<Field>& A, Int n );
template<typename Field>
void DruinskyToledo( ElementalMatrix<Field>& A, Int n );

// DynamicRegCounter
// -----------------
// An example of a well-conditioned sparse matrix that, for n >= 100, causes
// catastrophic failure of dynamic regularization since the first n pivots are
// all one, but the upper-left quadrant has a condition number behaving as
// O(2^n).
//
// The matrix is of the form
//
//   A = | 4 J_{1/2}(n)^T J_{1/2}(n)   I |,
//       |            I               -I |
//
// where J_{lambda}(n) is the n x n Jordan block with eigenvalue lambda.
//
// The upper-left corner corresponds to the matrix formed by JordanCholesky,
// whose Cholesky factor is 2 J_{1/2}(n), which has a unit diagonal but is
// incredibly ill-conditioned.
//
template<typename T>
void DynamicRegCounter( Matrix<T>& A, Int n );
template<typename T>
void DynamicRegCounter( ElementalMatrix<T>& A, Int n );
template<typename T>
void DynamicRegCounter( SparseMatrix<T>& A, Int n );
template<typename T>
void DynamicRegCounter( DistSparseMatrix<T>& A, Int n );

// Ehrenfest
// ---------
template<typename Field>
void Ehrenfest( Matrix<Field>& P, Int n );
template<typename Field>
void Ehrenfest( AbstractDistMatrix<Field>& P, Int n );

template<typename Field>
void Ehrenfest( Matrix<Field>& P, Matrix<Field>& PInf, Int n );
template<typename Field>
void Ehrenfest
( ElementalMatrix<Field>& P, ElementalMatrix<Field>& PInf, Int n );

template<typename Field>
void EhrenfestStationary( Matrix<Field>& PInf, Int n );
template<typename Field>
void EhrenfestStationary( AbstractDistMatrix<Field>& PInf, Int n );

template<typename Field>
void EhrenfestDecay( Matrix<Field>& A, Int n );
template<typename Field>
void EhrenfestDecay( ElementalMatrix<Field>& A, Int n );

// Extended Kahan
// --------------
template<typename Field>
void ExtendedKahan
( Matrix<Field>& A, Int k, Base<Field> phi, Base<Field> mu );
template<typename Field>
void ExtendedKahan
( ElementalMatrix<Field>& A, Int k, Base<Field> phi, Base<Field> mu );

// Gear matrix
// -----------
template<typename T>
void Gear( Matrix<T>& G, Int n, Int s, Int t );
template<typename T>
void Gear( AbstractDistMatrix<T>& G, Int n, Int s, Int t );

// Gaussian Elimination with Partial Pivoting Growth
// -------------------------------------------------
template<typename Field>
void GEPPGrowth( Matrix<Field>& A, Int n );
template<typename Field>
void GEPPGrowth( ElementalMatrix<Field>& A, Int n );

// Golub Klema Stewart matrix
// --------------------------
template<typename Field>
void GKS( Matrix<Field>& A, Int n );
template<typename Field>
void GKS( AbstractDistMatrix<Field>& A, Int n );

// Hanowa
// ------
template<typename T>
void Hanowa( Matrix<T>& A, Int n, T mu );
template<typename T>
void Hanowa( ElementalMatrix<T>& A, Int n, T mu );

// JordanCholesky
// --------------
// Set A = 4 J_{1/2}(n)^T J_{1/2}(n), which has is tridiagonal with its main
// diagonal equal to five everywhere but the top-left entry, and its sub and
// super-diagonals equal to 2.
template<typename T>
void JordanCholesky( Matrix<T>& A, Int n );
template<typename T>
void JordanCholesky( AbstractDistMatrix<T>& A, Int n );
template<typename T>
void JordanCholesky( SparseMatrix<T>& A, Int n );
template<typename T>
void JordanCholesky( DistSparseMatrix<T>& A, Int n );

// Kahan
// -----
template<typename Field>
void Kahan( Matrix<Field>& A, Int n, Field phi );
template<typename Field>
void Kahan( AbstractDistMatrix<Field>& A, Int n, Field phi );

// KMS
// ---
template<typename T>
void KMS( Matrix<T>& K, Int n, T rho );
template<typename T>
void KMS( AbstractDistMatrix<T>& K, Int n, T rho );

// Lauchli
// -------
template<typename T>
void Lauchli( Matrix<T>& A, Int n, T mu );
template<typename T>
void Lauchli( ElementalMatrix<T>& A, Int n, T mu );

// Legendre
// --------
template<typename Field>
void Legendre( Matrix<Field>& A, Int n );
template<typename Field>
void Legendre( AbstractDistMatrix<Field>& A, Int n );

// Lehmer
// ------
template<typename Field>
void Lehmer( Matrix<Field>& L, Int n );
template<typename Field>
void Lehmer( AbstractDistMatrix<Field>& L, Int n );

// Lotkin
// ------
template<typename Field>
void Lotkin( Matrix<Field>& A, Int n );
template<typename Field>
void Lotkin( AbstractDistMatrix<Field>& A, Int n );

// MinIJ
// -----
template<typename T>
void MinIJ( Matrix<T>& M, Int n );
template<typename T>
void MinIJ( AbstractDistMatrix<T>& M, Int n );

// Parter
// ------
template<typename Field>
void Parter( Matrix<Field>& P, Int n );
template<typename Field>
void Parter( AbstractDistMatrix<Field>& P, Int n );

// Pei
// ---
template<typename T>
void Pei( Matrix<T>& P, Int n, T alpha );
template<typename T>
void Pei( AbstractDistMatrix<T>& P, Int n, T alpha );

// Redheffer
// ---------
template<typename T>
void Redheffer( Matrix<T>& R, Int n );
template<typename T>
void Redheffer( AbstractDistMatrix<T>& R, Int n );

// Riffle
// ------
template<typename Field>
void Riffle( Matrix<Field>& P, Int n );
template<typename Field>
void Riffle( AbstractDistMatrix<Field>& P, Int n );

template<typename Field>
void Riffle
( Matrix<Field>& P, Matrix<Field>& PInf, Int n );
template<typename Field>
void Riffle
( ElementalMatrix<Field>& P, ElementalMatrix<Field>& PInf, Int n );

template<typename Field>
void RiffleStationary( Matrix<Field>& PInf, Int n );
template<typename Field>
void RiffleStationary( AbstractDistMatrix<Field>& PInf, Int n );

template<typename Field>
void RiffleDecay( Matrix<Field>& A, Int n );
template<typename Field>
void RiffleDecay( ElementalMatrix<Field>& A, Int n );

// Ris
// ---
template<typename Field>
void Ris( Matrix<Field>& R, Int n );
template<typename Field>
void Ris( AbstractDistMatrix<Field>& R, Int n );

// Wilkinson
// ---------
template<typename T>
void Wilkinson( Matrix<T>& A, Int k );
template<typename T>
void Wilkinson( AbstractDistMatrix<T>& A, Int k );

// Sparse Toeplitz
// ===============

// Bull's Head
// -----------
template<typename Real>
void BullsHead( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void BullsHead( AbstractDistMatrix<Complex<Real>>& A, Int n );

// Forsythe
// --------
template<typename T>
void Forsythe( Matrix<T>& J, Int n, T alpha, T lambda );
template<typename T>
void Forsythe( AbstractDistMatrix<T>& J, Int n, T alpha, T lambda );

// Grcar matrix
// ------------
template<typename T>
void Grcar( Matrix<T>& A, Int n, Int k=3 );
template<typename T>
void Grcar( AbstractDistMatrix<T>& A, Int n, Int k=3 );

// 1-2-1 matrix
// ------------
template<typename T>
void OneTwoOne( Matrix<T>& A, Int n );
template<typename T>
void OneTwoOne( AbstractDistMatrix<T>& A, Int n );

// Trefethen-Embree
// ----------------
template<typename Real>
void TrefethenEmbree( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void TrefethenEmbree( AbstractDistMatrix<Complex<Real>>& A, Int n );

// Triangle
// --------
template<typename Field>
void Triangle( Matrix<Field>& A, Int n );
template<typename Field>
void Triangle( AbstractDistMatrix<Field>& A, Int n );

// TriW
// ----
template<typename T>
void TriW( Matrix<T>& A, Int n, T alpha, Int k );
template<typename T>
void TriW( AbstractDistMatrix<T>& A, Int n, T alpha, Int k );

// Whale
// -----
template<typename Real>
void Whale( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void Whale( AbstractDistMatrix<Complex<Real>>& A, Int n );

// Lattice
// =======

// The transpose of the 2 x 2 block matrix described in sub-subsection 3.4.1
// in:
//
//   Jeffrey Hoffstein, Jill Pipher, Joseph H. Silverman,
//   "NTRU: A ring-based public key cryptosystem"
//
// NOTE: While 'q' should be an integer, we accept it as floating-point for now
template<typename Real>
void NTRUAttack( Matrix<Real>& A, const Matrix<Real>& h, Real alpha, Real q );

// Random
// ######

// Independent
// ===========

// Bernoulli
// ---------
template<typename T>
void Bernoulli( Matrix<T>& A, Int m, Int n, double p=0.5 );
template<typename T>
void Bernoulli( AbstractDistMatrix<T>& A, Int m, Int n, double p=0.5 );

// Gaussian
// --------
template<typename Field>
void MakeGaussian
( Matrix<Field>& A, Field mean=0, Base<Field> stddev=1 );
template<typename Field>
void MakeGaussian
( AbstractDistMatrix<Field>& A, Field mean=0, Base<Field> stddev=1 );
template<typename Field>
void MakeGaussian
( DistMultiVec<Field>& A, Field mean=0, Base<Field> stddev=1 );

template<typename Field>
void Gaussian
( Matrix<Field>& A, Int m, Int n,
  Field mean=0, Base<Field> stddev=1 );
template<typename Field>
void Gaussian
( AbstractDistMatrix<Field>& A, Int m, Int n,
  Field mean=0, Base<Field> stddev=1 );
template<typename Field>
void Gaussian
( DistMultiVec<Field>& A, Int m, Int n,
  Field mean=0, Base<Field> stddev=1 );

// Rademacher
// ----------
template<typename T>
void Rademacher( Matrix<T>& A, Int m, Int n );
template<typename T>
void Rademacher( AbstractDistMatrix<T>& A, Int m, Int n );

// Three-valued
// ------------
template<typename T>
void ThreeValued( Matrix<T>& A, Int m, Int n, double p=2./3. );
template<typename T>
void ThreeValued( AbstractDistMatrix<T>& A, Int m, Int n, double p=2./3. );

// Uniform
// -------
// Draw each entry from a uniform PDF over a closed ball.
template<typename T>
void MakeUniform( Matrix<T>& A, T center=0, Base<T> radius=1 );
template<typename T>
void MakeUniform( AbstractDistMatrix<T>& A, T center=0, Base<T> radius=1 );
template<typename T>
void MakeUniform( DistMultiVec<T>& X, T center=0, Base<T> radius=1 );

template<typename T>
void Uniform( Matrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 );
template<typename T>
void Uniform
( AbstractDistMatrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 );
template<typename T>
void Uniform( DistMultiVec<T>& X, Int m, Int n, T center=0, Base<T> radius=1 );

// Lattice bases
// =============

// AjtaiTypeBasis
// --------------
// Generalization of of the Ajtai-type bases of dimension d described in
// Subsection 3.4 of Nguyen and Stehle's "LLL on the Average". Note that we
// use the transpose of their scheme and generalize to complex arithmetic by
// preserving real diagonal entries but sampling the off-diagonal from the
// unit ball of radius B_{i,i}/2.
//
// NOTE: Be careful of overflow for modest values of n and alpha
template<typename Field>
void AjtaiTypeBasis( Matrix<Field>& A, Int n, Base<Field> alpha );
template<typename Field>
void AjtaiTypeBasis( AbstractDistMatrix<Field>& A, Int n, Base<Field> alpha );

// KnapsackTypeBasis
// -----------------
// See Subsection 3.4 of Nguyen and Stehle's "LLL on the Average".
// Note that we use the transpose of the convention of said paper.
// The result is an (n+1) x n matrix whose bottom row is sampled from the
// integers lying within the unit ball of radius 'radius'
template<typename Field>
void KnapsackTypeBasis
( Matrix<Field>& A, Int n, Base<Field> radius );
template<typename Field>
void KnapsackTypeBasis
( AbstractDistMatrix<Field>& A, Int n, Base<Field> radius );

// Miscellaneous (to be categorized)
// =================================

// Hatano-Nelson
// -------------
template<typename Field>
void HatanoNelson
( Matrix<Field>& A, Int n, Field center, Base<Field> radius, Field g,
  bool periodic=true );
template<typename Field>
void HatanoNelson
( ElementalMatrix<Field>& A, Int n, Field center, Base<Field> radius, Field g,
  bool periodic=true );

// Haar
// ----
template<typename Field>
void Haar( Matrix<Field>& A, Int n );
template<typename Field>
void Haar( ElementalMatrix<Field>& A, Int n );

template<typename Field>
void ImplicitHaar
( Matrix<Field>& A,
  Matrix<Field>& householderScalars,
  Matrix<Base<Field>>& signature, Int n );
template<typename Field>
void ImplicitHaar
( ElementalMatrix<Field>& A,
  ElementalMatrix<Field>& householderScalars,
  ElementalMatrix<Base<Field>>& signature, Int n );

// Hermitian uniform spectrum
// --------------------------
template<typename Field>
void HermitianUniformSpectrum
( Matrix<Field>& A, Int n, Base<Field> lower=0, Base<Field> upper=1 );
template<typename Field>
void HermitianUniformSpectrum
( ElementalMatrix<Field>& A, Int n, Base<Field> lower=0, Base<Field> upper=1 );

// Normal uniform spectrum
// -----------------------
template<typename Real>
void NormalUniformSpectrum
( Matrix<Complex<Real>>& A, Int n,
  Complex<Real> center=0, Real radius=1 );
template<typename Real>
void NormalUniformSpectrum
( ElementalMatrix<Complex<Real>>& A, Int n,
  Complex<Real> center=0, Real radius=1 );

// Uniform Helmholtz Green's
// -------------------------
template<typename Real>
void UniformHelmholtzGreens
( Matrix<Complex<Real>>& A, Int n, Real lambda );
template<typename Real>
void UniformHelmholtzGreens
( ElementalMatrix<Complex<Real>>& A, Int n, Real lambda );

// Wigner
// ------
template<typename Field>
void Wigner
( Matrix<Field>& A, Int n, Field mean=0, Base<Field> stddev=1 );
template<typename Field>
void Wigner
( ElementalMatrix<Field>& A, Int n, Field mean=0, Base<Field> stddev=1 );

} // namespace El

// TODO(poulson): Group these into a small number of includes of parent dir's
#include <El/matrices/deterministic/classical/Circulant.hpp>
#include <El/matrices/deterministic/lattice/NTRUAttack.hpp>

#endif // ifndef EL_MATRICES_HPP
