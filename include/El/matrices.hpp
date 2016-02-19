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
template<typename F>
void Fiedler( Matrix<F>& A, const vector<F>& c );
template<typename F>
void Fiedler( AbstractDistMatrix<F>& A, const vector<F>& c );

// Fourier
// -------
template<typename Real>
void Fourier( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void Fourier( AbstractDistMatrix<Complex<Real>>& A, Int n );

// Fourier-Identity
// ----------------
// A common example of a low-coherence n x 2n matrix, [F I]
template<typename Real>
void FourierIdentity( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void FourierIdentity( ElementalMatrix<Complex<Real>>& A, Int n );

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
template<typename F>
void Hilbert( Matrix<F>& A, Int n );
template<typename F>
void Hilbert( AbstractDistMatrix<F>& A, Int n );

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

// Walsh-Identity
// --------------
template<typename T>
void WalshIdentity( Matrix<T>& A, Int k, bool binary=false );
template<typename T>
void WalshIdentity( ElementalMatrix<T>& A, Int k, bool binary=false );

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
void FoxLi( ElementalMatrix<Complex<Real>>& A, Int n, Real omega );

// Partial differential equations
// ==============================

// Helmholtz
// ---------
template<typename F>
void Helmholtz( Matrix<F>& H, Int nx, F shift );
template<typename F>
void Helmholtz( AbstractDistMatrix<F>& H, Int nx, F shift );
template<typename F>
void Helmholtz( SparseMatrix<F>& H, Int nx, F shift );
template<typename F>
void Helmholtz( DistSparseMatrix<F>& H, Int nx, F shift );

template<typename F>
void Helmholtz( Matrix<F>& H, Int nx, Int ny, F shift );
template<typename F>
void Helmholtz( AbstractDistMatrix<F>& H, Int nx, Int ny, F shift );
template<typename F>
void Helmholtz( SparseMatrix<F>& H, Int nx, Int ny, F shift );
template<typename F>
void Helmholtz( DistSparseMatrix<F>& H, Int nx, Int ny, F shift );

template<typename F>
void Helmholtz( Matrix<F>& H, Int nx, Int ny, Int nz, F shift );
template<typename F>
void Helmholtz( AbstractDistMatrix<F>& H, Int nx, Int ny, Int nz, F shift );
template<typename F>
void Helmholtz( SparseMatrix<F>& H, Int nx, Int ny, Int nz, F shift );
template<typename F>
void Helmholtz( DistSparseMatrix<F>& H, Int nx, Int ny, Int nz, F shift );

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
template<typename F>
void Laplacian( Matrix<F>& L, Int nx );
template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int nx );
template<typename F>
void Laplacian( SparseMatrix<F>& L, Int nx );
template<typename F>
void Laplacian( DistSparseMatrix<F>& L, Int nx );

template<typename F>
void Laplacian( Matrix<F>& L, Int nx, Int ny );
template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny );
template<typename F>
void Laplacian( SparseMatrix<F>& L, Int nx, Int ny );
template<typename F>
void Laplacian( DistSparseMatrix<F>& L, Int nx, Int ny );

template<typename F>
void Laplacian( Matrix<F>& L, Int nx, Int ny, Int nz );
template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny, Int nz );
template<typename F>
void Laplacian( SparseMatrix<F>& L, Int nx, Int ny, Int nz );
template<typename F>
void Laplacian( DistSparseMatrix<F>& L, Int nx, Int ny, Int nz );

// Miscellaneous (to be categorized)
// =================================

// Demmel
// ------
template<typename F>
void Demmel( Matrix<F>& A, Int n );
template<typename F>
void Demmel( AbstractDistMatrix<F>& A, Int n );

// Druinsky-Toledo matrices
// ------------------------
// An example of Bunch-Kaufman A producing large element growth in 
// floating-point arithmetic. Please see Theorem 5 from:
//     http://www.alexdruinsky.com/pdfs/bkbound-revised.pdf
template<typename F>
void DruinskyToledo( Matrix<F>& A, Int n );
template<typename F>
void DruinskyToledo( ElementalMatrix<F>& A, Int n );

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
template<typename F>
void Ehrenfest( Matrix<F>& P, Int n );
template<typename F>
void Ehrenfest( AbstractDistMatrix<F>& P, Int n );

template<typename F>
void Ehrenfest( Matrix<F>& P, Matrix<F>& PInf, Int n );
template<typename F>
void Ehrenfest( ElementalMatrix<F>& P, ElementalMatrix<F>& PInf, Int n );

template<typename F>
void EhrenfestStationary( Matrix<F>& PInf, Int n );
template<typename F>
void EhrenfestStationary( AbstractDistMatrix<F>& PInf, Int n );

template<typename F>
void EhrenfestDecay( Matrix<F>& A, Int n );
template<typename F>
void EhrenfestDecay( ElementalMatrix<F>& A, Int n );

// Extended Kahan
// --------------
template<typename F>
void ExtendedKahan( Matrix<F>& A, Int k, Base<F> phi, Base<F> mu );
template<typename F>
void ExtendedKahan( ElementalMatrix<F>& A, Int k, Base<F> phi, Base<F> mu );

// Gear matrix
// -----------
template<typename T>
void Gear( Matrix<T>& G, Int n, Int s, Int t );
template<typename T>
void Gear( AbstractDistMatrix<T>& G, Int n, Int s, Int t );

// Gaussian Elimination with Partial Pivoting Growth
// -------------------------------------------------
template<typename F>
void GEPPGrowth( Matrix<F>& A, Int n );
template<typename F>
void GEPPGrowth( ElementalMatrix<F>& A, Int n );

// Golub Klema Stewart matrix
// --------------------------
template<typename F>
void GKS( Matrix<F>& A, Int n );
template<typename F>
void GKS( AbstractDistMatrix<F>& A, Int n );

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
template<typename F>
void Kahan( Matrix<F>& A, Int n, F phi );
template<typename F>
void Kahan( AbstractDistMatrix<F>& A, Int n, F phi );

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
template<typename F>
void Legendre( Matrix<F>& A, Int n );
template<typename F>
void Legendre( AbstractDistMatrix<F>& A, Int n );

// Lehmer
// ------
template<typename F>
void Lehmer( Matrix<F>& L, Int n );
template<typename F>
void Lehmer( AbstractDistMatrix<F>& L, Int n );

// Lotkin
// ------
template<typename F>
void Lotkin( Matrix<F>& A, Int n );
template<typename F>
void Lotkin( AbstractDistMatrix<F>& A, Int n );

// MinIJ
// -----
template<typename T>
void MinIJ( Matrix<T>& M, Int n );
template<typename T>
void MinIJ( AbstractDistMatrix<T>& M, Int n );

// Parter
// ------
template<typename F>
void Parter( Matrix<F>& P, Int n );
template<typename F>
void Parter( AbstractDistMatrix<F>& P, Int n );

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
template<typename F>
void Riffle( Matrix<F>& P, Int n );
template<typename F>
void Riffle( AbstractDistMatrix<F>& P, Int n );

template<typename F>
void Riffle
( Matrix<F>& P, Matrix<F>& PInf, Int n );
template<typename F>
void Riffle
( ElementalMatrix<F>& P, ElementalMatrix<F>& PInf, Int n );

template<typename F>
void RiffleStationary( Matrix<F>& PInf, Int n );
template<typename F>
void RiffleStationary( AbstractDistMatrix<F>& PInf, Int n );

template<typename F>
void RiffleDecay( Matrix<F>& A, Int n );
template<typename F>
void RiffleDecay( ElementalMatrix<F>& A, Int n );

// Ris
// ---
template<typename F>
void Ris( Matrix<F>& R, Int n );
template<typename F>
void Ris( AbstractDistMatrix<F>& R, Int n );

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
template<typename F>
void Triangle( Matrix<F>& A, Int n );
template<typename F>
void Triangle( AbstractDistMatrix<F>& A, Int n );

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
template<typename F>
void MakeGaussian( Matrix<F>& A, F mean=0, Base<F> stddev=1 );
template<typename F>
void MakeGaussian( AbstractDistMatrix<F>& A, F mean=0, Base<F> stddev=1 );
template<typename F>
void MakeGaussian( DistMultiVec<F>& A, F mean=0, Base<F> stddev=1 );

template<typename F>
void Gaussian( Matrix<F>& A, Int m, Int n, F mean=0, Base<F> stddev=1 );
template<typename F>
void Gaussian
( AbstractDistMatrix<F>& A, Int m, Int n, F mean=0, Base<F> stddev=1 );
template<typename F>
void Gaussian
( DistMultiVec<F>& A, Int m, Int n, F mean=0, Base<F> stddev=1 );

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
template<typename F>
void AjtaiTypeBasis( Matrix<F>& A, Int n, Base<F> alpha );
template<typename F>
void AjtaiTypeBasis( AbstractDistMatrix<F>& A, Int n, Base<F> alpha );

// KnapsackTypeBasis
// -----------------
// See Subsection 3.4 of Nguyen and Stehle's "LLL on the Average".
// Note that we use the transpose of the convention of said paper.
// The result is an (n+1) x n matrix whose bottom row is sampled from the 
// integers lying within the unit ball of radius 'radius'
template<typename F>
void KnapsackTypeBasis( Matrix<F>& A, Int n, Base<F> radius );
template<typename F>
void KnapsackTypeBasis( AbstractDistMatrix<F>& A, Int n, Base<F> radius );

// Miscellaneous (to be categorized)
// =================================

// Hatano-Nelson
// -------------
template<typename F>
void HatanoNelson
( Matrix<F>& A, Int n, F center, Base<F> radius, F g,
  bool periodic=true );
template<typename F>
void HatanoNelson
( ElementalMatrix<F>& A, Int n, F center, Base<F> radius, F g,
  bool periodic=true );

// Haar
// ----
template<typename F> 
void Haar( Matrix<F>& A, Int n );
template<typename F> 
void Haar( ElementalMatrix<F>& A, Int n );

template<typename F> 
void ImplicitHaar
( Matrix<F>& A,
  Matrix<F>& t,
  Matrix<Base<F>>& d, Int n );
template<typename F> 
void ImplicitHaar
( ElementalMatrix<F>& A,
  ElementalMatrix<F>& t,
  ElementalMatrix<Base<F>>& d, Int n );

// Hermitian uniform spectrum
// --------------------------
template<typename F>
void HermitianUniformSpectrum
( Matrix<F>& A, Int n, Base<F> lower=0, Base<F> upper=1 );
template<typename F>
void HermitianUniformSpectrum
( ElementalMatrix<F>& A, Int n, Base<F> lower=0, Base<F> upper=1 );

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
template<typename F>
void Wigner( Matrix<F>& A, Int n, F mean=0, Base<F> stddev=1 );
template<typename F>
void Wigner( ElementalMatrix<F>& A, Int n, F mean=0, Base<F> stddev=1 );

} // namespace El

// TODO: Group these into a small number of includes of parent dir's
#include <El/matrices/deterministic/classical/Circulant.hpp>
#include <El/matrices/deterministic/lattice/NTRUAttack.hpp>

#endif // ifndef EL_MATRICES_HPP
