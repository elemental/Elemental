/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MATRICES_HPP
#define EL_MATRICES_HPP

namespace El {

// Deterministic
// #############

// Bull's Head 
// ===========
template<typename Real>
void BullsHead( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void BullsHead( AbstractDistMatrix<Complex<Real>>& A, Int n );
template<typename Real>
void BullsHead( AbstractBlockDistMatrix<Complex<Real>>& A, Int n );

// Cauchy
// ======
template<typename F1,typename F2>
void Cauchy
( Matrix<F1>& A, const std::vector<F2>& x, const std::vector<F2>& y );
template<typename F1,typename F2>
void Cauchy
( AbstractDistMatrix<F1>& A,
  const std::vector<F2>& x, const std::vector<F2>& y );
template<typename F1,typename F2>
void Cauchy
( AbstractBlockDistMatrix<F1>& A,
  const std::vector<F2>& x, const std::vector<F2>& y );

// Cauchy-like
// ===========
template<typename F1,typename F2>
void CauchyLike
( Matrix<F1>& A, 
  const std::vector<F2>& r, const std::vector<F2>& s, 
  const std::vector<F2>& x, const std::vector<F2>& y );
template<typename F1,typename F2>
void CauchyLike
( AbstractDistMatrix<F1>& A,
  const std::vector<F2>& r, const std::vector<F2>& s, 
  const std::vector<F2>& x, const std::vector<F2>& y );
template<typename F1,typename F2>
void CauchyLike
( AbstractBlockDistMatrix<F1>& A,
  const std::vector<F2>& r, const std::vector<F2>& s, 
  const std::vector<F2>& x, const std::vector<F2>& y );

// Circulant
// =========
template<typename T>
void Circulant( Matrix<T>& A, const std::vector<T>& a );
template<typename T>
void Circulant( AbstractDistMatrix<T>& A, const std::vector<T>& a );
template<typename T>
void Circulant( AbstractBlockDistMatrix<T>& A, const std::vector<T>& a );

// Demmel
// ======
template<typename F>
void Demmel( Matrix<F>& A, Int n );
template<typename F>
void Demmel( AbstractDistMatrix<F>& A, Int n );
template<typename F>
void Demmel( AbstractBlockDistMatrix<F>& A, Int n );

// Diagonal
// ========
template<typename S,typename T>
void Diagonal( Matrix<S>& D, const std::vector<T>& d );
template<typename S,typename T>
void Diagonal( AbstractDistMatrix<S>& D, const std::vector<T>& d );
template<typename S,typename T>
void Diagonal( AbstractBlockDistMatrix<S>& D, const std::vector<T>& d );

// Ehrenfest
// =========
template<typename F>
void Ehrenfest( Matrix<F>& P, Int n );
template<typename F>
void Ehrenfest( AbstractDistMatrix<F>& P, Int n );
template<typename F>
void Ehrenfest( AbstractBlockDistMatrix<F>& P, Int n );

template<typename F>
void Ehrenfest( Matrix<F>& P, Matrix<F>& PInf, Int n );
template<typename F>
void Ehrenfest( AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& PInf, Int n );
template<typename F>
void Ehrenfest
( AbstractBlockDistMatrix<F>& P, AbstractBlockDistMatrix<F>& PInf, Int n );

template<typename F>
void EhrenfestStationary( Matrix<F>& PInf, Int n );
template<typename F>
void EhrenfestStationary( AbstractDistMatrix<F>& PInf, Int n );
template<typename F>
void EhrenfestStationary( AbstractBlockDistMatrix<F>& PInf, Int n );

template<typename F>
void EhrenfestDecay( Matrix<F>& A, Int n );
template<typename F,Dist U,Dist V>
void EhrenfestDecay( DistMatrix<F,U,V>& A, Int n );
/*
template<typename F,Dist U,Dist V>
void EhrenfestDecay( BlockDistMatrix<F,U,V>& A, Int n );
*/

// Extended Kahan
// ==============
template<typename F>
void ExtendedKahan( Matrix<F>& A, Int k, Base<F> phi, Base<F> mu );
template<typename F,Dist U,Dist V>
void ExtendedKahan
( DistMatrix<F,U,V>& A, Int k, Base<F> phi, Base<F> mu );

// Fiedler
// =======
template<typename F>
void Fiedler( Matrix<F>& A, const std::vector<F>& c );
template<typename F>
void Fiedler( AbstractDistMatrix<F>& A, const std::vector<F>& c );
template<typename F>
void Fiedler( AbstractBlockDistMatrix<F>& A, const std::vector<F>& c );

// Forsythe
// ========
template<typename T>
void Forsythe( Matrix<T>& J, Int n, T alpha, T lambda );
template<typename T>
void Forsythe( AbstractDistMatrix<T>& J, Int n, T alpha, T lambda );
template<typename T>
void Forsythe( AbstractBlockDistMatrix<T>& J, Int n, T alpha, T lambda );

// Fox-Li (aka the Landau matrix)
// ==============================
template<typename Real>
void FoxLi( Matrix<Complex<Real>>& A, Int n, Real omega );
template<typename Real,Dist U,Dist V>
void FoxLi( DistMatrix<Complex<Real>,U,V>& A, Int n, Real omega );

// Fourier
// =======
template<typename Real>
void Fourier( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void Fourier( AbstractDistMatrix<Complex<Real>>& A, Int n );
template<typename Real>
void Fourier( AbstractBlockDistMatrix<Complex<Real>>& A, Int n );

// Greatest Common Denominator matrix
// ==================================
template<typename T>
void GCDMatrix( Matrix<T>& G, Int m, Int n );
template<typename T>
void GCDMatrix( AbstractDistMatrix<T>& G, Int m, Int n );
template<typename T>
void GCDMatrix( AbstractBlockDistMatrix<T>& G, Int m, Int n );

// Gear matrix
// ===========
template<typename T>
void Gear( Matrix<T>& G, Int n, Int s, Int t );
template<typename T>
void Gear( AbstractDistMatrix<T>& G, Int n, Int s, Int t );
template<typename T>
void Gear( AbstractBlockDistMatrix<T>& G, Int n, Int s, Int t );

// Golub Klema Stewart matrix
// ==========================
template<typename F>
void GKS( Matrix<F>& A, Int n );
template<typename F>
void GKS( AbstractDistMatrix<F>& A, Int n );
template<typename F>
void GKS( AbstractBlockDistMatrix<F>& A, Int n );

// Grcar matrix
// ============
template<typename T>
void Grcar( Matrix<T>& A, Int n, Int k=3 );
template<typename T>
void Grcar( AbstractDistMatrix<T>& A, Int n, Int k=3 );
template<typename T>
void Grcar( AbstractBlockDistMatrix<T>& A, Int n, Int k=3 );

// Hankel
// ======
template<typename T>
void Hankel( Matrix<T>& A, Int m, Int n, const std::vector<T>& a );
template<typename T>
void Hankel( AbstractDistMatrix<T>& A, Int m, Int n, const std::vector<T>& a );
template<typename T>
void Hankel
( AbstractBlockDistMatrix<T>& A, Int m, Int n, const std::vector<T>& a );

// Hanowa
// ======
template<typename T>
void Hanowa( Matrix<T>& A, Int n, T mu );
template<typename T,Dist U,Dist V>
void Hanowa( DistMatrix<T,U,V>& A, Int n, T mu );
/*
template<typename T,Dist U,Dist V>
void Hanowa( BlockDistMatrix<T,U,V>& A, Int n, T mu );
*/

// Hatano-Nelson
// =============
template<typename F>
void HatanoNelson
( Matrix<F>& A, Int n, F center, Base<F> radius, F g, bool periodic=true );
template<typename F,Dist U,Dist V>
void HatanoNelson
( DistMatrix<F,U,V>& A, Int n, F center, Base<F> radius, F g,
  bool periodic=true );
template<typename F,Dist U,Dist V>
void HatanoNelson
( BlockDistMatrix<F,U,V>& A, Int n, F center, Base<F> radius, F g,
  bool periodic=true );

// Helmholtz
// =========
template<typename F>
void Helmholtz( Matrix<F>& H, Int nx, F shift );
template<typename F>
void Helmholtz( AbstractDistMatrix<F>& H, Int nx, F shift );
template<typename F>
void Helmholtz( Matrix<F>& H, Int nx, Int ny, F shift );
template<typename F>
void Helmholtz( AbstractDistMatrix<F>& H, Int nx, Int ny, F shift );
template<typename F>
void Helmholtz( Matrix<F>& H, Int nx, Int ny, Int nz, F shift );
template<typename F>
void Helmholtz( AbstractDistMatrix<F>& H, Int nx, Int ny, Int nz, F shift );

// Helmholtz PML
// =============
template<typename Real>
void HelmholtzPML
( Matrix<Complex<Real>>& H, Int nx, Complex<Real> omega,
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( DistMatrix<Complex<Real>>& H, Int nx, Complex<Real> omega,
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );

template<typename Real>
void HelmholtzPML
( Matrix<Complex<Real>>& H, Int nx, Int ny, Complex<Real> omega,
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( DistMatrix<Complex<Real>>& H, Int nx, Int ny, Complex<Real> omega,
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );

template<typename Real>
void HelmholtzPML
( Matrix<Complex<Real>>& H, Int nx, Int ny, Int nz, Complex<Real> omega,
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );
template<typename Real>
void HelmholtzPML
( DistMatrix<Complex<Real>>& H, Int nx, Int ny, Int nz, Complex<Real> omega,
  Int numPmlPoints=5, Real sigma=1.5, Real pmlExp=3 );

// Hermitian from EVD
// ==================
template<typename F>
void HermitianFromEVD
( UpperOrLower uplo, Matrix<F>& A,
  const Matrix<Base<F>>& w, const Matrix<F>& Z );
template<typename F>
void HermitianFromEVD
( UpperOrLower uplo, DistMatrix<F>& A,
  const DistMatrix<Base<F>,VR,STAR>& w, const DistMatrix<F>& Z );

// Hilbert
// =======
template<typename F>
void Hilbert( Matrix<F>& A, Int n );
template<typename F>
void Hilbert( AbstractDistMatrix<F>& A, Int n );
template<typename F>
void Hilbert( AbstractBlockDistMatrix<F>& A, Int n );

// Identity
// ========
template<typename T>
void MakeIdentity( Matrix<T>& I );
template<typename T> 
void MakeIdentity( AbstractDistMatrix<T>& I );
template<typename T> 
void MakeIdentity( AbstractBlockDistMatrix<T>& I );

template<typename T> 
void Identity( Matrix<T>& I, Int m, Int n );
template<typename T> 
void Identity( AbstractDistMatrix<T>& I, Int m, Int n );
template<typename T> 
void Identity( AbstractBlockDistMatrix<T>& I, Int m, Int n );

// Jordan
// ======
template<typename T>
void Jordan( Matrix<T>& J, Int n, T lambda );
template<typename T>
void Jordan( AbstractDistMatrix<T>& J, Int n, T lambda );
template<typename T>
void Jordan( AbstractBlockDistMatrix<T>& J, Int n, T lambda );

// Kahan
// =====
template<typename F>
void Kahan( Matrix<F>& A, Int n, F phi );
template<typename F>
void Kahan( AbstractDistMatrix<F>& A, Int n, F phi );
template<typename F>
void Kahan( AbstractBlockDistMatrix<F>& A, Int n, F phi );

// KMS
// ===
template<typename T>
void KMS( Matrix<T>& K, Int n, T rho );
template<typename T>
void KMS( AbstractDistMatrix<T>& K, Int n, T rho );
template<typename T>
void KMS( AbstractBlockDistMatrix<T>& K, Int n, T rho );

// Lehmer
// ======
template<typename F>
void Lehmer( Matrix<F>& L, Int n );
template<typename F>
void Lehmer( AbstractDistMatrix<F>& L, Int n );
template<typename F>
void Lehmer( AbstractBlockDistMatrix<F>& L, Int n );

// Laplacian
// =========
template<typename F>
void Laplacian( Matrix<F>& L, Int nx );
template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int nx );
template<typename F>
void Laplacian( AbstractBlockDistMatrix<F>& L, Int nx );

template<typename F>
void Laplacian( Matrix<F>& L, Int nx, Int ny );
template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny );
template<typename F>
void Laplacian( AbstractBlockDistMatrix<F>& L, Int nx, Int ny );

template<typename F>
void Laplacian( Matrix<F>& L, Int nx, Int ny, Int nz );
template<typename F>
void Laplacian( AbstractDistMatrix<F>& L, Int nx, Int ny, Int nz );
template<typename F>
void Laplacian( AbstractBlockDistMatrix<F>& L, Int nx, Int ny, Int nz );

// Lauchli
// =======
template<typename T>
void Lauchli( Matrix<T>& A, Int n, T mu );
template<typename T,Dist U,Dist V>
void Lauchli( DistMatrix<T,U,V>& A, Int n, T mu );
/*
template<typename T,Dist U,Dist V>
void Lauchli( BlockDistMatrix<T,U,V>& A, Int n, T mu );
*/

// Legendre
// ========
template<typename F>
void Legendre( Matrix<F>& A, Int n );
template<typename F>
void Legendre( AbstractDistMatrix<F>& A, Int n );
template<typename F>
void Legendre( AbstractBlockDistMatrix<F>& A, Int n );

// Lotkin
// ======
template<typename F>
void Lotkin( Matrix<F>& A, Int n );
template<typename F>
void Lotkin( AbstractDistMatrix<F>& A, Int n );
template<typename F>
void Lotkin( AbstractBlockDistMatrix<F>& A, Int n );

// MinIJ
// =====
template<typename T>
void MinIJ( Matrix<T>& M, Int n );
template<typename T>
void MinIJ( AbstractDistMatrix<T>& M, Int n );
template<typename T>
void MinIJ( AbstractBlockDistMatrix<T>& M, Int n );

// Normal from EVD
// ===============
template<typename Real>
void NormalFromEVD
(       Matrix<Complex<Real>>& A,
  const Matrix<Complex<Real>>& w,
  const Matrix<Complex<Real>>& Z );
template<typename Real>
void NormalFromEVD
(       DistMatrix<Complex<Real>>& A,
  const DistMatrix<Complex<Real>,VR,STAR>& w,
  const DistMatrix<Complex<Real>>& Z );

// Ones
// ====
template<typename T>
void Ones( Matrix<T>& A, Int m, Int n );
template<typename T>
void Ones( AbstractDistMatrix<T>& A, Int m, Int n );
template<typename T>
void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n );

// 1-2-1 matrix
// ============
template<typename T>
void OneTwoOne( Matrix<T>& A, Int n );
template<typename T>
void OneTwoOne( AbstractDistMatrix<T>& A, Int n );
template<typename T>
void OneTwoOne( AbstractBlockDistMatrix<T>& A, Int n );

// Parter
// ======
template<typename F>
void Parter( Matrix<F>& P, Int n );
template<typename F>
void Parter( AbstractDistMatrix<F>& P, Int n );
template<typename F>
void Parter( AbstractBlockDistMatrix<F>& P, Int n );

// Pei
// ===
template<typename T>
void Pei( Matrix<T>& P, Int n, T alpha );
template<typename T>
void Pei( AbstractDistMatrix<T>& P, Int n, T alpha );
template<typename T>
void Pei( AbstractBlockDistMatrix<T>& P, Int n, T alpha );

// Redheffer
// =========
template<typename T>
void Redheffer( Matrix<T>& R, Int n );
template<typename T>
void Redheffer( AbstractDistMatrix<T>& R, Int n );
template<typename T>
void Redheffer( AbstractBlockDistMatrix<T>& R, Int n );

// Riemann
// =======
template<typename T>
void Riemann( Matrix<T>& R, Int n );
template<typename T>
void Riemann( AbstractDistMatrix<T>& R, Int n );
template<typename T>
void Riemann( AbstractBlockDistMatrix<T>& R, Int n );

// Riffle
// ======
template<typename F>
void Riffle( Matrix<F>& P, Int n );
template<typename F>
void Riffle( AbstractDistMatrix<F>& P, Int n );
template<typename F>
void Riffle( AbstractBlockDistMatrix<F>& P, Int n );

template<typename F>
void Riffle
( Matrix<F>& P, Matrix<F>& PInf, Int n );
template<typename F>
void Riffle
( AbstractDistMatrix<F>& P, AbstractDistMatrix<F>& PInf, Int n );
template<typename F>
void Riffle
( AbstractBlockDistMatrix<F>& P, AbstractBlockDistMatrix<F>& PInf, Int n );

template<typename F>
void RiffleStationary( Matrix<F>& PInf, Int n );
template<typename F>
void RiffleStationary( AbstractDistMatrix<F>& PInf, Int n );
template<typename F>
void RiffleStationary( AbstractBlockDistMatrix<F>& PInf, Int n );

template<typename F>
void RiffleDecay( Matrix<F>& A, Int n );
template<typename F,Dist U,Dist V>
void RiffleDecay( DistMatrix<F,U,V>& A, Int n );
/*
template<typename F,Dist U,Dist V>
void RiffleDecay( BlockDistMatrix<F,U,V>& A, Int n );
*/

// Ris
// ===
template<typename F>
void Ris( Matrix<F>& R, Int n );
template<typename F>
void Ris( AbstractDistMatrix<F>& R, Int n );
template<typename F>
void Ris( AbstractBlockDistMatrix<F>& R, Int n );

// Toeplitz
// ========
template<typename S,typename T>
void Toeplitz
( Matrix<S>& A, Int m, Int n, const std::vector<T>& a );
template<typename S,typename T>
void Toeplitz
( AbstractDistMatrix<S>& A, Int m, Int n, const std::vector<T>& a );
template<typename S,typename T>
void Toeplitz
( AbstractBlockDistMatrix<S>& A, Int m, Int n, const std::vector<T>& a );

// Trefethen
// =========
template<typename Real>
void Trefethen( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void Trefethen( AbstractDistMatrix<Complex<Real>>& A, Int n );
template<typename Real>
void Trefethen( AbstractBlockDistMatrix<Complex<Real>>& A, Int n );

// Triangle
// ========
template<typename F>
void Triangle( Matrix<F>& A, Int n );
template<typename F>
void Triangle( AbstractDistMatrix<F>& A, Int n );
template<typename F>
void Triangle( AbstractBlockDistMatrix<F>& A, Int n );

// TriW
// ====
template<typename T>
void TriW( Matrix<T>& A, Int m, Int n, T alpha, Int k );
template<typename T>
void TriW( AbstractDistMatrix<T>& A, Int m, Int n, T alpha, Int k );
template<typename T>
void TriW( AbstractBlockDistMatrix<T>& A, Int m, Int n, T alpha, Int k );

// Walsh
// =====
template<typename T>
void Walsh( Matrix<T>& A, Int k, bool binary=false );
template<typename T>
void Walsh( AbstractDistMatrix<T>& A, Int k, bool binary=false );

// Whale
// =====
template<typename Real>
void Whale( Matrix<Complex<Real>>& A, Int n );
template<typename Real>
void Whale( AbstractDistMatrix<Complex<Real>>& A, Int n );
template<typename Real>
void Whale( AbstractBlockDistMatrix<Complex<Real>>& A, Int n );

// Wilkinson
// =========
template<typename T>
void Wilkinson( Matrix<T>& A, Int k );
template<typename T>
void Wilkinson( AbstractDistMatrix<T>& A, Int k );

// Zeros
// =====
template<typename T>
void Zeros( Matrix<T>& A, Int m, Int n );
template<typename T>
void Zeros( AbstractDistMatrix<T>& A, Int m, Int n );
template<typename T>
void Zeros( AbstractBlockDistMatrix<T>& A, Int m, Int n );

// Random
// ######

// Gaussian
// ========
template<typename T>
void MakeGaussian( Matrix<T>& A, T mean=0, Base<T> stddev=1 );
template<typename T>
void MakeGaussian( AbstractDistMatrix<T>& A, T mean=0, Base<T> stddev=1 );
template<typename T>
void MakeGaussian( AbstractBlockDistMatrix<T>& A, T mean=0, Base<T> stddev=1 );

template<typename T>
void Gaussian( Matrix<T>& A, Int m, Int n, T mean=0, Base<T> stddev=1 );
template<typename T>
void Gaussian
( AbstractDistMatrix<T>& A, Int m, Int n, T mean=0, Base<T> stddev=1 );
template<typename T>
void Gaussian
( AbstractBlockDistMatrix<T>& A, Int m, Int n, T mean=0, Base<T> stddev=1 );

// Haar
// ====
template<typename F> 
void Haar( Matrix<F>& A, Int n );
template<typename F> 
void Haar( DistMatrix<F>& A, Int n );

template<typename F> 
void ImplicitHaar( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d, Int n );
template<typename F> 
void ImplicitHaar
( DistMatrix<F>& A,
  DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d, Int n );

// Hermitian uniform spectrum
// ==========================
template<typename F>
void HermitianUniformSpectrum
( Matrix<F>& A, Int n, Base<F> lower=0, Base<F> upper=1 );
template<typename F>
void HermitianUniformSpectrum
( DistMatrix<F>& A, Int n, Base<F> lower=0, Base<F> upper=1 );
template<typename F,Dist U,Dist V>
void HermitianUniformSpectrum
( DistMatrix<F,U,V>& A, Int n, Base<F> lower=0, Base<F> upper=1 );

// Normal uniform spectrum
// =======================
template<typename Real>
void NormalUniformSpectrum
( Matrix<Complex<Real>>& A, Int n,
  Complex<Real> center=0, Real radius=1 );
template<typename Real>
void NormalUniformSpectrum
( DistMatrix<Complex<Real>>& A, Int n,
  Complex<Real> center=0, Real radius=1 );
template<typename Real,Dist U,Dist V>
void NormalUniformSpectrum
( DistMatrix<Complex<Real>,U,V>& A, Int n,
  Complex<Real> center=0, Real radius=1 );

// Uniform
// =======
// Draw each entry from a uniform PDF over a closed ball.
template<typename T>
void MakeUniform( Matrix<T>& A, T center=0, Base<T> radius=1 );
template<typename T>
void MakeUniform( AbstractDistMatrix<T>& A, T center=0, Base<T> radius=1 );
template<typename T>
void MakeUniform( AbstractBlockDistMatrix<T>& A, T center=0, Base<T> radius=1 );

template<typename T>
void Uniform( Matrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 );
template<typename T>
void Uniform
( AbstractDistMatrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 );
template<typename T>
void Uniform
( AbstractBlockDistMatrix<T>& A, Int m, Int n, T center=0, Base<T> radius=1 );

// Uniform Helmholtz Green's
// =========================
template<typename Real>
void UniformHelmholtzGreens
( Matrix<Complex<Real>>& A, Int n, Real lambda );
template<typename Real>
void UniformHelmholtzGreens
( AbstractDistMatrix<Complex<Real>>& A, Int n, Real lambda );
template<typename Real>
void UniformHelmholtzGreens
( AbstractBlockDistMatrix<Complex<Real>>& A, Int n, Real lambda );

// Wigner
// ======
template<typename T>
void Wigner( Matrix<T>& A, Int n, T mean=0, Base<T> stddev=1 );
template<typename T>
void Wigner( AbstractDistMatrix<T>& A, Int n, T mean=0, Base<T> stddev=1 );

} // namespace El

#include "./matrices/Egorov.hpp"

#endif // ifndef EL_MATRICES_HPP
