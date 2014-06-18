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

// Ones
// ====
template<typename T>
void Ones( Matrix<T>& A, Int m, Int n );
template<typename T>
void Ones( AbstractDistMatrix<T>& A, Int m, Int n );
template<typename T>
void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n );

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

// Wigner
// ======
template<typename T>
void Wigner( Matrix<T>& A, Int n, T mean=0, Base<T> stddev=1 );
template<typename T>
void Wigner( AbstractDistMatrix<T>& A, Int n, T mean=0, Base<T> stddev=1 );

} // namespace El

#include "./matrices/Egorov.hpp"
#include "./matrices/GCDMatrix.hpp"
#include "./matrices/Gear.hpp"
#include "./matrices/GKS.hpp"
#include "./matrices/Grcar.hpp"
#include "./matrices/Hankel.hpp"
#include "./matrices/Hanowa.hpp"
#include "./matrices/HatanoNelson.hpp"
#include "./matrices/Helmholtz.hpp"
#include "./matrices/HelmholtzPML.hpp"
#include "./matrices/HermitianFromEVD.hpp"
#include "./matrices/Hilbert.hpp"
#include "./matrices/Kahan.hpp"
#include "./matrices/KMS.hpp"
#include "./matrices/Laplacian.hpp"
#include "./matrices/Lauchli.hpp"
#include "./matrices/Legendre.hpp"
#include "./matrices/Lehmer.hpp"
#include "./matrices/Lotkin.hpp"
#include "./matrices/MinIJ.hpp"
#include "./matrices/NormalFromEVD.hpp"
#include "./matrices/OneTwoOne.hpp"
#include "./matrices/Parter.hpp"
#include "./matrices/Pei.hpp"
#include "./matrices/Redheffer.hpp"
#include "./matrices/Riemann.hpp"
#include "./matrices/Riffle.hpp"
#include "./matrices/Ris.hpp"
#include "./matrices/Toeplitz.hpp"
#include "./matrices/Trefethen.hpp"
#include "./matrices/Triangle.hpp"
#include "./matrices/TriW.hpp"
#include "./matrices/Walsh.hpp"
#include "./matrices/Whale.hpp"
#include "./matrices/Wilkinson.hpp"

// Random matrices
// ===============

// Uniform
#include "./matrices/HermitianUniformSpectrum.hpp"
#include "./matrices/NormalUniformSpectrum.hpp"
#include "./matrices/UniformHelmholtzGreens.hpp"

#endif // ifndef EL_MATRICES_HPP
