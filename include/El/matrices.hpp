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

} // namespace El

#include "./matrices/Diagonal.hpp"
#include "./matrices/Egorov.hpp"
#include "./matrices/Ehrenfest.hpp"
#include "./matrices/ExtendedKahan.hpp"
#include "./matrices/Fiedler.hpp"
#include "./matrices/Forsythe.hpp"
#include "./matrices/FoxLi.hpp"
#include "./matrices/Fourier.hpp"
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
#include "./matrices/Identity.hpp"
#include "./matrices/Jordan.hpp"
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

// Gaussian
#include "./matrices/Gaussian.hpp"
#include "./matrices/Wigner.hpp"
#include "./matrices/Haar.hpp"

#endif // ifndef EL_MATRICES_HPP
