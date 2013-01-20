/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_DECL_HPP
#define MATRICES_DECL_HPP

namespace elem {

//----------------------------------------------------------------------------//
// Deterministic                                                              //
//----------------------------------------------------------------------------//

// Generate an m x n Cauchy matrix, i.e., 
//
//   A(i,j) = 1/(x(i)-y(j)), where x(i) != y(j), 
//                           for all (i,j) such that 0 <= i < m, 0 <= j < n.
//
template<typename F>
void Cauchy
( const std::vector<F>& x, const std::vector<F>& y, Matrix<F>& A );
template<typename F,Distribution U,Distribution V>
void Cauchy
( const std::vector<F>& x, const std::vector<F>& y, DistMatrix<F,U,V>& A );

// Generate an m x n Cauchy-like matrix, i.e., 
//
//   A(i,j) = r(i)s(j)/(x(i)-y(j)), where x(i) != y(j), 
//                                  for all (i,j) such that 
//                                  0 <= i < m, 0 <= j < n.
//
template<typename F>
void CauchyLike
( const std::vector<F>& r, const std::vector<F>& s,
  const std::vector<F>& x, const std::vector<F>& y,
  Matrix<F>& A );
template<typename F,Distribution U,Distribution V>
void CauchyLike
( const std::vector<F>& r, const std::vector<F>& s,
  const std::vector<F>& x, const std::vector<F>& y,
  DistMatrix<F,U,V>& A );

// Generate an n x n circulant matrix, i.e., 
//
//   A(i,j) = a((i-j) mod n), for all (i,j) such that 0 <= i,j < n.
//
// Note that circulant matrices are special cases of Toeplitz matrices, but 
// that our indexing scheme for 'a' is different.
//
template<typename T>
void Circulant( const std::vector<T>& a, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Circulant( const std::vector<T>& a, DistMatrix<T,U,V>& A );

// Generate an n x n diagonal matrix, i.e., 
//
//   D(i,j) = d(j), if i = j, and zero otherwise.
//
template<typename T>
void Diagonal( const std::vector<T>& d, Matrix<T>& D );
template<typename T,Distribution U,Distribution V>
void Diagonal( const std::vector<T>& d, DistMatrix<T,U,V>& D );

// Generate an n x n Discrete Fourier Transform matrix, i.e., 
//
//   A(k,l) = exp(-2*pi*i*k*l/n) / sqrt(n).
//
template<typename Z>
void DiscreteFourier( int n, Matrix<Complex<Z> >& A );
template<typename Z,Distribution U,Distribution V>
void DiscreteFourier( int n, DistMatrix<Complex<Z>,U,V>& A );
// Turn the existing square matrix into a DFT matrix
template<typename Z>
void MakeDiscreteFourier( int n, Matrix<Complex<Z> >& A );
template<typename Z,Distribution U,Distribution V>
void MakeDiscreteFourier( int n, DistMatrix<Complex<Z>,U,V>& A );

// Generate an m x n Hankel matrix, i.e.,
//
//    A(i,j) = a(i+j), for all (i,j) such that 0 <= i < m, 0 <= j < n.
//
template<typename T>
void Hankel( int m, int n, const std::vector<T>& a, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Hankel( int m, int n, const std::vector<T>& a, DistMatrix<T,U,V>& A );

// Generate an n x n Hilbert matrix, i.e., 
//
//   H(i,j) = 1/(i+j+1), for all (i,j) such that 0 <= i, j < n.
//
template<typename F>
void Hilbert( int n, Matrix<F>& H );
template<typename F,Distribution U,Distribution V>
void Hilbert( int n, DistMatrix<F,U,V>& H );
// Turn the existing square matrix into a Hilbert matrix
template<typename F>
void MakeHilbert( int n, Matrix<F>& H );
template<typename F,Distribution U,Distribution V>
void MakeHilbert( int n, DistMatrix<F,U,V>& H );

// Generate an m x n identity-like matrix, i.e., 
//
//   I(i,j) = 1, if i = j, and zero otherwise.
//
template<typename T>
void Identity( int m, int n, Matrix<T>& I );
template<typename T,Distribution U,Distribution V>
void Identity( int m, int n, DistMatrix<T,U,V>& I );
// Turn the existing matrix into an identity-like matrix
template<typename T>
void MakeIdentity( Matrix<T>& I );
template<typename T,Distribution U,Distribution V>
void MakeIdentity( DistMatrix<T,U,V>& I );

// Generate an m x n matrix of all ones.
template<typename T>
void Ones( int m, int n, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Ones( int m, int n, DistMatrix<T,U,V>& A ); 
// Turn the existing matrix into a matrix of all ones.
template<typename T>
void MakeOnes( Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void MakeOnes( DistMatrix<T,U,V>& A );

// Generate a so-called "1-2-1" matrix of order k (and dimension k), 
// which is tridiagonal with a diagonal of all twos and sub- and super-diagonals
// of all ones.
template<typename T>
void OneTwoOne( int n, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void OneTwoOne( int n, DistMatrix<T,U,V>& A );
// Turn the existing matrix into a 1-2-1 matrix.
template<typename T>
void MakeOneTwoOne( int n, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void MakeOneTwoOne( int n, DistMatrix<T,U,V>& A );

// Generate an m x n Toeplitz matrix, i.e.,
//
//   A(i,j) = a(i-j+(n-1)), for all (i,j) such that 0 <= i < m, 0 <= j < n.
//
template<typename T>
void Toeplitz( int m, int n, const std::vector<T>& a, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Toeplitz( int m, int n, const std::vector<T>& a, DistMatrix<T,U,V>& A );

// Generate a Walsh matrix of order k (and dimension 2^k), where
//
//    W_1 = | 1  1 |, and
//          | 1 -1 |
//
//    W_k = | W(k-1)  W(k-1) |, for k >= 2.
//          | W(k-1) -W(k-1) |
//
// A binary Walsh matrix is the same as above, but with all -1 entries replaced
// with zeros.
//
template<typename T>
void Walsh( int k, Matrix<T>& W, bool binary=false );
template<typename T,Distribution U,Distribution V>
void Walsh( int k, DistMatrix<T,U,V>& W, bool binary=false );

// Generate a Wilkinson matrix of order k (and dimension 2k+1), 
// which is tridiagonal with
//             
//   diag(W_k) = [k,k-1,k-2,...,1,0,1,...,k-2,k-1,k],
// 
// and sub- and super-diagonals of all ones.
template<typename T>
void Wilkinson( int k, Matrix<T>& W );
template<typename T,Distribution U,Distribution V>
void Wilkinson( int k, DistMatrix<T,U,V>& W );

// Generate an m x n matrix of all zeros
template<typename T>
void Zeros( int m, int n, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Zeros( int m, int n, DistMatrix<T,U,V>& A );
// Turn the existing matrix into a matrix of all zeros
template<typename T>
void MakeZeros( Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void MakeZeros( DistMatrix<T,U,V>& A );

//----------------------------------------------------------------------------//
// Random                                                                     //
//----------------------------------------------------------------------------//

// Generate an m x n matrix of samples from the uniform PDF over the 
// closed unit ball.
template<typename T>
void Uniform
( int m, int n, Matrix<T>& A, 
  T center=0, typename Base<T>::type radius=1 );
template<typename T,Distribution U,Distribution V>
void Uniform
( int m, int n, DistMatrix<T,U,V>& A, 
  T center=0, typename Base<T>::type radius=1 );
// Turn the existing matrix into a uniform random matrix
template<typename T>
void MakeUniform
( Matrix<T>& A, T center=0, typename Base<T>::type radius=1 );
template<typename T,Distribution U,Distribution V>
void MakeUniform
( DistMatrix<T,U,V>& A, T center=0, typename Base<T>::type radius=1 );

// Choose the spectrum from a uniform distribution over the specified 
// half-open interval, (lower,upper], and then perform a random Householder
// similarity transformation
template<typename F>
void HermitianUniformSpectrum
( int n, Matrix<F>& A, 
  typename Base<F>::type lower=0, typename Base<F>::type upper=1 );
template<typename F,Distribution U,Distribution V>
void HermitianUniformSpectrum
( int n, DistMatrix<F,U,V>& A, 
  typename Base<F>::type lower=0, typename Base<F>::type upper=1 );
// Turn the existing matrix into a Hermitian uniform random matrix
template<typename F>
void MakeHermitianUniformSpectrum
( Matrix<F>& A, 
  typename Base<F>::type lower=0, typename Base<F>::type upper=1 );
template<typename F,Distribution U,Distribution V>
void MakeHermitianUniformSpectrum
( DistMatrix<F,U,V>& A, 
  typename Base<F>::type lower=0, typename Base<F>::type upper=1 );

// Choose the spectrum from a uniform distribution over the specified 
// ball, B_radius(center), and then perform a random Householder
// similarity transformation
template<typename R>
void NormalUniformSpectrum
( int n, Matrix<Complex<R> >& A, 
  Complex<R> center=0, R radius=1 );
template<typename R,Distribution U,Distribution V>
void NormalUniformSpectrum
( int n, DistMatrix<Complex<R>,U,V>& A, Complex<R> center=0, R radius=1 );
// Turn the existing matrix into a normal uniform random matrix
template<typename R>
void MakeNormalUniformSpectrum
( Matrix<Complex<R> >& A, Complex<R> center=0, R radius=1 );
template<typename R,Distribution U,Distribution V>
void MakeNormalUniformSpectrum
( DistMatrix<Complex<R>,U,V>& A, Complex<R> center=0, R radius=1 );

// TODO: Gaussian random matrices

} // namespace elem

#endif // ifndef MATRICES_DECL_HPP
