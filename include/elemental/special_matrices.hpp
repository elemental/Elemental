/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef ELEMENTAL_SPECIAL_MATRICES_HPP
#define ELEMENTAL_SPECIAL_MATRICES_HPP 1

namespace elem {

// TODO: Add support for random matrices

// Generate an m x n Cauchy matrix, i.e., 
//
//   A(i,j) = 1/(x(i)-y(j)), where x(i) != y(j), 
//                           for all (i,j) such that 0 <= i < m, 0 <= j < n.
//
template<typename F>
void Cauchy
( int m, int n, const std::vector<F>& x, const std::vector<F>& y,
  Matrix<F>& A );
template<typename F,Distribution U,Distribution V>
void Cauchy
( int m, int n, const std::vector<F>& x, const std::vector<F>& y,
  DistMatrix<F,U,V>& A );

// Generate an m x n Cauchy-like matrix, i.e., 
//
//   A(i,j) = r(i)s(j)/(x(i)-y(j)), where x(i) != y(j), 
//                                  for all (i,j) such that 
//                                  0 <= i < m, 0 <= j < n.
//
template<typename F>
void CauchyLike
( int m, int n, 
  const std::vector<F>& r, const std::vector<F>& s,
  const std::vector<F>& x, const std::vector<F>& y,
  Matrix<F>& A );
template<typename F,Distribution U,Distribution V>
void CauchyLike
( int m, int n, 
  const std::vector<F>& r, const std::vector<F>& s,
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
void Circulant( int n, const std::vector<T>& a, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Circulant( int n, const std::vector<T>& a, DistMatrix<T,U,V>& A );

// Generate an n x n diagonal matrix, i.e., 
//
//   D(i,j) = d(j), if i = j, and zero otherwise.
//
template<typename T>
void Diagonal( int n, const std::vector<T>& d, Matrix<T>& D );
template<typename T,Distribution U,Distribution V>
void Diagonal( int n, const std::vector<T>& d, DistMatrix<T,U,V>& D );

// Generate an m x n Hankel matrix, i.e.,
//
//    A(i,j) = a(i+j+1), for all (i,j) such that 0 <= i < m, 0 <= j < n.
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

// Generate an m x n identity-like matrix, i.e., 
//
//   I(i,j) = 1, if i = j, and zero otherwise.
//
template<typename T>
void Identity( int m, int n, Matrix<T>& I );
template<typename T,Distribution U,Distribution V>
void Identity( int m, int n, DistMatrix<T,U,V>& I );

// Force the matrix to be lower or upper trapezoidal, with the diagonal
// defined relative to either the top-left or bottom-right corner of the matrix
// (based on the 'side' parameter). The 'offset' parameter determines where the
// last nonzero diagonal is, with '0' meaning the main diagonal, '1' meaning 
// the superdiagonal, '-1' meaning the subdiagonal, and likewise for all other
// integer values.
template<typename T>
void MakeTrapezoidal
( LeftOrRight side, UpperOrLower uplo, int offset, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void MakeTrapezoidal
( LeftOrRight side, UpperOrLower uplo, int offset, DistMatrix<T,U,V>& A );

// Generate an m x n matrix of all ones.
template<typename T>
void Ones( int m, int n, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Ones( int m, int n, DistMatrix<T,U,V>& A ); 

// Generate a so-called "1-2-1" matrix of order k (and dimension k), 
// which is tridiagonal with a diagonal of all twos and sub- and super-diagonals
// of all ones.
template<typename T>
void OneTwoOne( int k, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void OneTwoOne( int k, DistMatrix<T,U,V>& A );

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
template<typename T>
void Walsh( int k, Matrix<T>& W );
template<typename T,Distribution U,Distribution V>
void Walsh( int k, DistMatrix<T,U,V>& W );

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

// Zero the contents of a matrix
template<typename T>
void Zero( Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Zero( DistMatrix<T,U,V>& A );

// Generate an m x n matrix of all zeros
template<typename T>
void Zeros( int m, int n, Matrix<T>& A );
template<typename T,Distribution U,Distribution V>
void Zeros( int m, int n, DistMatrix<T,U,V>& A );

} // namespace elem

#endif /* ELEMENTAL_SPECIAL_HPP */

