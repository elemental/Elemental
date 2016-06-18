/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/blas_like/level1.hpp>
#include <El/lapack_like/funcs.hpp>
#include <El/lapack_like/euclidean_min.hpp>
#include <El/control.hpp>

namespace El {

// W = | A^H  L |, where K and L are Hermitian.
//     | K   -A |
//
// The solution, X, to the equation
//   X K X - A^H X - X A = L
// is returned, as well as the number of Newton iterations for computing sgn(W).
//
// See Chapter 2 of Nicholas J. Higham's "Functions of Matrices"

template<typename F>
void Riccati
( Matrix<F>& W, Matrix<F>& X, SignCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
    Sign( W, ctrl );
    const Int n = W.Height()/2;
    Matrix<F> WTL, WTR,
              WBL, WBR;
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    // (ML, MR) = sgn(W) - I
    ShiftDiagonal( W, F(-1) );

    // Solve for X in ML X = -MR
    Matrix<F> ML, MR;
    PartitionRight( W, ML, MR, n );
    MR *= -1;
    ls::Overwrite( NORMAL, ML, MR, X );
}

template<typename F>
void Riccati
( ElementalMatrix<F>& WPre,
  ElementalMatrix<F>& X, 
  SignCtrl<Base<F>> ctrl )
{
    DEBUG_CSE

    DistMatrixReadProxy<F,F,MC,MR> WProx( WPre );
    auto& W = WProx.Get();

    const Grid& g = W.Grid();
    Sign( W, ctrl );
    const Int n = W.Height()/2;
    DistMatrix<F> WTL(g), WTR(g),
                  WBL(g), WBR(g);
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    // (ML, MR) = sgn(W) - I
    ShiftDiagonal( W, F(-1) );

    // Solve for X in ML X = -MR
    DistMatrix<F> ML(g), MR(g);
    PartitionRight( W, ML, MR, n );
    MR *= -1;
    ls::Overwrite( NORMAL, ML, MR, X );
}

template<typename F>
void Riccati
( UpperOrLower uplo, 
  const Matrix<F>& A,
  const Matrix<F>& K,
  const Matrix<F>& L,
        Matrix<F>& X,
  SignCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( K.Height() != K.Width() )
          LogicError("K must be square");
      if( L.Height() != L.Width() )
          LogicError("L must be square");
      if( A.Height() != K.Height() || A.Height() != L.Height() )
          LogicError("A, K, and L must be the same size");
    )
    const Int n = A.Height();
    Matrix<F> W, WTL, WTR,
                 WBL, WBR;
    W.Resize( 2*n, 2*n );
    Zero( W );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    Adjoint( A, WTL );
    WBR = A; WBR *= -1;
    WBL = K; MakeHermitian( uplo, WBL );
    WTR = L; MakeHermitian( uplo, WTR );

    Riccati( W, X, ctrl );
}

template<typename F>
void Riccati
( UpperOrLower uplo, 
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& K, 
  const ElementalMatrix<F>& L,
        ElementalMatrix<F>& X, 
  SignCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( K.Height() != K.Width() )
          LogicError("K must be square");
      if( L.Height() != L.Width() )
          LogicError("L must be square");
      if( A.Height() != K.Height() || A.Height() != L.Height() )
          LogicError("A, K, and L must be the same size");
      AssertSameGrids( A, K, L );
    )
    const Grid& g = A.Grid();
    const Int n = A.Height();
    DistMatrix<F> W(g), WTL(g), WTR(g),
                        WBL(g), WBR(g);
    W.Resize( 2*n, 2*n );
    Zero( W );
    PartitionDownDiagonal
    ( W, WTL, WTR,
         WBL, WBR, n );

    Adjoint( A, WTL );
    WBR = A; WBR *= -1;
    WBL = K; MakeHermitian( uplo, WBL );
    WTR = L; MakeHermitian( uplo, WTR );

    Riccati( W, X, ctrl );
}

#define PROTO(F) \
  template void Riccati \
  ( Matrix<F>& W, \
    Matrix<F>& X, \
    SignCtrl<Base<F>> ctrl ); \
  template void Riccati \
  ( ElementalMatrix<F>& W, \
    ElementalMatrix<F>& X, \
    SignCtrl<Base<F>> ctrl ); \
  template void Riccati \
  ( UpperOrLower uplo, \
    const Matrix<F>& A, \
    const Matrix<F>& K, \
    const Matrix<F>& L, \
          Matrix<F>& X, \
    SignCtrl<Base<F>> ctrl ); \
  template void Riccati \
  ( UpperOrLower uplo, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& K, \
    const ElementalMatrix<F>& L, \
          ElementalMatrix<F>& X, \
    SignCtrl<Base<F>> ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_QUAD
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
