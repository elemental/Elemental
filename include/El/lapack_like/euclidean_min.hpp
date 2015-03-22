/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_EUCLIDEANMIN_HPP
#define EL_EUCLIDEANMIN_HPP

#include "El/lapack_like/factor.hpp"

namespace El {

// Generalized (Gauss-Markov) Linear Model
// =======================================
// Solve 
//   min_{X,Y} || Y ||_F subject to D = A X + B Y
template<typename F>
void GLM
( Matrix<F>& A, Matrix<F>& B, Matrix<F>& D, Matrix<F>& Y );
template<typename F>
void GLM
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, 
  AbstractDistMatrix<F>& D, AbstractDistMatrix<F>& Y );
// TODO: Sparse-direct implementations

// Least Squares (or Minimum Length)
// =================================
// When height(op(A)) >= width(op(A)), solve
//
//    min_X || op(A) X - B ||_F,
//
// otherwise, solve 
//
//    min_X || X ||_F s.t. op(A) X = B.
//
template<typename F>
void LeastSquares
( Orientation orientation, 
  Matrix<F>& A, const Matrix<F>& B, 
                      Matrix<F>& X );
template<typename F>
void LeastSquares
( Orientation orientation, 
  AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, 
                                  AbstractDistMatrix<F>& X );

template<typename Real>
struct LeastSquaresCtrl {
    Real alpha;
    RegQSDCtrl<Real> qsdCtrl;
    bool equilibrate;
    bool progress;
    bool time;

    LeastSquaresCtrl()
    : alpha(Pow(lapack::MachineEpsilon<Real>(),Real(0.25))),
      equilibrate(true), progress(false), time(false)
    { }
};

template<typename F>
void LeastSquares
( Orientation orientation,
  const SparseMatrix<F>& A, const Matrix<F>& Y, 
                                  Matrix<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );
template<typename F>
void LeastSquares
( Orientation orientation,
  const DistSparseMatrix<F>& A, const DistMultiVec<F>& Y, 
                                      DistMultiVec<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );

// Equality-constrained Least Squarees
// ===================================
// Solve
//   min_X || A X - C ||_F subject to B X = D
template<typename F>
void LSE
( Matrix<F>& A, Matrix<F>& B, Matrix<F>& C, Matrix<F>& D, 
  Matrix<F>& X, bool computeResidual=false );
template<typename F>
void LSE
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, 
  AbstractDistMatrix<F>& C, AbstractDistMatrix<F>& D, 
  AbstractDistMatrix<F>& X, bool computeResidual=false );

// TODO: Sparse-direct implementations

// Ridge regression
// ================
// A special case of Tikhonov regularization where the regularization matrix
// G is of the form gamma I.

namespace RidgeAlgNS {
enum RidgeAlg {
    RIDGE_CHOLESKY,
    RIDGE_QR,
    RIDGE_SVD
};
}
using namespace RidgeAlgNS;

template<typename F>
void Ridge
( Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& B, 
        Base<F> gamma,      Matrix<F>& X, 
  RidgeAlg alg=RIDGE_CHOLESKY );
template<typename F>
void Ridge
( Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, 
  Base<F> gamma,                        AbstractDistMatrix<F>& X, 
  RidgeAlg alg=RIDGE_CHOLESKY );

template<typename F>
void Ridge
( Orientation orientation,
  const SparseMatrix<F>& A, const Matrix<F>& B, 
        Base<F> gamma,            Matrix<F>& X, 
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );
template<typename F>
void Ridge
( Orientation orientation,
  const DistSparseMatrix<F>& A, const DistMultiVec<F>& B, 
        Base<F> gamma,                DistMultiVec<F>& X, 
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );

// Tikhonov regularization
// =======================

// Either solve the regularized Least Squares problem
//
//    min_ X || [W;G] X - [B;0] ||_F
// 
// or the regularized Minimum Length problem
//
//    min_{X,S} || [X; S] ||_F 
//    s.t. [W, G] [X; S] = B,
//
// where W is defined as op(A), which is either A, A^T, or A^H.
//
namespace TikhonovAlgNS {
enum TikhonovAlg {
    TIKHONOV_CHOLESKY,
    TIKHONOV_QR
};
}
using namespace TikhonovAlgNS;

template<typename F>
void Tikhonov
( Orientation orientation,
  const Matrix<F>& A, const Matrix<F>& B, 
  const Matrix<F>& G,       Matrix<F>& X, 
  TikhonovAlg alg=TIKHONOV_CHOLESKY );
template<typename F>
void Tikhonov
( Orientation orientation,
  const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, 
  const AbstractDistMatrix<F>& G,       AbstractDistMatrix<F>& X, 
  TikhonovAlg alg=TIKHONOV_CHOLESKY );

template<typename F>
void Tikhonov
( Orientation orientation,
  const SparseMatrix<F>& A, const Matrix<F>& B,
  const SparseMatrix<F>& G,       Matrix<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );
template<typename F>
void Tikhonov
( Orientation orientation,
  const DistSparseMatrix<F>& A, const DistMultiVec<F>& B,
  const DistSparseMatrix<F>& G,       DistMultiVec<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );

// TODO: Generalized Tikhonov regularization

// TODO: Total Least Squares

} // namespace El

#endif // ifndef EL_EUCLIDEANMIN_HPP
