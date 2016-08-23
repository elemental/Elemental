/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_EUCLIDEANMIN_HPP
#define EL_EUCLIDEANMIN_HPP

#include <El/lapack_like/factor.hpp>

namespace El {

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
  const Matrix<F>& A,
  const Matrix<F>& B, 
        Matrix<F>& X );
template<typename F>
void LeastSquares
( Orientation orientation, 
  const ElementalMatrix<F>& A,
  const ElementalMatrix<F>& B,
        ElementalMatrix<F>& X );

template<typename Real>
struct LeastSquaresCtrl 
{
    bool scaleTwoNorm=true;
    Int basisSize=15; // only used if 'scaleTwoNorm' is true

    // Note: while 'alpha' should ideally be roughly equal to the minimum 
    //       singular value of A (possibly scaled down to unit two-norm),
    //       Saunders has recommended in at least one publication that 
    //       alpha ~= 1e-4 is a decent default. After experimenting with the
    //       estimation of the minimum singular value via Lanczos on A^H A
    //       failed (the literature agrees), I fell back to this default value.
    Real alpha=Pow(limits::Epsilon<Real>(),Real(0.25));

    // Temporary and permanent regularization for the first, positive block of
    // the augmented system
    Real reg0Tmp = Pow(limits::Epsilon<Real>(),Real(0.25));
    Real reg0Perm = Pow(limits::Epsilon<Real>(),Real(0.4));

    // Temporary and permanent regularization for the second, negative block of
    // the augmented system
    Real reg1Tmp = Pow(limits::Epsilon<Real>(),Real(0.25));
    Real reg1Perm = Pow(limits::Epsilon<Real>(),Real(0.4));

    RegSolveCtrl<Real> solveCtrl;
    bool equilibrate=true;
    bool progress=false;
    bool time=false;
};

template<typename F>
void LeastSquares
( Orientation orientation,
  const SparseMatrix<F>& A,
  const Matrix<F>& Y, 
        Matrix<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );
template<typename F>
void LeastSquares
( Orientation orientation,
  const DistSparseMatrix<F>& A,
  const DistMultiVec<F>& Y, 
        DistMultiVec<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );

// Dense versions which overwrite their input
// ------------------------------------------
namespace ls {

template<typename F>
void Overwrite
( Orientation orientation, 
  Matrix<F>& A, const Matrix<F>& B, 
                      Matrix<F>& X );
template<typename F>
void Overwrite
( Orientation orientation, 
  ElementalMatrix<F>& A, const ElementalMatrix<F>& B, 
                                  ElementalMatrix<F>& X );

} // namespace ls

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
  const ElementalMatrix<F>& A, const ElementalMatrix<F>& B, 
  Base<F> gamma,                     ElementalMatrix<F>& X, 
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
  const ElementalMatrix<F>& A, const ElementalMatrix<F>& B, 
  const ElementalMatrix<F>& G,       ElementalMatrix<F>& X, 
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

// Equality-constrained Least Squarees
// ===================================
// Solve
//   min_X || A X - C ||_F subject to B X = D
template<typename F>
void LSE
( const Matrix<F>& A, const Matrix<F>& B, 
  const Matrix<F>& C, const Matrix<F>& D, 
        Matrix<F>& X );
template<typename F>
void LSE
( const ElementalMatrix<F>& A, const ElementalMatrix<F>& B, 
  const ElementalMatrix<F>& C, const ElementalMatrix<F>& D, 
        ElementalMatrix<F>& X );

template<typename F>
void LSE
( const SparseMatrix<F>& A, const SparseMatrix<F>& B,
  const Matrix<F>& C,       const Matrix<F>& D,
        Matrix<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );
template<typename F>
void LSE
( const DistSparseMatrix<F>& A, const DistSparseMatrix<F>& B,
  const DistMultiVec<F>& C,     const DistMultiVec<F>& D,
        DistMultiVec<F>& X,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );

// Dense versions which overwrite inputs where possible
// ----------------------------------------------------
namespace lse {

// A and B are overwritten with their factorizations and both C and D
// are modified
template<typename F>
void Overwrite
( Matrix<F>& A, Matrix<F>& B, Matrix<F>& C, Matrix<F>& D, 
  Matrix<F>& X, bool computeResidual=false );
template<typename F>
void Overwrite
( ElementalMatrix<F>& A, ElementalMatrix<F>& B, 
  ElementalMatrix<F>& C, ElementalMatrix<F>& D, 
  ElementalMatrix<F>& X, bool computeResidual=false );

} // namespace lse

// Generalized (Gauss-Markov) Linear Model
// =======================================
// Solve 
//   min_{X,Y} || Y ||_F subject to A X + B Y = D

template<typename F>
void GLM
( const Matrix<F>& A, const Matrix<F>& B, 
  const Matrix<F>& D, 
        Matrix<F>& X,       Matrix<F>& Y );
template<typename F>
void GLM
( const ElementalMatrix<F>& A, const ElementalMatrix<F>& B, 
  const ElementalMatrix<F>& D, 
        ElementalMatrix<F>& X,       ElementalMatrix<F>& Y );

template<typename F>
void GLM
( const SparseMatrix<F>& A, const SparseMatrix<F>& B,
  const Matrix<F>& D,             
        Matrix<F>& X,             Matrix<F>& Y,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );
template<typename F>
void GLM
( const DistSparseMatrix<F>& A, const DistSparseMatrix<F>& B,
  const DistMultiVec<F>& D,           
        DistMultiVec<F>& X,           DistMultiVec<F>& Y,
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );

// Dense versions which overwrite the input where possible
// -------------------------------------------------------
namespace glm {

// A and B are overwritten with their factorizations and X is returned in D
template<typename F>
void Overwrite
( Matrix<F>& A, Matrix<F>& B, Matrix<F>& D, Matrix<F>& Y );
template<typename F>
void Overwrite
( ElementalMatrix<F>& A, ElementalMatrix<F>& B, 
  ElementalMatrix<F>& D, ElementalMatrix<F>& Y );

} // namespace glm

// TODO: Generalized Tikhonov regularization

// TODO: Total Least Squares

} // namespace El

#endif // ifndef EL_EUCLIDEANMIN_HPP
