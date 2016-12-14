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
template<typename Field>
void LeastSquares
( Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& B,
        Matrix<Field>& X );
template<typename Field>
void LeastSquares
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Field>& X );

template<typename Real>
struct SQSDCtrl
{
    bool scaleTwoNorm=true;
    Int basisSize=15; // only used if 'scaleTwoNorm' is true

    bool canOverwrite=false;

    // TODO(poulson): Means of rescaling F similar to 'alpha' in
    // LeastSquaresCtrl. If F was the identity matrix and G was the zero matrix,
    // then the condition number could be potentially be unnecessarily squared.

    // Temporary and permanent regularization for the first, positive block
    Real reg0Tmp = Pow(limits::Epsilon<Real>(),Real(0.25));
    Real reg0Perm = Pow(limits::Epsilon<Real>(),Real(0.4));

    // Temporary and permanent regularization for the second, negative block
    Real reg1Tmp = Pow(limits::Epsilon<Real>(),Real(0.25));
    Real reg1Perm = Pow(limits::Epsilon<Real>(),Real(0.4));

    RegSolveCtrl<Real> solveCtrl;
    bool equilibrate=true;
    bool progress=false;
    bool time=false;
};

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

    SQSDCtrl<Real> sqsdCtrl;

    bool equilibrate=true;
    bool progress=false;
    bool time=false;

    LeastSquaresCtrl()
    {
        sqsdCtrl.scaleTwoNorm = false;
        sqsdCtrl.canOverwrite = true;
        sqsdCtrl.equilibrate = false;
    }
};

template<typename Field>
void LeastSquares
( Orientation orientation,
  const SparseMatrix<Field>& A,
  const Matrix<Field>& Y,
        Matrix<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );
template<typename Field>
void LeastSquares
( Orientation orientation,
  const DistSparseMatrix<Field>& A,
  const DistMultiVec<Field>& Y,
        DistMultiVec<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );

// Dense versions which overwrite their input
// ------------------------------------------
namespace ls {

template<typename Field>
void Overwrite
( Orientation orientation,
        Matrix<Field>& A,
  const Matrix<Field>& B,
        Matrix<Field>& X );
template<typename Field>
void Overwrite
( Orientation orientation,
        AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& B,
        AbstractDistMatrix<Field>& X );

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

template<typename Field>
void Ridge
( Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& B,
        Base<Field> gamma,
        Matrix<Field>& X,
  RidgeAlg alg=RIDGE_CHOLESKY );
template<typename Field>
void Ridge
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& B,
  Base<Field> gamma,
        AbstractDistMatrix<Field>& X,
  RidgeAlg alg=RIDGE_CHOLESKY );

template<typename Field>
void Ridge
( Orientation orientation,
  const SparseMatrix<Field>& A,
  const Matrix<Field>& B,
        Base<Field> gamma,
        Matrix<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );
template<typename Field>
void Ridge
( Orientation orientation,
  const DistSparseMatrix<Field>& A,
  const DistMultiVec<Field>& B,
        Base<Field> gamma,
        DistMultiVec<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );

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

template<typename Field>
void Tikhonov
( Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& B,
  const Matrix<Field>& G,
        Matrix<Field>& X,
  TikhonovAlg alg=TIKHONOV_CHOLESKY );
template<typename Field>
void Tikhonov
( Orientation orientation,
  const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& B,
  const AbstractDistMatrix<Field>& G,
        AbstractDistMatrix<Field>& X,
  TikhonovAlg alg=TIKHONOV_CHOLESKY );

template<typename Field>
void Tikhonov
( Orientation orientation,
  const SparseMatrix<Field>& A,
  const Matrix<Field>& B,
  const SparseMatrix<Field>& G,
        Matrix<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );
template<typename Field>
void Tikhonov
( Orientation orientation,
  const DistSparseMatrix<Field>& A,
  const DistMultiVec<Field>& B,
  const DistSparseMatrix<Field>& G,
        DistMultiVec<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );

// Equality-constrained Least Squarees
// ===================================
// Solve
//   min_X || A X - C ||_F subject to B X = D
template<typename Field>
void LSE
( const Matrix<Field>& A,
  const Matrix<Field>& B,
  const Matrix<Field>& C,
  const Matrix<Field>& D,
        Matrix<Field>& X );
template<typename Field>
void LSE
( const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& B,
  const AbstractDistMatrix<Field>& C,
  const AbstractDistMatrix<Field>& D,
        AbstractDistMatrix<Field>& X );

template<typename Field>
void LSE
( const SparseMatrix<Field>& A,
  const SparseMatrix<Field>& B,
  const Matrix<Field>& C,
  const Matrix<Field>& D,
        Matrix<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );
template<typename Field>
void LSE
( const DistSparseMatrix<Field>& A,
  const DistSparseMatrix<Field>& B,
  const DistMultiVec<Field>& C,
  const DistMultiVec<Field>& D,
        DistMultiVec<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );

// Dense versions which overwrite inputs where possible
// ----------------------------------------------------
namespace lse {

// A and B are overwritten with their factorizations and both C and D
// are modified
template<typename Field>
void Overwrite
( Matrix<Field>& A,
  Matrix<Field>& B,
  Matrix<Field>& C,
  Matrix<Field>& D,
  Matrix<Field>& X,
  bool computeResidual=false );
template<typename Field>
void Overwrite
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& B,
  AbstractDistMatrix<Field>& C,
  AbstractDistMatrix<Field>& D,
  AbstractDistMatrix<Field>& X,
  bool computeResidual=false );

} // namespace lse

// General (Gauss-Markov) Linear Model
// ===================================
// Solve
//   min_{X,Y} || Y ||_F subject to A X + B Y = D

template<typename Field>
void GLM
( const Matrix<Field>& A,
  const Matrix<Field>& B,
  const Matrix<Field>& D,
        Matrix<Field>& X,
        Matrix<Field>& Y );
template<typename Field>
void GLM
( const AbstractDistMatrix<Field>& A,
  const AbstractDistMatrix<Field>& B,
  const AbstractDistMatrix<Field>& D,
        AbstractDistMatrix<Field>& X,
        AbstractDistMatrix<Field>& Y );

template<typename Field>
void GLM
( const SparseMatrix<Field>& A,
  const SparseMatrix<Field>& B,
  const Matrix<Field>& D,
        Matrix<Field>& X,
        Matrix<Field>& Y,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );
template<typename Field>
void GLM
( const DistSparseMatrix<Field>& A,
  const DistSparseMatrix<Field>& B,
  const DistMultiVec<Field>& D,
        DistMultiVec<Field>& X,
        DistMultiVec<Field>& Y,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );

// Dense versions which overwrite the input where possible
// -------------------------------------------------------
namespace glm {

// A and B are overwritten with their factorizations and X is returned in D
template<typename Field>
void Overwrite
( Matrix<Field>& A,
  Matrix<Field>& B,
  Matrix<Field>& D,
  Matrix<Field>& Y );
template<typename Field>
void Overwrite
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& B,
  AbstractDistMatrix<Field>& D,
  AbstractDistMatrix<Field>& Y );

} // namespace glm

// TODO: Generalized Tikhonov regularization

// TODO: Total Least Squares

} // namespace El

#endif // ifndef EL_EUCLIDEANMIN_HPP
