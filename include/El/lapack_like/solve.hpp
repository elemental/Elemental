/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SOLVE_HPP
#define EL_SOLVE_HPP

#include <El/lapack_like/factor.hpp>
#include <El/lapack_like/euclidean_min.hpp>

namespace El {

// Linear
// ======
template<typename Field>
void LinearSolve
( const Matrix<Field>& A,
        Matrix<Field>& B );
template<typename Field>
void LinearSolve
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
  bool scalapack=false );

template<typename Field>
void LinearSolve
( const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );
template<typename Field>
void LinearSolve
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  const LeastSquaresCtrl<Base<Field>>& ctrl=LeastSquaresCtrl<Base<Field>>() );

namespace lin_solve {

template<typename Field>
void Overwrite( Matrix<Field>& A, Matrix<Field>& B );
template<typename Field>
void Overwrite( AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B );

} // namespace lin_solve

// Hermitian
// =========
template<typename Field>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation,
  const Matrix<Field>& A,
        Matrix<Field>& B,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );

template<typename Field>
void HermitianSolve
( const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  bool tryLDL=false,
  const BisectCtrl& ctrl=BisectCtrl() );
template<typename Field>
void HermitianSolve
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  bool tryLDL=false,
  const BisectCtrl& ctrl=BisectCtrl() );

namespace herm_solve {

template<typename Field>
void Overwrite
( UpperOrLower uplo, Orientation orientation,
  Matrix<Field>& A,
  Matrix<Field>& B,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void Overwrite
( UpperOrLower uplo, Orientation orientation,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& B,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );

} // namespace herm_solve

// Symmetric
// =========
template<typename Field>
void SymmetricSolve
( UpperOrLower uplo,
  Orientation orientation,
  const Matrix<Field>& A,
        Matrix<Field>& B,
  bool conjugate=false,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void SymmetricSolve
( UpperOrLower uplo,
  Orientation orientation,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
  bool conjugate=false,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );

template<typename Field>
void SymmetricSolve
( const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  bool conjugate=false,
  bool tryLDL=false,
  const BisectCtrl& ctrl=BisectCtrl() );
template<typename Field>
void SymmetricSolve
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  bool conjugate=false,
  bool tryLDL=false,
  const BisectCtrl& ctrl=BisectCtrl() );

namespace symm_solve {

template<typename Field>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation,
  Matrix<Field>& A,
  Matrix<Field>& B,
  bool conjugate=false,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );
template<typename Field>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& B,
  bool conjugate=false,
  const LDLPivotCtrl<Base<Field>>& ctrl=LDLPivotCtrl<Base<Field>>() );

} // namespace symm_solve

// Symmetric Quasi Semi-Definite (SQSD)
// ====================================
// Solve J X = B, where
//
//   J = | F,    A |,
//       | A^T, -G |
//
// and F and G are Symmetric Positive Semi-Definite (and F is n0 x n0).
//
// Dense matrices need only be explicitly filled in either their upper or
// lower triangle, but sparse matrices must be explicitly symmetric.
//
template<typename Field>
void SQSDSolve
( Int n0,
  UpperOrLower uplo,
  const Matrix<Field>& J,
        Matrix<Field>& B );
template<typename Field>
void SQSDSolve
( Int n0,
  UpperOrLower uplo,
  const AbstractDistMatrix<Field>& J,
        AbstractDistMatrix<Field>& B );

template<typename Field>
void SQSDSolve
( Int n0,
  const SparseMatrix<Field>& J,
        Matrix<Field>& B,
  const SQSDCtrl<Base<Field>>& ctrl=SQSDCtrl<Base<Field>>() );
template<typename Field>
void SQSDSolve
( Int n0,
  const DistSparseMatrix<Field>& J,
        DistMultiVec<Field>& B,
  const SQSDCtrl<Base<Field>>& ctrl=SQSDCtrl<Base<Field>>() );

// Hermitian Positive-Definite
// ===========================
template<typename Field>
void HPDSolve
( UpperOrLower uplo,
  Orientation orientation,
  const Matrix<Field>& A,
        Matrix<Field>& B );
template<typename Field>
void HPDSolve
( UpperOrLower uplo,
  Orientation orientation,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B );

template<typename Field>
void HPDSolve
( const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  const BisectCtrl& ctrl=BisectCtrl() );
template<typename Field>
void HPDSolve
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  const BisectCtrl& ctrl=BisectCtrl() );

namespace hpd_solve {

template<typename Field>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation,
  Matrix<Field>& A,
  Matrix<Field>& B );
template<typename Field>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& B );

} // namespace hpd_solve

// Multi-shift Hessenberg
// ======================
template<typename Field>
void MultiShiftHessSolve
( UpperOrLower uplo,
  Orientation orientation,
  Field alpha,
  const Matrix<Field>& H,
  const Matrix<Field>& shifts,
        Matrix<Field>& X );
template<typename Field>
void MultiShiftHessSolve
( UpperOrLower uplo,
  Orientation orientation,
  Field alpha,
  const AbstractDistMatrix<Field>& H,
  const AbstractDistMatrix<Field>& shifts,
        AbstractDistMatrix<Field>& X );

} // namespace El

#include <El/lapack_like/solve/FGMRES.hpp>
#include <El/lapack_like/solve/LGMRES.hpp>
#include <El/lapack_like/solve/Refined.hpp>

#endif // ifndef EL_SOLVE_HPP
