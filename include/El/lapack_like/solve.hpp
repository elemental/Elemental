/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SOLVE_HPP
#define EL_SOLVE_HPP

#include "El/lapack_like/factor.hpp"
#include "El/lapack_like/euclidean_min.hpp"

namespace El {

// Linear
// ======
template<typename F>
void LinearSolve( const Matrix<F>& A, Matrix<F>& B );
template<typename F>
void LinearSolve( const ElementalMatrix<F>& A, ElementalMatrix<F>& B );
template<typename F>
void LinearSolve
( const DistMatrix<F,MC,MR,BLOCK>& A, DistMatrix<F,MC,MR,BLOCK>& B );

template<typename F>
void LinearSolve
( const SparseMatrix<F>& A, Matrix<F>& B, 
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );
template<typename F>
void LinearSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& B, 
  const LeastSquaresCtrl<Base<F>>& ctrl=LeastSquaresCtrl<Base<F>>() );

namespace lin_solve {

template<typename F>
void Overwrite( Matrix<F>& A, Matrix<F>& B );
template<typename F>
void Overwrite( ElementalMatrix<F>& A, ElementalMatrix<F>& B );

} // namespace lin_solve

// Hermitian
// =========
template<typename F>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, Matrix<F>& B, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation,
  const ElementalMatrix<F>& A, ElementalMatrix<F>& B, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

template<typename F>
void HermitianSolve
( const SparseMatrix<F>& A, Matrix<F>& B,
  bool tryLDL=false,
  const BisectCtrl& ctrl=BisectCtrl() );
template<typename F>
void HermitianSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& B,
  bool tryLDL=false,
  const BisectCtrl& ctrl=BisectCtrl() );

namespace herm_solve {

template<typename F>
void Overwrite
( UpperOrLower uplo, Orientation orientation, 
  Matrix<F>& A, Matrix<F>& B, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void Overwrite
( UpperOrLower uplo, Orientation orientation,
  ElementalMatrix<F>& A, ElementalMatrix<F>& B, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

} // namespace herm_solve

// Symmetric
// =========
template<typename F>
void SymmetricSolve
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, Matrix<F>& B, 
  bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void SymmetricSolve
( UpperOrLower uplo, Orientation orientation,
  const ElementalMatrix<F>& A, ElementalMatrix<F>& B, 
  bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

template<typename F>
void SymmetricSolve
( const SparseMatrix<F>& A, Matrix<F>& B, 
  bool conjugate=false, bool tryLDL=false,
  const BisectCtrl& ctrl=BisectCtrl() );
template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& B, 
  bool conjugate=false, bool tryLDL=false,
  const BisectCtrl& ctrl=BisectCtrl() );

namespace symm_solve {

template<typename F>
void Overwrite
( UpperOrLower uplo, Orientation orientation, 
  Matrix<F>& A, Matrix<F>& B, 
  bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void Overwrite
( UpperOrLower uplo, Orientation orientation,
  ElementalMatrix<F>& A, ElementalMatrix<F>& B, 
  bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

} // namespace symm_solve

// Hermitian Positive-Definite
// ===========================
template<typename F>
void HPDSolve
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, Matrix<F>& B );
template<typename F>
void HPDSolve
( UpperOrLower uplo, Orientation orientation,
  const ElementalMatrix<F>& A, ElementalMatrix<F>& B );

template<typename F>
void HPDSolve
( const SparseMatrix<F>& A, Matrix<F>& B, 
  const BisectCtrl& ctrl=BisectCtrl() );
template<typename F>
void HPDSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& B, 
  const BisectCtrl& ctrl=BisectCtrl() );

namespace hpd_solve {

template<typename F>
void Overwrite
( UpperOrLower uplo, Orientation orientation, 
  Matrix<F>& A, Matrix<F>& B );
template<typename F>
void Overwrite
( UpperOrLower uplo, Orientation orientation,
  ElementalMatrix<F>& A, ElementalMatrix<F>& B );

} // namespace hpd_solve

// Multi-shift Hessenberg
// ======================
template<typename F>
void MultiShiftHessSolve
( UpperOrLower uplo, Orientation orientation,
  F alpha, const Matrix<F>& H, const Matrix<F>& shifts,
  Matrix<F>& X );
template<typename F>
void MultiShiftHessSolve
( UpperOrLower uplo, Orientation orientation,
  F alpha, const ElementalMatrix<F>& H, const ElementalMatrix<F>& shifts,
  ElementalMatrix<F>& X );

} // namespace El

#include "./solve/FGMRES.hpp"
#include "./solve/LGMRES.hpp"
#include "./solve/Refined.hpp"

#endif // ifndef EL_SOLVE_HPP
