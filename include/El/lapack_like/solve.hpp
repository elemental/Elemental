/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_SOLVE_HPP
#define EL_SOLVE_HPP

#include "El/lapack_like/factor.hpp"

namespace El {

// Linear solvers
// ==============

// B := A \ B for a general square A
// ---------------------------------
template<typename F>
void LinearSolve( Matrix<F>& A, Matrix<F>& B );
template<typename F>
void LinearSolve( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );
template<typename F>
void LinearSolve( const SparseMatrix<F>& A, Matrix<F>& B );
template<typename F>
void LinearSolve( const DistSparseMatrix<F>& A, DistMultiVec<F>& B );

// B := A \ B for a Hermitian A
// ----------------------------
template<typename F>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation, 
  Matrix<F>& A, Matrix<F>& B, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation,
  AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

template<typename F>
void HermitianSolve
( const SparseMatrix<F>& A, Matrix<F>& X,
  const BisectCtrl& ctrl=BisectCtrl() );
template<typename F>
void HermitianSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X,
  const BisectCtrl& ctrl=BisectCtrl() );

// B := B \ A for a symmetric A
// ----------------------------
template<typename F>
void SymmetricSolve
( UpperOrLower uplo, Orientation orientation, 
  Matrix<F>& A, Matrix<F>& B, bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );
template<typename F>
void SymmetricSolve
( UpperOrLower uplo, Orientation orientation,
  AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B, bool conjugate=false, 
  const LDLPivotCtrl<Base<F>>& ctrl=LDLPivotCtrl<Base<F>>() );

template<typename F>
void SymmetricSolve
( const SparseMatrix<F>& A, Matrix<F>& X, bool conjugate=false,
  const BisectCtrl& ctrl=BisectCtrl() );
template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A, DistMultiVec<F>& X, bool conjugate=false,
  const BisectCtrl& ctrl=BisectCtrl() );

// B := A \ B for an HPD A
// -----------------------
template<typename F>
void HPDSolve
( UpperOrLower uplo, Orientation orientation, 
  Matrix<F>& A, Matrix<F>& B );
template<typename F>
void HPDSolve
( UpperOrLower uplo, Orientation orientation,
  AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B );

// Solve for B in op(H) B - B op(D) = X, where H is Hessenberg
// -----------------------------------------------------------
template<typename F>
void MultiShiftHessSolve
( UpperOrLower uplo, Orientation orientation,
  F alpha, const Matrix<F>& H, const Matrix<F>& shifts,
  Matrix<F>& X );
template<typename F>
void MultiShiftHessSolve
( UpperOrLower uplo, Orientation orientation,
  F alpha, const AbstractDistMatrix<F>& H, const AbstractDistMatrix<F>& shifts,
  AbstractDistMatrix<F>& X );

} // namespace El

#endif // ifndef EL_SOLVE_HPP
