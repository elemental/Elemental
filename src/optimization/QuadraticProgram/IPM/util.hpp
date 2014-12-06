/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace quad_prog {

// Full system
// -----------
template<typename Real>
void KKT
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& s, const Matrix<Real>& x,
        Matrix<Real>& J );
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x,
  AbstractDistMatrix<Real>& J );

template<typename Real>
void KKTRHS
( const Matrix<Real>& rmu, const Matrix<Real>& rc,
  const Matrix<Real>& rb, Matrix<Real>& y );
template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc,
  const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& y );

template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const Matrix<Real>& y,
  Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl );
template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const AbstractDistMatrix<Real>& yPre,
  AbstractDistMatrix<Real>& ds, AbstractDistMatrix<Real>& dx, 
  AbstractDistMatrix<Real>& dl );

// Augmented system
// ----------------
template<typename Real>
void AugmentedKKT
( const Matrix<Real>& Q, const Matrix<Real>& A,
  const Matrix<Real>& s, const Matrix<Real>& x,
  Matrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x,
  AbstractDistMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A,
  const Matrix<Real>& s, const Matrix<Real>& x,
  SparseMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& s, const DistMultiVec<Real>& x,
  DistSparseMatrix<Real>& J, bool onlyLower=true );

template<typename Real>
void AugmentedKKTRHS
( const Matrix<Real>& x,
  const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb,
  Matrix<Real>& y );
template<typename Real>
void AugmentedKKTRHS
( const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, 
  const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& y );
template<typename Real>
void AugmentedKKTRHS
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, 
  const DistMultiVec<Real>& rb, DistMultiVec<Real>& y );

template<typename Real>
void ExpandAugmentedSolution
( const Matrix<Real>& s, const Matrix<Real>& x,
  const Matrix<Real>& rmu, const Matrix<Real>& y,
  Matrix<Real>& ds, Matrix<Real>& dx, Matrix<Real>& dl );
template<typename Real>
void ExpandAugmentedSolution
( const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& y,
  AbstractDistMatrix<Real>& ds, AbstractDistMatrix<Real>& dx, 
  AbstractDistMatrix<Real>& dl );
template<typename Real>
void ExpandAugmentedSolution
( const DistMultiVec<Real>& s, const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& y,
  DistMultiVec<Real>& ds, DistMultiVec<Real>& dx, DistMultiVec<Real>& dl );

// Line search
// -----------
template<typename Real>
Real IPFLineSearch
( const Matrix<Real>& Q,  const Matrix<Real>& A,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& s,  const Matrix<Real>& x,  const Matrix<Real>& l,
  const Matrix<Real>& ds, const Matrix<Real>& dx, const Matrix<Real>& dl,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const AbstractDistMatrix<Real>& Q,  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& s,  const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& l,
  const AbstractDistMatrix<Real>& ds, const AbstractDistMatrix<Real>& dx,
  const AbstractDistMatrix<Real>& dl,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& s,  const Matrix<Real>& x,  const Matrix<Real>& l,
  const Matrix<Real>& ds, const Matrix<Real>& dx, const Matrix<Real>& dl,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const DistSparseMatrix<Real>& Q, 
  const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& l,
  const DistMultiVec<Real>& ds,
  const DistMultiVec<Real>& dx,
  const DistMultiVec<Real>& dl,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );

} // namespace quad_prog
} // namespace El
