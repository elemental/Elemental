/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace qp {
namespace primal {

// Full system
// -----------
template<typename Real>
void KKT
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& x, const Matrix<Real>& z,
        Matrix<Real>& J );
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  AbstractDistMatrix<Real>& J );

template<typename Real>
void KKTRHS
( const Matrix<Real>& rmu, const Matrix<Real>& rc,
  const Matrix<Real>& rb, Matrix<Real>& rhs );
template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc,
  const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& rhs );

template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const Matrix<Real>& rhs,
  Matrix<Real>& dx, Matrix<Real>& dy, Matrix<Real>& dz );
template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const AbstractDistMatrix<Real>& rhs,
  AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, 
  AbstractDistMatrix<Real>& dz );

// Augmented system
// ----------------
template<typename Real>
void AugmentedKKT
( const Matrix<Real>& Q, const Matrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
  Matrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  AbstractDistMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
  SparseMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
  DistSparseMatrix<Real>& J, bool onlyLower=true );

template<typename Real>
void AugmentedKKTRHS
( const Matrix<Real>& x,
  const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb,
  Matrix<Real>& rhs );
template<typename Real>
void AugmentedKKTRHS
( const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, 
  const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& rhs );
template<typename Real>
void AugmentedKKTRHS
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, 
  const DistMultiVec<Real>& rb, DistMultiVec<Real>& rhs );

template<typename Real>
void ExpandAugmentedSolution
( const Matrix<Real>& x, const Matrix<Real>& z,
  const Matrix<Real>& rmu, const Matrix<Real>& rhs,
  Matrix<Real>& dx, Matrix<Real>& dy, Matrix<Real>& dz );
template<typename Real>
void ExpandAugmentedSolution
( const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rhs,
  AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, 
  AbstractDistMatrix<Real>& dz );
template<typename Real>
void ExpandAugmentedSolution
( const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rhs,
  DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, DistMultiVec<Real>& dz );

// Line search
// -----------
template<typename Real>
Real IPFLineSearch
( const Matrix<Real>& Q,  const Matrix<Real>& A,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& x,  const Matrix<Real>& y,  const Matrix<Real>& z,
  const Matrix<Real>& dx, const Matrix<Real>& dy, const Matrix<Real>& dz,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const AbstractDistMatrix<Real>& Q,  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& dx, const AbstractDistMatrix<Real>& dy,
  const AbstractDistMatrix<Real>& dz,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& x,  const Matrix<Real>& y,  const Matrix<Real>& z,
  const Matrix<Real>& dx, const Matrix<Real>& dy, const Matrix<Real>& dz,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const DistSparseMatrix<Real>& Q, 
  const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& y,
  const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& dx,
  const DistMultiVec<Real>& dy,
  const DistMultiVec<Real>& dz,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );

} // namespace primal
} // namespace qp
} // namespace El
