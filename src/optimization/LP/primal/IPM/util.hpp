/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace lp {
namespace primal {

// Full system
// -----------
template<typename Real>
void KKT
( const Matrix<Real>& A, 
  const Matrix<Real>& x, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  AbstractDistMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& x, const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower=true );

template<typename Real>
void KKTRHS
( const Matrix<Real>& rmu, const Matrix<Real>& rc,
  const Matrix<Real>& rb, const Matrix<Real>& z,
        Matrix<Real>& d );
template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc,
  const AbstractDistMatrix<Real>& rb, const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& d );
template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb, const DistMultiVec<Real>& z,
        DistMultiVec<Real>& d );

template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const Matrix<Real>& d,
  Matrix<Real>& dx, Matrix<Real>& dy, Matrix<Real>& dz );
template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const AbstractDistMatrix<Real>& d,
  AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, 
  AbstractDistMatrix<Real>& dz );
template<typename Real>
void ExpandKKTSolution
( Int m, Int n, const DistMultiVec<Real>& d,
  DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, 
  DistMultiVec<Real>& dz );

// Augmented system
// ----------------
template<typename Real>
void AugmentedKKT
( const Matrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
  Matrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  AbstractDistMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
  SparseMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
  DistSparseMatrix<Real>& J, bool onlyLower=true );

template<typename Real>
void AugmentedKKTRHS
( const Matrix<Real>& x,
  const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb,
  Matrix<Real>& d );
template<typename Real>
void AugmentedKKTRHS
( const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, 
  const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& d );
template<typename Real>
void AugmentedKKTRHS
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc, 
  const DistMultiVec<Real>& rb, DistMultiVec<Real>& d );

template<typename Real>
void ExpandAugmentedSolution
( const Matrix<Real>& x, const Matrix<Real>& z,
  const Matrix<Real>& rmu, const Matrix<Real>& d,
  Matrix<Real>& dx, Matrix<Real>& dy, Matrix<Real>& dz );
template<typename Real>
void ExpandAugmentedSolution
( const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& d,
  AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, 
  AbstractDistMatrix<Real>& dz );
template<typename Real>
void ExpandAugmentedSolution
( const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& d,
  DistMultiVec<Real>& dx, DistMultiVec<Real>& dy, DistMultiVec<Real>& dz );

// Normal system
// -------------
template<typename Real>
void NormalKKT
( const Matrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
  Matrix<Real>& J, bool onlyLower=true );
template<typename Real>
void NormalKKT
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  AbstractDistMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void NormalKKT
( const SparseMatrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
  SparseMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void NormalKKT
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
  DistSparseMatrix<Real>& J, bool onlyLower=true );

template<typename Real>
void NormalKKTRHS
( const Matrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
  const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb,
        Matrix<Real>& d );
template<typename Real>
void NormalKKTRHS
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& x, const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc, 
  const AbstractDistMatrix<Real>& rb, AbstractDistMatrix<Real>& d );
template<typename Real>
void NormalKKTRHS
( const SparseMatrix<Real>& A,
  const Matrix<Real>& x, const Matrix<Real>& z,
  const Matrix<Real>& rmu, const Matrix<Real>& rc, const Matrix<Real>& rb,
        Matrix<Real>& d );
template<typename Real>
void NormalKKTRHS
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& x, const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rmu, const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb,
        DistMultiVec<Real>& d );

template<typename Real>
void ExpandNormalSolution
( const Matrix<Real>& A,   const Matrix<Real>& c,
  const Matrix<Real>& x,   const Matrix<Real>& z,
  const Matrix<Real>& rmu, const Matrix<Real>& rc,
        Matrix<Real>& dx, 
  const Matrix<Real>& dy, 
        Matrix<Real>& dz );
template<typename Real>
void ExpandNormalSolution
( const AbstractDistMatrix<Real>& A,   const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& x,   const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& rmu, const AbstractDistMatrix<Real>& rc,
        AbstractDistMatrix<Real>& dx, 
  const AbstractDistMatrix<Real>& dy, 
        AbstractDistMatrix<Real>& dz );
template<typename Real>
void ExpandNormalSolution
( const SparseMatrix<Real>& A, const Matrix<Real>& c,
  const Matrix<Real>& x,       const Matrix<Real>& z,
  const Matrix<Real>& rmu,     const Matrix<Real>& rc,
        Matrix<Real>& dx, 
  const Matrix<Real>& dy, 
        Matrix<Real>& dz );
template<typename Real>
void ExpandNormalSolution
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rmu,   const DistMultiVec<Real>& rc,
        DistMultiVec<Real>& dx,
  const DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz );

// Line search
// -----------
template<typename Real>
Real IPFLineSearch
( const Matrix<Real>& A,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& x,  const Matrix<Real>& y,  const Matrix<Real>& z,
  const Matrix<Real>& dx, const Matrix<Real>& dy, const Matrix<Real>& dz,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& dx, const AbstractDistMatrix<Real>& dy,
  const AbstractDistMatrix<Real>& dz,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& x,  const Matrix<Real>& y,  const Matrix<Real>& z,
  const Matrix<Real>& dx, const Matrix<Real>& dy, const Matrix<Real>& dz,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x,  const DistMultiVec<Real>& y,
  const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& dx, const DistMultiVec<Real>& dy,
  const DistMultiVec<Real>& dz,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );

} // namespace primal
} // namespace lp
} // namespace El
