/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "../../../QP/direct/IPM/util.hpp"

namespace El {
namespace lp {
namespace direct {

// Initialize
// ==========
template<typename Real>
void Initialize
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
        SparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit, bool dualInit, bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl );
template<typename Real>
void Initialize
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
        DistSparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit, bool dualInit, bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl );

// Full system
// ===========
template<typename Real>
void KKT
( const Matrix<Real>& A,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const DistMatrix<Real>& A,
  const DistMatrix<Real>& x,
  const DistMatrix<Real>& z,
        DistMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const SparseMatrix<Real>& A,
        Real gamma,
        Real delta,
        Real beta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A,
        Real gamma,
        Real delta,
        Real beta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower=true );

using qp::direct::KKTRHS;
using qp::direct::ExpandSolution;

// Augmented system
// ================
template<typename Real>
void AugmentedKKT
( const Matrix<Real>& A,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const DistMatrix<Real>& A,
  const DistMatrix<Real>& x,
  const DistMatrix<Real>& z,
        DistMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& A,
        Real gamma,
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& A,
        Real gamma,
        Real delta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower=true );

using qp::direct::AugmentedKKTRHS;
using qp::direct::ExpandAugmentedSolution;

// Normal system
// =============
template<typename Real>
void NormalKKT
( const Matrix<Real>& A,
        Real gamma,
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower=false );
template<typename Real>
void NormalKKT
( const DistMatrix<Real>& A,
        Real gamma,
        Real delta,
  const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& z,
        DistMatrix<Real>& J,
  bool onlyLower=false );
template<typename Real>
void NormalKKT
( const SparseMatrix<Real>& A,
        Real gamma,
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void NormalKKT
( const DistSparseMatrix<Real>& A,
        Real gamma,
        Real delta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower=true );

template<typename Real>
void NormalKKTRHS
( const Matrix<Real>& A,
        Real gamma,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
  const Matrix<Real>& rc,
  const Matrix<Real>& rb,
  const Matrix<Real>& rmu,
        Matrix<Real>& d );
template<typename Real>
void NormalKKTRHS
( const DistMatrix<Real>& A,
        Real gamma,
  const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& z,
  const DistMatrix<Real>& rc,
  const DistMatrix<Real>& rb,
  const DistMatrix<Real>& rmu,
        DistMatrix<Real>& d );
template<typename Real>
void NormalKKTRHS
( const SparseMatrix<Real>& A,
        Real gamma,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
  const Matrix<Real>& rc,
  const Matrix<Real>& rb,
  const Matrix<Real>& rmu,
        Matrix<Real>& d );
template<typename Real>
void NormalKKTRHS
( const DistSparseMatrix<Real>& A,
        Real gamma,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb,
  const DistMultiVec<Real>& rmu,
        DistMultiVec<Real>& d );

template<typename Real>
void ExpandNormalSolution
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        Real gamma,
  const DirectLPSolution<Matrix<Real>>& solution,
  const DirectLPResidual<Matrix<Real>>& residual,
        DirectLPSolution<Matrix<Real>>& correction );
template<typename Real>
void ExpandNormalSolution
( const DistMatrix<Real>& A,
        Real gamma,
  const AbstractDistMatrix<Real>& x,
  const AbstractDistMatrix<Real>& z,
  const DistMatrix<Real>& rc,
  const DistMatrix<Real>& rmu,
        DistMatrix<Real>& dx,
  const DistMatrix<Real>& dy,
        DistMatrix<Real>& dz );
template<typename Real>
void ExpandNormalSolution
( const SparseMatrix<Real>& A,
        Real gamma,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
  const Matrix<Real>& rc,
  const Matrix<Real>& rmu,
        Matrix<Real>& dx,
  const Matrix<Real>& dy,
        Matrix<Real>& dz );
template<typename Real>
void ExpandNormalSolution
( const DistSparseMatrix<Real>& A,
        Real gamma,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rmu,
        DistMultiVec<Real>& dx,
  const DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz );

} // namespace direct
} // namespace lp
} // namespace El
