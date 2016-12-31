/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "../../../QP/affine/IPM/util.hpp"

namespace El {
namespace lp {
namespace affine {

// Initialize
// ==========
template<typename Real>
void Initialize
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const SparseMatrix<Real>& JStatic,
  const Matrix<Real>& regTmp,
        SparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit, bool dualInit, bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl );
template<typename Real>
void Initialize
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const DistSparseMatrix<Real>& JStatic,
  const DistMultiVec<Real>& regTmp,
        DistSparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit, bool dualInit, bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl );

// Full system
// ===========
// We explicitly accept s and z rather than AffineLPSolution because it is
// common to initialize IPM's by solving the KKT system with s = z = ones(k,1).
template<typename Real>
void KKT
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const DistMatrix<Real>& A,
  const DistMatrix<Real>& G,
  const DistMatrix<Real>& s,
  const DistMatrix<Real>& z,
        DistMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void StaticKKT
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        SparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void StaticKKT
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        DistSparseMatrix<Real>& J,
  bool onlyLower=true );

using qp::affine::FinishKKT;

using qp::affine::KKTRHS;
using qp::affine::ExpandCoreSolution;
using qp::affine::ExpandSolution;

} // namespace affine
} // namespace lp
} // namespace El
