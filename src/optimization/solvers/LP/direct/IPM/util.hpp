/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "../../../QP/direct/IPM/util.hpp"

namespace El {
namespace lp {
namespace direct {

// Initialize
// ==========
template<typename Real>
void Initialize
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y, 
        Matrix<Real>& z,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& b,
  const ElementalMatrix<Real>& c,
        ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& y, 
        ElementalMatrix<Real>& z,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y, 
        Matrix<Real>& z,
        vector<Int>& map,
        vector<Int>& invMap,
        ldl::Separator& rootSep,
        ldl::NodeInfo& info,
  bool primalInit, bool dualInit, bool standardShift, 
  const RegSolveCtrl<Real>& solveCtrl );
template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z, 
        DistMap& map,
        DistMap& invMap,
        ldl::DistSeparator& rootSep,
        ldl::DistNodeInfo& info,
        vector<Int>& mappedSources, 
        vector<Int>& mappedTargets,
        vector<Int>& colOffs,
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
( const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& J,
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
( const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& J,
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
( const ElementalMatrix<Real>& A,
        Real gamma,
        Real delta,
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& J,
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
( const ElementalMatrix<Real>& A,
        Real gamma,
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
  const ElementalMatrix<Real>& rc,
  const ElementalMatrix<Real>& rb, 
  const ElementalMatrix<Real>& rmu,
        ElementalMatrix<Real>& d );
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
( const Matrix<Real>& A,
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
( const ElementalMatrix<Real>& A,
        Real gamma,
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
  const ElementalMatrix<Real>& rc,
  const ElementalMatrix<Real>& rmu,
        ElementalMatrix<Real>& dx, 
  const ElementalMatrix<Real>& dy, 
        ElementalMatrix<Real>& dz );
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
