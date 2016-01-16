/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace qp {
namespace direct {

// Initialize
// ==========
template<typename Real>
void Initialize
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y, 
        Matrix<Real>& z,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& b,
  const ElementalMatrix<Real>& c,
        ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& y, 
        ElementalMatrix<Real>& z,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
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
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
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
( const Matrix<Real>& Q,
  const Matrix<Real>& A, 
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Real>& A, 
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A, 
        Real gamma,
        Real delta,
        Real beta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A, 
        Real gamma,
        Real delta,
        Real beta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower=true );

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,
  const Matrix<Real>& rb, 
  const Matrix<Real>& rmu,
  const Matrix<Real>& z,
        Matrix<Real>& d );
template<typename Real>
void KKTRHS
( const ElementalMatrix<Real>& rc,
  const ElementalMatrix<Real>& rb, 
  const ElementalMatrix<Real>& rmu,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& d );
template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& d );

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const Matrix<Real>& d,
        Matrix<Real>& dx,
        Matrix<Real>& dy, 
        Matrix<Real>& dz );
template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const ElementalMatrix<Real>& d,
        ElementalMatrix<Real>& dx,
        ElementalMatrix<Real>& dy, 
        ElementalMatrix<Real>& dz );
template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const DistMultiVec<Real>& d,
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz );

// Augmented system
// ================
template<typename Real>
void AugmentedKKT
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const ElementalMatrix<Real>& Q,
  const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
        Real gamma,
        Real delta,
  const Matrix<Real>& x,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void AugmentedKKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
        Real gamma,
        Real delta,
  const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower=true );

template<typename Real>
void AugmentedKKTRHS
( const Matrix<Real>& x,
  const Matrix<Real>& rc,
  const Matrix<Real>& rb, 
  const Matrix<Real>& rmu,
        Matrix<Real>& d );
template<typename Real>
void AugmentedKKTRHS
( const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& rc,
  const ElementalMatrix<Real>& rb, 
  const ElementalMatrix<Real>& rmu,
        ElementalMatrix<Real>& d );
template<typename Real>
void AugmentedKKTRHS
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb, 
  const DistMultiVec<Real>& rmu,
        DistMultiVec<Real>& d );

template<typename Real>
void ExpandAugmentedSolution
( const Matrix<Real>& x,
  const Matrix<Real>& z,
  const Matrix<Real>& rmu,
  const Matrix<Real>& d,
        Matrix<Real>& dx,
        Matrix<Real>& dy, 
        Matrix<Real>& dz );
template<typename Real>
void ExpandAugmentedSolution
( const ElementalMatrix<Real>& x,
  const ElementalMatrix<Real>& z,
  const ElementalMatrix<Real>& rmu,
  const ElementalMatrix<Real>& d,
        ElementalMatrix<Real>& dx,
        ElementalMatrix<Real>& dy, 
        ElementalMatrix<Real>& dz );
template<typename Real>
void ExpandAugmentedSolution
( const DistMultiVec<Real>& x,
  const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& d,
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy, 
        DistMultiVec<Real>& dz );

} // namespace direct
} // namespace qp
} // namespace El
