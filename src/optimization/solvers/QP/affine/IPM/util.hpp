/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace qp {
namespace affine {

// Initialize
// ==========
template<typename Real>
void Initialize
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const SparseMatrix<Real>& JStatic,
  const Matrix<Real>& regTmp,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const vector<Int>& map,
  const vector<Int>& invMap,
  const ldl::Separator& rootSep,
  const ldl::NodeInfo& info,
  bool primalInit, bool dualInit, bool standardShift, 
  const RegSolveCtrl<Real>& solveCtrl );
template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& JStatic,
  const DistMultiVec<Real>& regTmp,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const DistMap& map,
  const DistMap& invMap,
  const ldl::DistSeparator& rootSep,
  const ldl::DistNodeInfo& info,
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
  const Matrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& s,
  const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void StaticKKT
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        SparseMatrix<Real>& J,
  bool onlyLower );
template<typename Real>
void FinishKKT
( Int m, Int n,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void StaticKKT
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        DistSparseMatrix<Real>& J,
  bool onlyLower );
template<typename Real>
void FinishKKT
( Int m, Int n,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J );

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,
  const Matrix<Real>& rb,
  const Matrix<Real>& rh,
  const Matrix<Real>& rmu, 
  const Matrix<Real>& z,
        Matrix<Real>& d );
template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rc,
  const AbstractDistMatrix<Real>& rb,
  const AbstractDistMatrix<Real>& rh,
  const AbstractDistMatrix<Real>& rmu, 
  const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& d );
template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb,
  const DistMultiVec<Real>& rh,
  const DistMultiVec<Real>& rmu, 
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& d );

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const Matrix<Real>& d, 
        Matrix<Real>& dx,
        Matrix<Real>& dy, 
        Matrix<Real>& dz );
template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const AbstractDistMatrix<Real>& d,
        AbstractDistMatrix<Real>& dx,
        AbstractDistMatrix<Real>& dy,
        AbstractDistMatrix<Real>& dz );
template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const DistMultiVec<Real>& d, 
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz );

template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const Matrix<Real>& d,
  const Matrix<Real>& rmu, 
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        Matrix<Real>& dx,
        Matrix<Real>& dy, 
        Matrix<Real>& dz,
        Matrix<Real>& ds );
template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& rmu,
  const AbstractDistMatrix<Real>& s,
  const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& dx,
        AbstractDistMatrix<Real>& dy,
        AbstractDistMatrix<Real>& dz,
        AbstractDistMatrix<Real>& ds );
template<typename Real>
void ExpandSolution
( Int m, Int n, 
  const DistMultiVec<Real>& d,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& z,
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz,
        DistMultiVec<Real>& ds );

} // namespace affine
} // namespace qp
} // namespace El
