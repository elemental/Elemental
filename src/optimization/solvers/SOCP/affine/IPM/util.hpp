/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "../../../LP/affine/IPM/util.hpp"

namespace El {
namespace socp {
namespace affine {

// Initialize
// ==========
template<typename Real>
void Initialize
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const DistMatrix<Real>& A,
  const DistMatrix<Real>& G,
  const DistMatrix<Real>& b,
  const DistMatrix<Real>& c,
  const DistMatrix<Real>& h,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        DistMatrix<Real>& x,
        DistMatrix<Real>& y,
        DistMatrix<Real>& z,
        DistMatrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift, Int cutoff );
template<typename Real>
void Initialize
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl );
template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  bool primalInit, bool dualInit, bool standardShift, Int cutoffPar,
  const RegSolveCtrl<Real>& solveCtrl );

// Full system
// ===========
template<typename Real>
void KKT
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& J,
  bool onlyLower=true );
template<typename Real>
void KKT
( const DistMatrix<Real>& A,
  const DistMatrix<Real>& G,
  const DistMatrix<Real>& w,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        DistMatrix<Real>& J,
  bool onlyLower=true, Int cutoff=1000 );
template<typename Real>
void KKT
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& origToSparseOrders,
  const Matrix<Int>& origToSparseFirstInds,
        Int kSparse,
        SparseMatrix<Real>& J,
  bool onlyLower=false );
template<typename Real>
void StaticKKT
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& origToSparseOrders,
  const Matrix<Int>& origToSparseFirstInds,
        Int kSparse,
        SparseMatrix<Real>& J,
  bool onlyLower=false );
template<typename Real>
void FinishKKT
( Int m, Int n,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& origToSparseOrders,
  const Matrix<Int>& origToSparseFirstInds,
        Int kSparse,
        SparseMatrix<Real>& J,
  bool onlyLower=false );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& origToSparseOrders,
  const DistMultiVec<Int>& origToSparseFirstInds,
        Int kSparse,
        DistSparseMatrix<Real>& J,
  bool onlyLower=false, Int cutoffPar=1000 );
template<typename Real>
void StaticKKT
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& origToSparseOrders,
  const DistMultiVec<Int>& origToSparseFirstInds,
        Int kSparse,
        DistSparseMatrix<Real>& J,
  bool onlyLower=false );
template<typename Real>
void FinishKKT
( Int m, Int n,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& origToSparseOrders,
  const DistMultiVec<Int>& origToSparseFirstInds,
        Int kSparse,
        DistSparseMatrix<Real>& J,
  bool onlyLower=false, Int cutoffPar=1000 );

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,
  const Matrix<Real>& rb,
  const Matrix<Real>& rh,
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& d );
template<typename Real>
void KKTRHS
( const DistMatrix<Real>& rc,
  const DistMatrix<Real>& rb,
  const DistMatrix<Real>& rh,
  const DistMatrix<Real>& rmu,
  const DistMatrix<Real>& wRoot,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& d,
  Int cutoff=1000 );
template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,
  const Matrix<Real>& rb,
  const Matrix<Real>& rh,
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& origToSparseFirstInds,
        Int kSparse,
        Matrix<Real>& d );
template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb,
  const DistMultiVec<Real>& rh,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& wRoot,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& origToSparseFirstInds,
        Int kSparse,
        DistMultiVec<Real>& d,
  Int cutoffPar=1000 );

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const Matrix<Real>& d,
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& dx,
        Matrix<Real>& dy,
        Matrix<Real>& dz,
        Matrix<Real>& ds );
template<typename Real>
void ExpandSolution
( Int m, Int n,
  const ElementalMatrix<Real>& d,
  const ElementalMatrix<Real>& rmu,
  const ElementalMatrix<Real>& wRoot,
  const ElementalMatrix<Int>& orders,
  const ElementalMatrix<Int>& firstInds,
        ElementalMatrix<Real>& dx,
        ElementalMatrix<Real>& dy,
        ElementalMatrix<Real>& dz,
        ElementalMatrix<Real>& ds,
  Int cutoff=1000 );
template<typename Real>
void ExpandSolution
( Int m, Int n,
  const Matrix<Real>& d,
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& sparseOrders,
  const Matrix<Int>& sparseFirstInds,
  const Matrix<Int>& sparseToOrigOrders,
  const Matrix<Int>& sparseToOrigFirstInds,
        Matrix<Real>& dx,
        Matrix<Real>& dy,
        Matrix<Real>& dz,
        Matrix<Real>& ds );
template<typename Real>
void ExpandSolution
( Int m, Int n,
  const DistMultiVec<Real>& d,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& wRoot,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& sparseOrders,
  const DistMultiVec<Int>& sparseFirstInds,
  const DistMultiVec<Int>& sparseToOrigOrders,
  const DistMultiVec<Int>& sparseToOrigFirstInds,
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz,
        DistMultiVec<Real>& ds,
  Int cutoffPar=1000 );

} // namespace affine
} // namespace socp
} // namespace El
