/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "../../../QP/affine/IPM/util.hpp"

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
  const Matrix<Int>& labels,
        Matrix<Real>& x, 
        Matrix<Real>& y,
        Matrix<Real>& z, 
        Matrix<Real>& s,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b, 
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  const AbstractDistMatrix<Int>& labels,
        AbstractDistMatrix<Real>& x, 
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,

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
  const Matrix<Int>& labels,
        Matrix<Real>& x, 
        Matrix<Real>& y,
        Matrix<Real>& z, 
        Matrix<Real>& s,
        vector<Int>& map, 
        vector<Int>& invMap,
        ldl::Separator& rootSep, 
        ldl::NodeInfo& info,
  bool primalInit, bool dualInit, bool standardShift, 
  const RegQSDCtrl<Real>& qsdCtrl );
template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& A, 
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& labels,
        DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z, 
        DistMultiVec<Real>& s,
        DistMap& map,
        DistMap& invMap,
        ldl::DistSeparator& rootSep, 
        ldl::DistNodeInfo& info,
  bool primalInit, bool dualInit, bool standardShift, Int cutoff,
  const RegQSDCtrl<Real>& qsdCtrl );

// Full system
// ===========
template<typename Real>
void KKT
( const Matrix<Real>& A, 
  const Matrix<Real>& G,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& labels,
        Matrix<Real>& J, 
  bool onlyLower=true );
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& w,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  const AbstractDistMatrix<Int>& labels,
        AbstractDistMatrix<Real>& J, 
  bool onlyLower=true, Int cutoff=1000 );
template<typename Real>
void KKT
( const SparseMatrix<Real>& A, 
  const SparseMatrix<Real>& G,
  const Matrix<Real>& w,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& labels,
        SparseMatrix<Real>& J, 
  bool onlyLower=true );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A, 
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& w,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& labels,
        DistSparseMatrix<Real>& J, 
  bool onlyLower=true, Int cutoff=1000 );

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc, 
  const Matrix<Real>& rb,
  const Matrix<Real>& rh, 
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& labels,
        Matrix<Real>& d );
template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rc,
  const AbstractDistMatrix<Real>& rb,
  const AbstractDistMatrix<Real>& rh,
  const AbstractDistMatrix<Real>& rmu,
  const AbstractDistMatrix<Real>& wRoot,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  const AbstractDistMatrix<Int>& labels,
        AbstractDistMatrix<Real>& d,
  Int cutoff=1000 );
template<typename Real>
void KKTRHS
( const DistMultiVec<Real>& rc,
  const DistMultiVec<Real>& rb,
  const DistMultiVec<Real>& rh,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& wRoot,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& labels,
        DistMultiVec<Real>& d,
  Int cutoff=1000 );

template<typename Real>
void ExpandSolution
( Int m, Int n,
  const Matrix<Real>& d,
  const Matrix<Real>& rmu,
  const Matrix<Real>& wRoot,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& labels,
        Matrix<Real>& dx,
        Matrix<Real>& dy,
        Matrix<Real>& dz,
        Matrix<Real>& ds );
template<typename Real>
void ExpandSolution
( Int m, Int n,
  const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& rmu,
  const AbstractDistMatrix<Real>& wRoot,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  const AbstractDistMatrix<Int>& labels,
        AbstractDistMatrix<Real>& dx,
        AbstractDistMatrix<Real>& dy,
        AbstractDistMatrix<Real>& dz,
        AbstractDistMatrix<Real>& ds,
  Int cutoff=1000 );
template<typename Real>
void ExpandSolution
( Int m, Int n,
  const DistMultiVec<Real>& d,
  const DistMultiVec<Real>& rmu,
  const DistMultiVec<Real>& wRoot,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& labels,
        DistMultiVec<Real>& dx,
        DistMultiVec<Real>& dy,
        DistMultiVec<Real>& dz,
        DistMultiVec<Real>& ds,
  Int cutoff=1000 );

} // namespace affine
} // namespace socp
} // namespace El
