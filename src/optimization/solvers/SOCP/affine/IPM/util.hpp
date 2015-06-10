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
        Matrix<Real>& x, 
        Matrix<Real>& y,
        Matrix<Real>& z, 
        Matrix<Real>& s,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b, 
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x, 
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
  bool primalInit, bool dualInit, bool standardShift );
template<typename Real>
void Initialize
( const SparseMatrix<Real>& A, 
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b, 
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x, 
        Matrix<Real>& y,
        Matrix<Real>& z, 
        Matrix<Real>& s,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
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
        DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z, 
        DistMultiVec<Real>& s,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMap& map,
        DistMap& invMap,
        ldl::DistSeparator& rootSep, 
        ldl::DistNodeInfo& info,
  bool primalInit, bool dualInit, bool standardShift,
  const RegQSDCtrl<Real>& qsdCtrl );

// Full system
// ===========
template<typename Real>
void KKT
( const Matrix<Real>& A, 
  const Matrix<Real>& G,
  const Matrix<Real>& s, 
  const Matrix<Real>& z,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& s, 
  const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const SparseMatrix<Real>& A, 
  const SparseMatrix<Real>& G,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        SparseMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A, 
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& z,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistSparseMatrix<Real>& J, bool onlyLower=true );

using qp::affine::KKTRHS;
using qp::affine::ExpandCoreSolution;
using qp::affine::ExpandSolution;

} // namespace affine
} // namespace socp
} // namespace El
