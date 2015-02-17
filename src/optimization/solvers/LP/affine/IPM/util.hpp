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
namespace lp {
namespace affine {

// Initialize
// ==========
template<typename Real>
void Initialize
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,       Matrix<Real>& s,
  bool primalInitialized, bool dualInitialized,
  bool standardShift );
template<typename Real>
void Initialize
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s,
  bool primalInitialized, bool dualInitialized,
  bool standardShift );
template<typename Real>
void Initialize
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z,             Matrix<Real>& s,
        vector<Int>& map,            vector<Int>& invMap,
        Separator& rootSep,          SymmNodeInfo& info,
  bool primalInitialized, bool dualInitialized,
  bool standardShift,     bool progress );
template<typename Real>
void Initialize
( const DistSparseMatrix<Real>& A,  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,      const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,            DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,            DistMultiVec<Real>& s,
        DistMap& map,                     DistMap& invMap,
        DistSeparator& rootSep,           DistSymmNodeInfo& info,
  bool primalInitialized, bool dualInitialized,
  bool standardShift,     bool progress );

// Full system
// ===========
template<typename Real>
void KKT
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& s, const Matrix<Real>& z,
        Matrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& s, const AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& s,       const Matrix<Real>& z,
        SparseMatrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& s,     const DistMultiVec<Real>& z,
        DistSparseMatrix<Real>& J, bool onlyLower=true );

using qp::affine::KKTRHS;
using qp::affine::ExpandCoreSolution;
using qp::affine::ExpandSolution;

// Line search
// ===========
template<typename Real>
Real IPFLineSearch
( const Matrix<Real>& A,  const Matrix<Real>& G,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Real>& x,  const Matrix<Real>& y,  
  const Matrix<Real>& z,  const Matrix<Real>& s,
  const Matrix<Real>& dx, const Matrix<Real>& dy, 
  const Matrix<Real>& dz, const Matrix<Real>& ds,
  Real upperBound,
  Real bTol, Real cTol, Real hTol,
  const lp::IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const AbstractDistMatrix<Real>& A,  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
  const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Real>& z,  const AbstractDistMatrix<Real>& s,
  const AbstractDistMatrix<Real>& dx, const AbstractDistMatrix<Real>& dy,
  const AbstractDistMatrix<Real>& dz, const AbstractDistMatrix<Real>& ds,
  Real upperBound,
  Real bTol, Real cTol, Real hTol,
  const lp::IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Real>& x,       const Matrix<Real>& y,  
  const Matrix<Real>& z,       const Matrix<Real>& s,
  const Matrix<Real>& dx,      const Matrix<Real>& dy, 
  const Matrix<Real>& dz,      const Matrix<Real>& ds,
  Real upperBound,
  Real bTol, Real cTol, Real hTol,
  const lp::IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& y,
  const DistMultiVec<Real>& z,     const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& dx,    const DistMultiVec<Real>& dy,
  const DistMultiVec<Real>& dz,    const DistMultiVec<Real>& ds,
  Real upperBound,
  Real bTol, Real cTol, Real hTol,
  const lp::IPFLineSearchCtrl<Real>& ctrl );

} // namespace affine
} // namespace lp
} // namespace El
