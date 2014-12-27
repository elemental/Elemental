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
namespace dual {

// Full system
// -----------
template<typename Real>
void KKT
( const Matrix<Real>& A, const Matrix<Real>& G, 
  const Matrix<Real>& w, 
  Matrix<Real>& J, bool onlyLower=true );
template<typename Real>
void KKT
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& w, 
  AbstractDistMatrix<Real>& J, bool onlyLower=true );

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc, const Matrix<Real>& rb, 
  const Matrix<Real>& rh, const Matrix<Real>& rmu,
  const Matrix<Real>& w,  const Matrix<Real>& l,
        Matrix<Real>& d );
template<typename Real>
void KKTRHS
( const AbstractDistMatrix<Real>& rc, const AbstractDistMatrix<Real>& rb, 
  const AbstractDistMatrix<Real>& rh, const AbstractDistMatrix<Real>& rmu,
  const AbstractDistMatrix<Real>& w,  const AbstractDistMatrix<Real>& l,
        AbstractDistMatrix<Real>& d );

template<typename Real>
void ExpandKKTSolution
( Int m, Int n, Int k, const Matrix<Real>& d,
  const Matrix<Real>& rmu,
  const Matrix<Real>& w, const Matrix<Real>& l,
  Matrix<Real>& dx, Matrix<Real>& dy, 
  Matrix<Real>& dz, Matrix<Real>& ds );
template<typename Real>
void ExpandKKTSolution
( Int m, Int n, Int k, const AbstractDistMatrix<Real>& d,
  const AbstractDistMatrix<Real>& rmu,
  const AbstractDistMatrix<Real>& w, const AbstractDistMatrix<Real>& l,
  AbstractDistMatrix<Real>& dx, AbstractDistMatrix<Real>& dy, 
  AbstractDistMatrix<Real>& dz, AbstractDistMatrix<Real>& ds );

// Line search
// -----------
template<typename Real>
Real IPFLineSearch
( const Matrix<Real>& A,  const Matrix<Real>& G,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Real>& x,  const Matrix<Real>& y,  
  const Matrix<Real>& z,  const Matrix<Real>& s,
  const Matrix<Real>& dx, const Matrix<Real>& dy, 
  const Matrix<Real>& dz, const Matrix<Real>& ds,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const AbstractDistMatrix<Real>& A,  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
  const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Real>& z,  const AbstractDistMatrix<Real>& s,
  const AbstractDistMatrix<Real>& dx, const AbstractDistMatrix<Real>& dy,
  const AbstractDistMatrix<Real>& dz, const AbstractDistMatrix<Real>& ds,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Real>& x,  const Matrix<Real>& y,  
  const Matrix<Real>& z,  const Matrix<Real>& s,
  const Matrix<Real>& dx, const Matrix<Real>& dy, 
  const Matrix<Real>& dz, const Matrix<Real>& ds,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );
template<typename Real>
Real IPFLineSearch
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Real>& x,  const DistMultiVec<Real>& y,
  const DistMultiVec<Real>& z,  const DistMultiVec<Real>& s,
  const DistMultiVec<Real>& dx, const DistMultiVec<Real>& dy,
  const DistMultiVec<Real>& dz, const DistMultiVec<Real>& ds,
  Real bTol, Real cTol,
  const IPFLineSearchCtrl<Real>& ctrl );

} // namespace dual
} // namespace lp
} // namespace El
