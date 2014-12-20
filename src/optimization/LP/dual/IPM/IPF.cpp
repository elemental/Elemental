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

// The following solves a linear program in "dual" conic form:
//
//   min c^T x
//   s.t. A x = b, G x <= h,
//
// as opposed to the more specific "primal" conic form:
//
//   min c^T x
//   s.t. A x = b, x >= 0,
//
// using a simple Infeasible Path Following (IPF) scheme. This routine
// should only be used for academic purposes, as the Mehrotra alternative
// typically requires an order of magnitude fewer iterations.
template<typename Real>
void IPF
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::IPF"))    
    LogicError("This routine is not yet written");
}

template<typename Real>
void IPF
( const AbstractDistMatrix<Real>& APre, 
  const AbstractDistMatrix<Real>& GPre,
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  AbstractDistMatrix<Real>& xPre, AbstractDistMatrix<Real>& y, 
  AbstractDistMatrix<Real>& zPre, const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::IPF"))    

    ProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = 0;
    proxCtrl.rowAlign = 0;
    auto APtr = ReadProxy<Real,MC,MR>(&APre,proxCtrl);      auto& A = *APtr;
    auto GPtr = ReadProxy<Real,MC,MR>(&GPre,proxCtrl);      auto& G = *GPtr;
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,proxCtrl); auto& x = *xPtr;
    auto zPtr = ReadWriteProxy<Real,MC,MR>(&zPre,proxCtrl); auto& z = *zPtr;

    LogicError("This routine is not yet written");
}

template<typename Real>
void IPF
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::IPF"))    
    LogicError("Sequential sparse-direct solvers not yet supported");
}

template<typename Real>
void IPF
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  DistMultiVec<Real>& x, DistMultiVec<Real>& y, DistMultiVec<Real>& z,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::IPF"))    
    LogicError("This routine is not yet written");
}

#define PROTO(Real) \
  template void IPF \
  ( const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    AbstractDistMatrix<Real>& x, AbstractDistMatrix<Real>& y, \
    AbstractDistMatrix<Real>& z, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    DistMultiVec<Real>& x, DistMultiVec<Real>& y, DistMultiVec<Real>& z, \
    const IPFCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace dual
} // namespace lp
} // namespace El
