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

// TODO: Experiment with running a few iterations of IPF until the residuals
//       are sufficiently low before switching to Mehrotra's scheme

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
// using a Mehrotra Predictor-Corrector scheme.
//
template<typename Real>
void Mehrotra
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::Mehrotra"))    
    LogicError("This routine is not yet written");
}

template<typename Real>
void Mehrotra
( const AbstractDistMatrix<Real>& APre, 
  const AbstractDistMatrix<Real>& GPre,
  const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c,
  AbstractDistMatrix<Real>& xPre, AbstractDistMatrix<Real>& y,
  AbstractDistMatrix<Real>& zPre,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::Mehrotra"))    

    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    auto APtr = ReadProxy<Real,MC,MR>(&APre,control);      auto& A = *APtr;
    auto GPtr = ReadProxy<Real,MC,MR>(&GPre,control);      auto& G = *GPtr;
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,control); auto& x = *xPtr;
    auto zPtr = ReadWriteProxy<Real,MC,MR>(&zPre,control); auto& z = *zPtr;

    LogicError("This routine is not yet written");
}

template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::Mehrotra"))    
    LogicError("Sequential sparse-direct solvers not yet supported");
}

template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& A, 
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  DistMultiVec<Real>& x, DistMultiVec<Real>& y, DistMultiVec<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::Mehrotra"))    
    LogicError("This routine is not yet written");
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    AbstractDistMatrix<Real>& x, AbstractDistMatrix<Real>& y, \
    AbstractDistMatrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    DistMultiVec<Real>& x, DistMultiVec<Real>& y, DistMultiVec<Real>& z, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace dual
} // namespace lp
} // namespace El
