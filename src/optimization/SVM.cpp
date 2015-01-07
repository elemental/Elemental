/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

// TODO: Document the soft-margin SVM of Cortes and Vapnik.

namespace El {

template<typename Real>
void SVM
( const Matrix<Real>& G, const Matrix<Real>& q, 
        Real gamma,            Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    LogicError("This routine is not yet written");
}

template<typename Real>
void SVM
( const AbstractDistMatrix<Real>& G, const AbstractDistMatrix<Real>& q, 
        Real gamma,                        AbstractDistMatrix<Real>& w, 
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("svm::ADMM"))
    LogicError("This routine is not yet written");
}

template<typename Real>
void SVM
( const SparseMatrix<Real>& G, const Matrix<Real>& q, 
        Real gamma,                  Matrix<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    LogicError("This routine is not yet written");
}

template<typename Real>
void SVM
( const DistSparseMatrix<Real>& G, const DistMultiVec<Real>& q, 
        Real gamma,                      DistMultiVec<Real>& x,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("SVM"))
    LogicError("This routine is not yet written");
}

#define PROTO(Real) \
  template void SVM \
  ( const Matrix<Real>& G, const Matrix<Real>& q, \
          Real gamma,            Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void SVM \
  ( const AbstractDistMatrix<Real>& G, const AbstractDistMatrix<Real>& q, \
          Real gamma,                        AbstractDistMatrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void SVM \
  ( const SparseMatrix<Real>& G, const Matrix<Real>& q, \
          Real gamma,                  Matrix<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void SVM \
  ( const DistSparseMatrix<Real>& G, const DistMultiVec<Real>& q, \
          Real gamma,                      DistMultiVec<Real>& x, \
    const qp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namepace elem
