/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "./SVM/IPM.hpp"

namespace El {

template<typename Real>
void SVM
( const Matrix<Real>& A,
  const Matrix<Real>& d,
        Real lambda,
        Matrix<Real>& x,
  const SVMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    svm::IPM( A, d, lambda, x, ctrl.ipmCtrl );
}

template<typename Real>
void SVM
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& d,
        Real lambda,
        AbstractDistMatrix<Real>& x,
  const SVMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    svm::IPM( A, d, lambda, x, ctrl.ipmCtrl );
}

template<typename Real>
void SVM
( const SparseMatrix<Real>& A,
  const Matrix<Real>& d,
        Real lambda,
        Matrix<Real>& x,
  const SVMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    svm::IPM( A, d, lambda, x, ctrl.ipmCtrl );
}

template<typename Real>
void SVM
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& d,
        Real lambda,
        DistMultiVec<Real>& x,
  const SVMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    svm::IPM( A, d, lambda, x, ctrl.ipmCtrl );
}

#define PROTO(Real) \
  template void SVM \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& d, \
          Real lambda, \
          Matrix<Real>& x, \
    const SVMCtrl<Real>& ctrl ); \
  template void SVM \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& d, \
          Real lambda, \
          AbstractDistMatrix<Real>& x, \
    const SVMCtrl<Real>& ctrl ); \
  template void SVM \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& d, \
          Real lambda, \
          Matrix<Real>& x, \
    const SVMCtrl<Real>& ctrl ); \
  template void SVM \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& d, \
          Real lambda, \
          DistMultiVec<Real>& x, \
    const SVMCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
