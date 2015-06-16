/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "./BP/ADMM.hpp"
#include "./BP/IPM.hpp"

namespace El {

template<typename Real>
void BP
( const Matrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("BP"))
    if( ctrl.useIPM )
    {
        if( ctrl.useSOCP )
            bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
        else
            bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
    }
    else
        bp::ADMM( A, b, x, ctrl.admmCtrl );
}

template<typename Real>
void BP
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, 
        AbstractDistMatrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("BP"))
    if( ctrl.useIPM )
    {
        if( ctrl.useSOCP )
            bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
        else
            bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
    }
    else
        bp::ADMM( A, b, x, ctrl.admmCtrl );
}

template<typename Real>
void BP
( const SparseMatrix<Real>& A, const Matrix<Real>& b, 
        Matrix<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("BP"))
    if( !ctrl.useIPM )
        LogicError("ADMM-based BP not yet supported for sparse matrices");
    if( ctrl.useSOCP )
        bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
    else
        bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
}

template<typename Real>
void BP
( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, 
        DistMultiVec<Real>& x,
  const BPCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("BP"))
    if( !ctrl.useIPM )
        LogicError("ADMM-based BP not yet supported for sparse matrices");
    if( ctrl.useSOCP )
        bp::SOCPIPM( A, b, x, ctrl.socpIPMCtrl );
    else
        bp::LPIPM( A, b, x, ctrl.lpIPMCtrl );
}

#define PROTO(Real) \
  template void BP \
  ( const Matrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const BPCtrl<Real>& ctrl ); \
  template void BP \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
          AbstractDistMatrix<Real>& x, \
    const BPCtrl<Real>& ctrl ); \
  template void BP \
  ( const SparseMatrix<Real>& A, const Matrix<Real>& b, \
          Matrix<Real>& x, \
    const BPCtrl<Real>& ctrl ); \
  template void BP \
  ( const DistSparseMatrix<Real>& A, const DistMultiVec<Real>& b, \
          DistMultiVec<Real>& x, \
    const BPCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
