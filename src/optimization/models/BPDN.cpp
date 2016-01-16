/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "./BPDN/ADMM.hpp"
#include "./BPDN/IPM.hpp"

namespace El {

template<typename Real>
void BPDN
( const Matrix<Real>& A, 
  const Matrix<Real>& b, 
        Real lambda,
        Matrix<Real>& x,
  const BPDNCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("BPDN"))
    if( ctrl.useIPM )
        bpdn::IPM( A, b, lambda, x, ctrl.ipmCtrl );
    else
        bpdn::ADMM( A, b, lambda, x, ctrl.admmCtrl );
}

template<typename Real>
void BPDN
( const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& b, 
        Real lambda,
        ElementalMatrix<Real>& x,
  const BPDNCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("BPDN"))
    if( ctrl.useIPM )
        bpdn::IPM( A, b, lambda, x, ctrl.ipmCtrl );
    else
        bpdn::ADMM( A, b, lambda, x, ctrl.admmCtrl );
}

template<typename Real>
void BPDN
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b, 
        Real lambda, 
        Matrix<Real>& x,
  const BPDNCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("BPDN"))
    if( !ctrl.useIPM )
        LogicError("ADMM-based BPDN not yet supported for sparse matrices");
    bpdn::IPM( A, b, lambda, x, ctrl.ipmCtrl );
}

template<typename Real>
void BPDN
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b, 
        Real lambda,
        DistMultiVec<Real>& x,
  const BPDNCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("BPDN"))
    if( !ctrl.useIPM )
        LogicError("ADMM-based BPDN not yet supported for sparse matrices");
    bpdn::IPM( A, b, lambda, x, ctrl.ipmCtrl );
}

#define PROTO(Real) \
  template void BPDN \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const BPDNCtrl<Real>& ctrl ); \
  template void BPDN \
  ( const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& b, \
          Real lambda, \
          ElementalMatrix<Real>& x, \
    const BPDNCtrl<Real>& ctrl ); \
  template void BPDN \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
          Real lambda, \
          Matrix<Real>& x, \
    const BPDNCtrl<Real>& ctrl ); \
  template void BPDN \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
          Real lambda, \
          DistMultiVec<Real>& x, \
    const BPDNCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
