/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "./QP/direct/IPM.hpp"
#include "./QP/affine/IPM.hpp"

namespace El {

// Direct conic form
// =================
template<typename Real>
void QP
( const Matrix<Real>& Q, 
  const Matrix<Real>& A,
  const Matrix<Real>& b, 
  const Matrix<Real>& c, 
        Matrix<Real>& x, 
        Matrix<Real>& y,
        Matrix<Real>& z, 
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.approach == QP_MEHROTRA )
        qp::direct::Mehrotra( Q, A, b, c, x, y, z, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void QP
( const ElementalMatrix<Real>& Q, 
  const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& b, 
  const ElementalMatrix<Real>& c, 
        ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& y,
        ElementalMatrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.approach == QP_MEHROTRA )
        qp::direct::Mehrotra( Q, A, b, c, x, y, z, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void QP
( const SparseMatrix<Real>& Q, 
  const SparseMatrix<Real>& A,
  const Matrix<Real>& b, 
  const Matrix<Real>& c, 
        Matrix<Real>& x, 
        Matrix<Real>& y,
        Matrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.approach == QP_MEHROTRA )
        qp::direct::Mehrotra( Q, A, b, c, x, y, z, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c, 
        DistMultiVec<Real>& x, 
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.approach == QP_MEHROTRA )
        qp::direct::Mehrotra( Q, A, b, c, x, y, z, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

// Affine conic form
// =================
template<typename Real>
void QP
( const Matrix<Real>& Q, 
  const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c, 
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s, 
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.approach == QP_MEHROTRA )
        qp::affine::Mehrotra( Q, A, G, b, c, h, x, y, z, s, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void QP
( const ElementalMatrix<Real>& Q, 
  const ElementalMatrix<Real>& A,
  const ElementalMatrix<Real>& G,
  const ElementalMatrix<Real>& b,
  const ElementalMatrix<Real>& c, 
  const ElementalMatrix<Real>& h,
        ElementalMatrix<Real>& x,
        ElementalMatrix<Real>& y,
        ElementalMatrix<Real>& z,
        ElementalMatrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.approach == QP_MEHROTRA )
        qp::affine::Mehrotra( Q, A, G, b, c, h, x, y, z, s, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void QP
( const SparseMatrix<Real>& Q, 
  const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c, 
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.approach == QP_MEHROTRA )
        qp::affine::Mehrotra( Q, A, G, b, c, h, x, y, z, s, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q, 
  const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c, 
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_CSE
    if( ctrl.approach == QP_MEHROTRA )
        qp::affine::Mehrotra( Q, A, G, b, c, h, x, y, z, s, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

#define PROTO(Real) \
  template void QP \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& b, \
    const ElementalMatrix<Real>& c, \
          ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const qp::direct::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& G, \
    const ElementalMatrix<Real>& b, \
    const ElementalMatrix<Real>& c, \
    const ElementalMatrix<Real>& h, \
          ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& s, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const qp::affine::Ctrl<Real>& ctrl ); \
  template void QP \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMultiVec<Real>& s, \
    const qp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
