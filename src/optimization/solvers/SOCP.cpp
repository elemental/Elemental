/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
//#include "./SOCP/direct/IPM.hpp"
#include "./SOCP/affine/IPM.hpp"

namespace El {

template<typename Real>
void SOCP
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
  const Matrix<Int>& labels,
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("SOCP"))
    if( ctrl.approach == SOCP_MEHROTRA )
        socp::affine::Mehrotra
        ( A, G, b, c, h, x, y, z, s, orders, firstInds, labels,
          ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void SOCP
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
  const AbstractDistMatrix<Int>& labels,
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("SOCP"))
    if( ctrl.approach == SOCP_MEHROTRA )
        socp::affine::Mehrotra
        ( A, G, b, c, h, x, y, z, s, orders, firstInds, labels, 
          ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void SOCP
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
  const Matrix<Int>& labels,
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("SOCP"))
    if( ctrl.approach == SOCP_MEHROTRA )
        socp::affine::Mehrotra
        ( A, G, b, c, h, x, y, z, s, orders, firstInds, labels,
          ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void SOCP
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
  const DistMultiVec<Int>& labels,
  const socp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("SOCP"))
    if( ctrl.approach == SOCP_MEHROTRA )
        socp::affine::Mehrotra
        ( A, G, b, c, h, x, y, z, s, orders, firstInds, labels, 
          ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

#define PROTO(Real) \
  template void SOCP \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& labels, \
    const socp::affine::Ctrl<Real>& ctrl ); \
  template void SOCP \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
          AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& s, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    const AbstractDistMatrix<Int>& labels, \
    const socp::affine::Ctrl<Real>& ctrl ); \
  template void SOCP \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& labels, \
    const socp::affine::Ctrl<Real>& ctrl ); \
  template void SOCP \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMultiVec<Real>& s, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    const DistMultiVec<Int>& labels, \
    const socp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
