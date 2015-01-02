/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "./LP/direct/IPM.hpp"
#include "./LP/affine/IPM.hpp"

namespace El {

template<typename Real>
void LP
( const Matrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c, 
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_ADMM )
        lp::direct::ADMM( A, b, c, x, ctrl.admmCtrl );
    else if( ctrl.approach == LP_IPF )
        lp::direct::IPF( A, b, c, x, y, z, ctrl.ipfCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::direct::Mehrotra( A, b, c, x, y, z, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void LP
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c, 
  const Matrix<Real>& h,
        Matrix<Real>& x,       Matrix<Real>& y,
        Matrix<Real>& z,       Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_IPF )
        lp::affine::IPF( A, G, b, c, h, x, y, z, s, ctrl.ipfCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::affine::Mehrotra( A, G, b, c, h, x, y, z, s, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void LP
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z, 
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_ADMM )
        lp::direct::ADMM( A, b, c, x, ctrl.admmCtrl );
    else if( ctrl.approach == LP_IPF )
        lp::direct::IPF( A, b, c, x, y, z, ctrl.ipfCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::direct::Mehrotra( A, b, c, x, y, z, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void LP
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_IPF )
        lp::affine::IPF( A, G, b, c, h, x, y, z, s, ctrl.ipfCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::affine::Mehrotra( A, G, b, c, h, x, y, z, s, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void LP
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,       const Matrix<Real>& c, 
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_IPF )
        lp::direct::IPF( A, b, c, x, y, z, ctrl.ipfCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::direct::Mehrotra( A, b, c, x, y, z, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void LP
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c, 
  const Matrix<Real>& h,
        Matrix<Real>& x,             Matrix<Real>& y,
        Matrix<Real>& z,             Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_IPF )
        lp::affine::IPF( A, G, b, c, h, x, y, z, s, ctrl.ipfCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::affine::Mehrotra( A, G, b, c, h, x, y, z, s, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void LP
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, 
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y,
        DistMultiVec<Real>& z, 
  const lp::direct::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_IPF )
        lp::direct::IPF( A, b, c, x, y, z, ctrl.ipfCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::direct::Mehrotra( A, b, c, x, y, z, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

template<typename Real>
void LP
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, 
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,           DistMultiVec<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("LP"))
    if( ctrl.approach == LP_IPF )
        lp::affine::IPF( A, G, b, c, h, x, y, z, s, ctrl.ipfCtrl );
    else if( ctrl.approach == LP_MEHROTRA )
        lp::affine::Mehrotra( A, G, b, c, h, x, y, z, s, ctrl.mehrotraCtrl );
    else
        LogicError("Unsupported solver");
}

#define PROTO(Real) \
  template void LP \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x,       Matrix<Real>& y, \
          Matrix<Real>& z, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,       Matrix<Real>& y, \
          Matrix<Real>& z,       Matrix<Real>& s, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
          AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
          Matrix<Real>& x,             Matrix<Real>& y, \
          Matrix<Real>& z, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,             Matrix<Real>& y, \
          Matrix<Real>& z,             Matrix<Real>& s, \
    const lp::affine::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x,           DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const lp::direct::Ctrl<Real>& ctrl ); \
  template void LP \
  ( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x,           DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z,           DistMultiVec<Real>& s, \
    const lp::affine::Ctrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
