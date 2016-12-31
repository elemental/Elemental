/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_SOLVERS_SOCP_HPP
#define EL_OPTIMIZATION_SOLVERS_SOCP_HPP

#include <El/optimization/solvers/util.hpp>

namespace El {

namespace SOCPApproachNS {
enum SOCPApproach {
  SOCP_ADMM,     // NOTE: Not yet supported
  SOCP_MEHROTRA
};
} // namespace SOCPApproachNS
using namespace SOCPApproachNS;

namespace socp {
namespace direct {

// Attempt to solve a pair of Second-Order Cone Programs in "direct" conic form:
//
//   min c^T x,
//   s.t. A x = b, x in K,
//
//   max -b^T y
//   s.t. A^T y - z + c = 0, z in K,
//
// where the cone K is a product of second-order cones.
//

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    SOCPApproach approach=SOCP_MEHROTRA;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl()
    {
        mehrotraCtrl.system = AUGMENTED_KKT;
        mehrotraCtrl.minTol = Pow(limits::Epsilon<Real>(),Real(0.25));
        mehrotraCtrl.targetTol = Pow(limits::Epsilon<Real>(),Real(0.5));
    }
};

} // namespace direct

namespace affine {

// Attempt to solve a pair of Second-Order Cone Programs in "affine" conic form:
//
//   min c^T x,
//   s.t. A x = b, G x + s = h, s in K,
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z in K,
//
// where the cone K is a product of second-order cones.
//

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    SOCPApproach approach=SOCP_MEHROTRA;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl()
    {
        mehrotraCtrl.minTol = Pow(limits::Epsilon<Real>(),Real(0.25));
        mehrotraCtrl.targetTol = Pow(limits::Epsilon<Real>(),Real(0.5));
    }
};

} // namespace affine
} // namespace socp

// Direct conic form
// -----------------
template<typename Real>
void SOCP
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );

// Affine conic form
// -----------------
template<typename Real>
void SOCP
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_SOLVERS_SOCP_HPP
