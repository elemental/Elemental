/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace nnls {

// Transform each problem
//
//   min || A x - b ||_2
//   s.t. x >= 0
//
// into the explicit QP
//
//   min (1/2) x^T (A^T A) x + (-A^T b)^T x
//   s.t. x >= 0
//
// and solve the sequence of problems simultaneously with ADMM.
//

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Int
ADMM
( const Matrix<Real>& A,
  const Matrix<Real>& B,
        Matrix<Real>& X,
  const ADMMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real maxReal = limits::Max<Real>();

    Matrix<Real> Q, C;
    Herk( LOWER, ADJOINT, Real(1), A, Q );
    Gemm( ADJOINT, NORMAL, Real(-1), A, B, C );

    return qp::box::ADMM( Q, C, Real(0), maxReal, X, ctrl );
}

template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Int
ADMM
( const AbstractDistMatrix<Real>& APre,
  const AbstractDistMatrix<Real>& B,
        AbstractDistMatrix<Real>& X,
  const ADMMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real maxReal = limits::Max<Real>();

    DistMatrixReadProxy<Real,Real,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    DistMatrix<Real> Q(A.Grid()), C(A.Grid());
    Herk( LOWER, ADJOINT, Real(1), A, Q );
    Gemm( ADJOINT, NORMAL, Real(-1), A, B, C );

    return qp::box::ADMM( Q, C, Real(0), maxReal, X, ctrl );
}

} // namespace nnls
} // namespace El
