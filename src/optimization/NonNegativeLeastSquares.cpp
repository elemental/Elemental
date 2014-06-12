/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_ZEROS_INC

namespace El {

template<typename Real>
Int NonNegativeLeastSquares
( const Matrix<Real>& A, const Matrix<Real>& y, Matrix<Real>& z, 
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("NonNegativeLeastSquares"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    const Real maxReal = std::numeric_limits<Real>::max();

    Matrix<Real> P, q;
    Herk( LOWER, ADJOINT, Real(1), A, P );
    Gemv( ADJOINT, Real(-1), A, y, q );
    
    Matrix<Real> x, u;
    return QuadraticProgram
    ( P, q, Real(0), maxReal, x, z, u, rho, alpha, maxIter, absTol, relTol, inv,
     progress );
}

template<typename Real>
Int NonNegativeLeastSquares
( const DistMatrix<Real>& A, const DistMatrix<Real>& y, DistMatrix<Real>& z,
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("NonNegativeLeastSquares"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    const Real maxReal = std::numeric_limits<Real>::max();

    DistMatrix<Real> P(A.Grid()), q(A.Grid());
    Herk( LOWER, ADJOINT, Real(1), A, P );
    Gemv( ADJOINT, Real(-1), A, y, q );
    
    DistMatrix<Real> x(A.Grid()), u(A.Grid());
    return QuadraticProgram
    ( P, q, Real(0), maxReal, x, z, u, rho, alpha, maxIter, absTol, relTol, inv,
     progress );
}

#define PROTO(Real) \
  template Int NonNegativeLeastSquares \
  ( const Matrix<Real>& A, const Matrix<Real>& y, Matrix<Real>& z, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress ); \
  template Int NonNegativeLeastSquares \
  ( const DistMatrix<Real>& P, const DistMatrix<Real>& q, DistMatrix<Real>& z, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress );

PROTO(float)
PROTO(double)

} // namespace El
