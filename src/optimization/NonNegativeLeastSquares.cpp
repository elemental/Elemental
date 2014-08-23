/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename Real>
Int NonNegativeLeastSquares
( const Matrix<Real>& A, const Matrix<Real>& Y, Matrix<Real>& Z, 
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("NonNegativeLeastSquares"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    const Real maxReal = std::numeric_limits<Real>::max();

    Matrix<Real> P, S;
    Herk( LOWER, ADJOINT, Real(1), A, P );
    Gemm( ADJOINT, NORMAL, Real(-1), A, Y, S );
    
    return QuadraticProgram
    ( P, S, Real(0), maxReal, Z, rho, alpha, maxIter, absTol, relTol, inv,
     progress );
}

template<typename Real>
Int NonNegativeLeastSquares
( const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& Y, 
        AbstractDistMatrix<Real>& Z,
  Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, 
  bool inv, bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("NonNegativeLeastSquares"))
    if( IsComplex<Real>::val ) 
        LogicError("The datatype was assumed to be real");
    const Real maxReal = std::numeric_limits<Real>::max();

    auto APtr = ReadProxy( &APre );
    auto& A = *APtr;

    DistMatrix<Real> P(A.Grid()), S(A.Grid());
    Herk( LOWER, ADJOINT, Real(1), A, P );
    Gemm( ADJOINT, NORMAL, Real(-1), A, Y, S );
    
    return QuadraticProgram
    ( P, S, Real(0), maxReal, Z, rho, alpha, maxIter, absTol, relTol, inv,
     progress );
}

#define PROTO(Real) \
  template Int NonNegativeLeastSquares \
  ( const Matrix<Real>& A, const Matrix<Real>& Y, Matrix<Real>& Z, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress ); \
  template Int NonNegativeLeastSquares \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& Y, \
          AbstractDistMatrix<Real>& Z, \
    Real rho, Real alpha, Int maxIter, Real absTol, Real relTol, \
    bool inv, bool progress );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
