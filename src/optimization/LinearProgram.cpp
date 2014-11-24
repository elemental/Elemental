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
void LinearProgram
( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, 
  Matrix<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("LinearProgram"))
    // TODO: Switch to an Interior Point Method?
    lin_prog::ADMM( A, b, c, x );
}

template<typename Real>
void LinearProgram
( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,       AbstractDistMatrix<Real>& x )
{
    DEBUG_ONLY(CallStackEntry cse("LinearProgram"))
    // TODO: Switch to an Interior Point Method?
    lin_prog::ADMM( A, b, c, x );
}

#define PROTO(Real) \
  template void LinearProgram \
  ( const Matrix<Real>& A, const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& x ); \
  template void LinearProgram \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c,       AbstractDistMatrix<Real>& x );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
