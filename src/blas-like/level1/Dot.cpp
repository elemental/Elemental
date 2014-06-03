/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F> 
F Dot( const Matrix<F>& A, const Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dot"))
    return HilbertSchmidt( A, B );
}

template<typename F>
F Dot( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dot"))
    return HilbertSchmidt( A, B );
}

#define PROTO(F) \
  template F Dot( const Matrix<F>& A, const Matrix<F>& B ); \
  template F Dot \
  ( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B ); 

PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
