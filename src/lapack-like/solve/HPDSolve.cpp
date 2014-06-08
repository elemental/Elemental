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
void HPDSolve
( UpperOrLower uplo, Orientation orientation, 
  Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("HPDSolve"))
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

template<typename F>
void HPDSolve
( UpperOrLower uplo, Orientation orientation, 
  DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("HPDSolve"))
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

#define PROTO(F) \
  template void HPDSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<F>& A, Matrix<F>& B ); \
  template void HPDSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    DistMatrix<F>& A, DistMatrix<F>& B );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
