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
void HermitianSolve
( UpperOrLower uplo, Orientation orientation, Matrix<F>& A, Matrix<F>& B, 
  LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSolve"))
    SymmetricSolve( uplo, orientation, A, B, true, pivotType );
}

template<typename F>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation, 
  DistMatrix<F>& A, DistMatrix<F>& B, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianSolve"))
    SymmetricSolve( uplo, orientation, A, B, true, pivotType );
}

#define PROTO(F) \
  template void HermitianSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<F>& A, Matrix<F>& B, LDLPivotType pivotType ); \
  template void HermitianSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    DistMatrix<F>& A, DistMatrix<F>& B, LDLPivotType pivotType );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
