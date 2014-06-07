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
void HermitianInverse( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInverse"))
    SymmetricInverse( uplo, A, true, pivotType );
}

template<typename F>
void HermitianInverse
( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianInverse"))
    SymmetricInverse( uplo, A, true, pivotType );
}

template<typename F>
void LocalHermitianInverse
( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, LDLPivotType pivotType )
{
    DEBUG_ONLY(CallStackEntry cse("LocalHermitianInverse"))
    SymmetricInverse( uplo, A.Matrix(), true, pivotType );
}

#define PROTO(F) \
  template void HermitianInverse \
  ( UpperOrLower uplo, Matrix<F>& A, LDLPivotType pivotType ); \
  template void HermitianInverse \
  ( UpperOrLower uplo, DistMatrix<F>& A, LDLPivotType pivotType ); \
  template void LocalHermitianInverse \
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A, LDLPivotType pivotType );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
