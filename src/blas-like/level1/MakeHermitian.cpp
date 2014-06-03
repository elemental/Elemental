/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename T>
void MakeHermitian( UpperOrLower uplo, Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHermitian"))
    MakeSymmetric( uplo, A, true );
}

template<typename T>
void MakeHermitian( UpperOrLower uplo, AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeHermitian"))
    MakeSymmetric( uplo, A, true );
}

#define PROTO(F) \
  template void MakeHermitian( UpperOrLower uplo, Matrix<F>& A ); \
  template void MakeHermitian( UpperOrLower uplo, AbstractDistMatrix<F>& A );

PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
