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
void Zero( Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Zero"))
    const Int height = A.Height();
    const Int width = A.Width();
    EL_PARALLEL_FOR
    for( Int j=0; j<width; ++j )
        MemZero( A.Buffer(0,j), height );
}

template<typename T>
void Zero( AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Zero"))
    Zero( A.Matrix() );
}

#define PROTO(T) \
  template void Zero( Matrix<T>& A ); \
  template void Zero( AbstractDistMatrix<T>& A );

PROTO(Int);
PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
