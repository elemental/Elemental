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
void Ones( Matrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    Fill( A, T(1) );
}

template<typename T>
void Ones( AbstractDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    Fill( A, T(1) );
}

template<typename T>
void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    Fill( A, T(1) );
}

#define PROTO(T) \
  template void Ones( Matrix<T>& A, Int m, Int n ); \
  template void Ones( AbstractDistMatrix<T>& A, Int m, Int n ); \
  template void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n );

PROTO(Int);
#ifndef EL_DISABLE_FLOAT
PROTO(float);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<float>);
#endif // ifndef EL_DISABLE_COMPLEX
#endif // ifndef EL_DISABLE_FLOAT

PROTO(double);
#ifndef EL_DISABLE_COMPLEX
PROTO(Complex<double>);
#endif // ifndef EL_DISABLE_COMPLEX

} // namespace El
