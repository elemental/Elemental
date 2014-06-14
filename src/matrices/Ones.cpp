/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_ONES_INC

namespace El {

template<typename T> 
void MakeOnes( Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, T(1) );
}

template<typename T>
void MakeOnes( AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    MakeOnes( A.Matrix() );
}

template<typename T>
void MakeOnes( AbstractBlockDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    MakeOnes( A.Matrix() );
}

template<typename T>
void Ones( Matrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

template<typename T>
void Ones( AbstractDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

template<typename T>
void Ones( AbstractBlockDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

#define PROTO(T) \
  template void MakeOnes( Matrix<T>& A ); \
  template void MakeOnes( AbstractDistMatrix<T>& A ); \
  template void MakeOnes( AbstractBlockDistMatrix<T>& A ); \
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
