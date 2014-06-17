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
void Fill( Matrix<T>& A, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Fill"))
    const Int height = A.Height();
    const Int width = A.Width();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            A.Set( i, j, alpha );
}

template<typename T>
void Fill( AbstractDistMatrix<T>& A, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Fill"))
    Fill( A.Matrix(), alpha );
}

template<typename T>
void Fill( AbstractBlockDistMatrix<T>& A, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Fill"))
    Fill( A.Matrix(), alpha );
}

#define PROTO(T) \
  template void Fill( Matrix<T>& A, T alpha ); \
  template void Fill( AbstractDistMatrix<T>& A, T alpha ); \
  template void Fill( AbstractBlockDistMatrix<T>& A, T alpha );

PROTO(Int);
PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
