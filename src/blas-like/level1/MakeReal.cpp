/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include EL_MAKEREAL_INC

namespace El {

template<typename T>
void MakeReal( Matrix<Complex<T>>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeReal"))
    Complex<T>* ABuffer = A.Buffer();
    const Int height = A.Height();
    const Int width = A.Width();
    const Int ldim = A.LDim();
    for( Int j=0; j<width; ++j )
        for( Int i=0; i<height; ++i )
            ABuffer[i+j*ldim] = RealPart(ABuffer[i+j*ldim]);
}

template<typename T>
void MakeReal( AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeReal"))
    MakeReal( A.Matrix() );
}

#define PROTO_REAL(T) \
  template void MakeReal( AbstractDistMatrix<T>& A );

#define PROTO_CPX(T) \
  template void MakeReal( Matrix<T>& A ); \
  template void MakeReal( AbstractDistMatrix<T>& A );

PROTO_REAL(float);
PROTO_REAL(double);
PROTO_CPX(Complex<float>);
PROTO_CPX(Complex<double>);

} // namespace El
