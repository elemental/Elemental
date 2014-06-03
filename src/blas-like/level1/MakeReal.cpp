/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename Real>
void MakeReal( Matrix<Real>& A )
{ }

template<typename Real>
void MakeReal( Matrix<Complex<Real>>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeReal"))
    Complex<Real>* ABuffer = A.Buffer();
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

#define PROTO(T) \
  template void MakeReal( Matrix<T>& A ); \
  template void MakeReal( AbstractDistMatrix<T>& A ); 

PROTO(Int);
PROTO(float);
PROTO(double);
PROTO(Complex<float>);
PROTO(Complex<double>);

} // namespace El
