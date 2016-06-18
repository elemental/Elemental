/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_MAKEREAL_HPP
#define EL_BLAS_MAKEREAL_HPP

namespace El {

template<typename Real>
void MakeReal( Matrix<Real>& A )
{ }

template<typename Real>
void MakeReal( Matrix<Complex<Real>>& A )
{
    DEBUG_CSE
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
    DEBUG_CSE
    MakeReal( A.Matrix() );
}

#ifdef EL_INSTANTIATE_BLAS_LEVEL1
# define EL_EXTERN
#else
# define EL_EXTERN extern
#endif

#define PROTO(T) \
  EL_EXTERN template void MakeReal( Matrix<T>& A ); \
  EL_EXTERN template void MakeReal( AbstractDistMatrix<T>& A ); 

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

#undef EL_EXTERN

} // namespace El

#endif // ifndef EL_BLAS_MAKEREAL_HPP
