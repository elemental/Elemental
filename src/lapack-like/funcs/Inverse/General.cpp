/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

#include "./General/LUPartialPiv.hpp"

namespace El {

template<typename F> 
void Inverse( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Inverse"))
    inverse::LUPartialPiv( A );
}

template<typename F> 
void Inverse( DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Inverse"))
    inverse::LUPartialPiv( A );
}

template<typename F>
void LocalInverse( DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LocalInverse"))
    Inverse( A.Matrix() );
}

#define PROTO(F) \
  template void Inverse( Matrix<F>& A ); \
  template void Inverse( DistMatrix<F>& A ); \
  template void LocalInverse( DistMatrix<F,STAR,STAR>& A ); \
  template void inverse::AfterLUPartialPiv \
  ( Matrix<F>& A, const Matrix<Int>& pPerm ); \
  template void inverse::AfterLUPartialPiv \
  ( DistMatrix<F>& A, const DistMatrix<Int,VC,STAR>& pPerm ); 

PROTO(float) 
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
