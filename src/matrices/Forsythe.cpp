/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Forsythe( Matrix<T>& J, Int n, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Forsythe"))
    Jordan( J, n, lambda );
    if( n > 0 )
        J.Set( n-1, 0, alpha );
}

template<typename T>
void Forsythe( AbstractDistMatrix<T>& J, Int n, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Forsythe"))
    Jordan( J, n, lambda );
    if( n > 0 )
        J.Set( n-1, 0, alpha );
}

template<typename T>
void Forsythe( AbstractBlockDistMatrix<T>& J, Int n, T alpha, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Forsythe"))
    Jordan( J, n, lambda );
    if( n > 0 )
        J.Set( n-1, 0, alpha );
}

#define PROTO(T) \
  template void Forsythe( Matrix<T>& J, Int n, T alpha, T lambda ); \
  template void Forsythe \
  ( AbstractDistMatrix<T>& J, Int n, T alpha, T lambda ); \
  template void Forsythe \
  ( AbstractBlockDistMatrix<T>& J, Int n, T alpha, T lambda );

PROTO(Int)
PROTO(float) 
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
