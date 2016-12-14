/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El-lite.hpp>
#include <El/matrices.hpp>

namespace El {

template<typename T>
void Forsythe( Matrix<T>& J, Int n, T alpha, T lambda )
{
    EL_DEBUG_CSE
    Jordan( J, n, lambda );
    if( n > 0 )
        J.Set( n-1, 0, alpha );
}

template<typename T>
void Forsythe( AbstractDistMatrix<T>& J, Int n, T alpha, T lambda )
{
    EL_DEBUG_CSE
    Jordan( J, n, lambda );
    if( n > 0 )
        J.Set( n-1, 0, alpha );
}

#define PROTO(T) \
  template void Forsythe( Matrix<T>& J, Int n, T alpha, T lambda ); \
  template void Forsythe \
  ( AbstractDistMatrix<T>& J, Int n, T alpha, T lambda );

#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGINT
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
