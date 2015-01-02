/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./General/LUPartialPiv.hpp"

namespace El {

template<typename F> 
void Inverse( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Inverse"))
    inverse::LUPartialPiv( A );
}

template<typename F> 
void Inverse( AbstractDistMatrix<F>& A )
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
  template void Inverse( AbstractDistMatrix<F>& A ); \
  template void LocalInverse( DistMatrix<F,STAR,STAR>& A ); \
  template void inverse::AfterLUPartialPiv \
  ( Matrix<F>& A, const Matrix<Int>& p ); \
  template void inverse::AfterLUPartialPiv \
  ( AbstractDistMatrix<F>& A, const AbstractDistMatrix<Int>& p ); 

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
