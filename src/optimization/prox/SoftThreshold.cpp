/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename F>
F SoftThreshold( F alpha, Base<F> tau )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( tau < 0 )
            LogicError("Negative threshold does not make sense");
    )
    const Base<F> scale = Abs(alpha);
    return ( scale <= tau ? F(0) : alpha-(alpha/scale)*tau );
}

template<typename F>
void SoftThreshold( Matrix<F>& A, Base<F> tau, bool relative )
{
    DEBUG_CSE
    if( relative )
        tau *= MaxNorm(A);
    auto softThresh = [&]( F alpha ) { return SoftThreshold(alpha,tau); };
    EntrywiseMap( A, function<F(F)>(softThresh) );
}

template<typename F>
void SoftThreshold( AbstractDistMatrix<F>& A, Base<F> tau, bool relative )
{
    DEBUG_CSE
    if( relative )
        tau *= MaxNorm(A);
    auto softThresh = [&]( F alpha ) { return SoftThreshold(alpha,tau); };
    EntrywiseMap( A, function<F(F)>(softThresh) );
}

#define PROTO(F) \
  template F SoftThreshold( F alpha, Base<F> tau ); \
  template void SoftThreshold \
  ( Matrix<F>& A, Base<F> tau, bool relative ); \
  template void SoftThreshold \
  ( AbstractDistMatrix<F>& A, Base<F> tau, bool relative );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
