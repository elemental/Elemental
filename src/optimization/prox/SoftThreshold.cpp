/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
Field SoftThreshold( Field alpha, Base<Field> tau )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( tau < 0 )
            LogicError("Negative threshold does not make sense");
    )
    const Base<Field> scale = Abs(alpha);
    return ( scale <= tau ? Field(0) : alpha-(alpha/scale)*tau );
}

template<typename Field>
void SoftThreshold( Matrix<Field>& A, Base<Field> tau, bool relative )
{
    DEBUG_CSE
    if( relative )
        tau *= MaxNorm(A);
    auto softThresh = [&]( Field alpha ) { return SoftThreshold(alpha,tau); };
    EntrywiseMap( A, function<Field(Field)>(softThresh) );
}

template<typename Field>
void SoftThreshold
( AbstractDistMatrix<Field>& A, Base<Field> tau, bool relative )
{
    DEBUG_CSE
    if( relative )
        tau *= MaxNorm(A);
    auto softThresh = [&]( Field alpha ) { return SoftThreshold(alpha,tau); };
    EntrywiseMap( A, function<Field(Field)>(softThresh) );
}

#define PROTO(Field) \
  template Field SoftThreshold( Field alpha, Base<Field> tau ); \
  template void SoftThreshold \
  ( Matrix<Field>& A, Base<Field> tau, bool relative ); \
  template void SoftThreshold \
  ( AbstractDistMatrix<Field>& A, Base<Field> tau, bool relative );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
