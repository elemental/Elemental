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
Field SoftThreshold( const Field& alpha, const Base<Field>& tau )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( tau < 0 )
          LogicError("Negative threshold does not make sense");
    )
    const Base<Field> scale = Abs(alpha);
    return ( scale <= tau ? Field(0) : alpha-(alpha/scale)*tau );
}

template<typename Field>
void SoftThreshold( Matrix<Field>& A, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    Base<Field> tauMod = tau;
    if( relative )
        tauMod *= MaxNorm(A);
    auto softThresh =
      [&]( const Field& alpha ) { return SoftThreshold(alpha,tauMod); };
    EntrywiseMap( A, MakeFunction(softThresh) );
}

template<typename Field>
void SoftThreshold
( AbstractDistMatrix<Field>& A, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    Base<Field> tauMod = tau;
    if( relative )
        tauMod *= MaxNorm(A);
    auto softThresh =
      [&]( const Field& alpha ) { return SoftThreshold(alpha,tauMod); };
    EntrywiseMap( A, MakeFunction(softThresh) );
}

#define PROTO(Field) \
  template Field SoftThreshold \
  ( const Field& alpha, const Base<Field>& tau ); \
  template void SoftThreshold \
  ( Matrix<Field>& A, const Base<Field>& tau, bool relative ); \
  template void SoftThreshold \
  ( AbstractDistMatrix<Field>& A, const Base<Field>& tau, bool relative );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
