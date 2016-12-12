/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

// The Frobenius norm prox returns the solution to
//     arg min || A ||_F + tau/2 || A - A0 ||_F^2
//        A

namespace El {

template<typename Field>
void FrobeniusProx( Matrix<Field>& A, const Base<Field>& tau )
{
    EL_DEBUG_CSE
    const Base<Field> frobNorm = FrobeniusNorm( A );
    if( frobNorm > 1/tau )
        A *= 1-1/(tau*frobNorm);
    else
        Zero( A );
}

template<typename Field>
void FrobeniusProx( AbstractDistMatrix<Field>& A, const Base<Field>& tau )
{
    EL_DEBUG_CSE
    const Base<Field> frobNorm = FrobeniusNorm( A );
    if( frobNorm > 1/tau )
        A *= 1-1/(tau*frobNorm);
    else
        Zero( A );
}

#define PROTO(Field) \
  template void FrobeniusProx \
  ( Matrix<Field>& A, const Base<Field>& tau ); \
  template void FrobeniusProx \
  ( AbstractDistMatrix<Field>& A, const Base<Field>& tau );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
