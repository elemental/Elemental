/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./General/LUPartialPiv.hpp"

namespace El {

template<typename Field>
void Inverse( Matrix<Field>& A )
{
    EL_DEBUG_CSE
    inverse::LUPartialPiv( A );
}

template<typename Field>
void Inverse( AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    inverse::LUPartialPiv( A );
}

template<typename Field>
void LocalInverse( DistMatrix<Field,STAR,STAR>& A )
{
    EL_DEBUG_CSE
    Inverse( A.Matrix() );
}

#define PROTO(Field) \
  template void Inverse( Matrix<Field>& A ); \
  template void Inverse( AbstractDistMatrix<Field>& A ); \
  template void LocalInverse( DistMatrix<Field,STAR,STAR>& A ); \
  template void inverse::AfterLUPartialPiv \
  ( Matrix<Field>& A, const Permutation& P ); \
  template void inverse::AfterLUPartialPiv \
  ( AbstractDistMatrix<Field>& A, const DistPermutation& P );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
