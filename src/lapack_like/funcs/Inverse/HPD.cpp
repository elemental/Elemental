/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./HPD/CholeskyLVar2.hpp"
#include "./HPD/CholeskyUVar2.hpp"

namespace El {

template<typename Field>
void HPDInverse( UpperOrLower uplo, Matrix<Field>& A )
{
    EL_DEBUG_CSE
    if( uplo == LOWER )
        hpd_inv::CholeskyLVar2( A );
    else
        hpd_inv::CholeskyUVar2( A );
}

template<typename Field>
void HPDInverse( UpperOrLower uplo, AbstractDistMatrix<Field>& A )
{
    EL_DEBUG_CSE
    if( uplo == LOWER )
        hpd_inv::CholeskyLVar2( A );
    else
        hpd_inv::CholeskyUVar2( A );
}

template<typename Field>
void LocalHPDInverse( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A )
{
    EL_DEBUG_CSE
    HPDInverse( uplo, A.Matrix() );
}

#define PROTO(Field) \
  template void HPDInverse( UpperOrLower uplo, Matrix<Field>& A ); \
  template void HPDInverse( UpperOrLower uplo, AbstractDistMatrix<Field>& A ); \
  template void LocalHPDInverse \
  ( UpperOrLower uplo, DistMatrix<Field,STAR,STAR>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
