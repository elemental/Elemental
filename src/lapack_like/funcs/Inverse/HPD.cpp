/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./HPD/CholeskyLVar2.hpp"
#include "./HPD/CholeskyUVar2.hpp"

namespace El {

template<typename F>
void HPDInverse( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HPDInverse"))
    if( uplo == LOWER )
        hpd_inv::CholeskyLVar2( A );
    else
        hpd_inv::CholeskyUVar2( A );
}

template<typename F>
void HPDInverse( UpperOrLower uplo, ElementalMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HPDInverse"))
    if( uplo == LOWER )
        hpd_inv::CholeskyLVar2( A );
    else
        hpd_inv::CholeskyUVar2( A );
}

template<typename F>
void LocalHPDInverse( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A )
{
    DEBUG_ONLY(CSE cse("LocalHPDInverse"))
    HPDInverse( uplo, A.Matrix() );
}

#define PROTO(F) \
  template void HPDInverse( UpperOrLower uplo, Matrix<F>& A ); \
  template void HPDInverse( UpperOrLower uplo, ElementalMatrix<F>& A ); \
  template void LocalHPDInverse \
  ( UpperOrLower uplo, DistMatrix<F,STAR,STAR>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
