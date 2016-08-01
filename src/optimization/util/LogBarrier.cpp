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
Base<F> LogBarrier( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_CSE
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
{
    DEBUG_CSE
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

template<typename F> 
Base<F> LogBarrier( UpperOrLower uplo, const ElementalMatrix<F>& A )
{
    DEBUG_CSE
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename F> 
Base<F> LogBarrier
( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite )
{
    DEBUG_CSE
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

#define PROTO(F) \
  template Base<F> LogBarrier( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> LogBarrier \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A ); \
  template Base<F> LogBarrier \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template Base<F> LogBarrier \
  ( UpperOrLower uplo, ElementalMatrix<F>& A, bool canOverwrite );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
