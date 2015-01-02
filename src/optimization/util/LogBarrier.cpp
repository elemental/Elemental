/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

template<typename F> 
Base<F> LogBarrier( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename F> 
Base<F> LogBarrier
( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool canOverwrite )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

#define PROTO(F) \
  template Base<F> LogBarrier( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> LogBarrier \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> LogBarrier \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template Base<F> LogBarrier \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, bool canOverwrite );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
