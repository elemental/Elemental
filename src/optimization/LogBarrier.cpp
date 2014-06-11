/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename F>
Base<F> LogBarrier( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

template<typename F> 
Base<F> LogBarrier( UpperOrLower uplo, const DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A );
    return -safeDet.kappa*safeDet.n;
}

template<typename F> 
Base<F> LogBarrier
( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite=false )
{
    DEBUG_ONLY(CallStackEntry cse("LogBarrier"))
    SafeProduct<Base<F>> safeDet = SafeHPDDeterminant( uplo, A, canOverwrite );
    return -safeDet.kappa*safeDet.n;
}

#define PROTO(F) \
  template Base<F> LogBarrier( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> LogBarrier( UpperOrLower uplo, const DistMatrix<F>& A ); \
  template Base<F> LogBarrier \
  ( UpperOrLower uplo, Matrix<F>& A, bool canOverwrite ); \
  template Base<F> LogBarrier \
  ( UpperOrLower uplo, DistMatrix<F>& A, bool canOverwrite );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
