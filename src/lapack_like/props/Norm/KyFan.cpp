/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
Base<F> KyFanNorm( const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CSE cse("KyFanNorm"))
    return KyFanSchattenNorm( A, k, Base<F>(1) );
}

template<typename F>
Base<F> HermitianKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CSE cse("HermitianKyFanNorm"))
    return HermitianKyFanSchattenNorm( uplo, A, k, Base<F>(1) );
}

template<typename F>
Base<F> SymmetricKyFanNorm( UpperOrLower uplo, const Matrix<F>& A, Int k )
{
    DEBUG_ONLY(CSE cse("SymmetricKyFanNorm"))
    return SymmetricKyFanSchattenNorm( uplo, A, k, Base<F>(1) );
}

template<typename F> 
Base<F> KyFanNorm( const ElementalMatrix<F>& A, Int k )
{
    DEBUG_ONLY(CSE cse("KyFanNorm"))
    return KyFanSchattenNorm( A, k, Base<F>(1) );
}

template<typename F>
Base<F> HermitianKyFanNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Int k )
{
    DEBUG_ONLY(CSE cse("HermitianKyFanNorm"))
    return HermitianKyFanSchattenNorm( uplo, A, k, Base<F>(1) );
}

template<typename F>
Base<F> SymmetricKyFanNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Int k )
{
    DEBUG_ONLY(CSE cse("SymmetricKyFanNorm"))
    return SymmetricKyFanSchattenNorm( uplo, A, k, Base<F>(1) );
}

#define PROTO(F) \
  template Base<F> KyFanNorm( const Matrix<F>& A, Int k ); \
  template Base<F> KyFanNorm( const ElementalMatrix<F>& A, Int k ); \
  template Base<F> HermitianKyFanNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Int k ); \
  template Base<F> HermitianKyFanNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A, Int k ); \
  template Base<F> SymmetricKyFanNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Int k ); \
  template Base<F> SymmetricKyFanNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A, Int k );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
