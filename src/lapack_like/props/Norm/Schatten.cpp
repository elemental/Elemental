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
Base<F> SchattenNorm( const Matrix<F>& A, Base<F> p )
{
    DEBUG_CSE
    const Int minDim = Min(A.Height(),A.Width());
    return KyFanSchattenNorm( A, minDim, p );
}

template<typename F>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_CSE
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename F>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const Matrix<F>& A, Base<F> p )
{
    DEBUG_CSE
    const Int minDim = A.Height();
    return SymmetricKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename F> 
Base<F> SchattenNorm( const ElementalMatrix<F>& A, Base<F> p )
{
    DEBUG_CSE
    const Int minDim = Min(A.Height(),A.Width());
    return KyFanSchattenNorm( A, minDim, p );
}

template<typename F>
Base<F> HermitianSchattenNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Base<F> p )
{
    DEBUG_CSE
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

template<typename F>
Base<F> SymmetricSchattenNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A, Base<F> p )
{
    DEBUG_CSE
    const Int minDim = A.Height();
    return HermitianKyFanSchattenNorm( uplo, A, minDim, p );
}

#define PROTO(F) \
  template Base<F> SchattenNorm( const Matrix<F>& A, Base<F> p ); \
  template Base<F> SchattenNorm( const ElementalMatrix<F>& A, Base<F> p ); \
  template Base<F> HermitianSchattenNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> p ); \
  template Base<F> HermitianSchattenNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricSchattenNorm \
  ( UpperOrLower uplo, const Matrix<F>& A, Base<F> p ); \
  template Base<F> SymmetricSchattenNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A, Base<F> p );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
