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
Base<F> NuclearNorm( const Matrix<F>& A )
{
    DEBUG_CSE
    return SchattenNorm( A, Base<F>(1) );
}

template<typename F>
Base<F> HermitianNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_CSE
    return HermitianSchattenNorm( uplo, A, Base<F>(1) );
}

template<typename F>
Base<F> SymmetricNuclearNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_CSE
    return SymmetricSchattenNorm( uplo, A, Base<F>(1) );
}

template<typename F> 
Base<F> NuclearNorm( const ElementalMatrix<F>& A )
{
    DEBUG_CSE
    return SchattenNorm( A, Base<F>(1) );
}

template<typename F>
Base<F> HermitianNuclearNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A )
{
    DEBUG_CSE
    return HermitianSchattenNorm( uplo, A, Base<F>(1) );
}

template<typename F>
Base<F> SymmetricNuclearNorm
( UpperOrLower uplo, const ElementalMatrix<F>& A )
{
    DEBUG_CSE
    return SymmetricSchattenNorm( uplo, A, Base<F>(1) );
}

#define PROTO(F) \
  template Base<F> NuclearNorm( const Matrix<F>& A ); \
  template Base<F> NuclearNorm( const ElementalMatrix<F>& A ); \
  template Base<F> HermitianNuclearNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianNuclearNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A ); \
  template Base<F> SymmetricNuclearNorm \
  ( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricNuclearNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
