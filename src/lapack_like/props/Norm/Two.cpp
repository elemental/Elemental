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
Base<F> TwoNorm( const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("TwoNorm"))
    Matrix<F> B( A );
    Matrix<Base<F>> s;
    SVD( B, s );
    return InfinityNorm( s );
}

template<typename F>
Base<F> HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HermitianTwoNorm"))
    Matrix<F> B( A );
    Matrix<Base<F>> s;
    HermitianSVD( uplo, B, s );
    return InfinityNorm( s );
}

template<typename F>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricTwoNorm"))
    Matrix<F> B( A );
    Matrix<Base<F>> s;
    MakeSymmetric( uplo, B );
    SVD( B, s );
    return MaxNorm( s );
}

template<typename F> 
Base<F> TwoNorm( const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("TwoNorm"))
    DistMatrix<F> B( A );
    DistMatrix<Base<F>,VR,STAR> s( A.Grid() );
    SVD( B, s );
    return InfinityNorm( s );
}

template<typename F>
Base<F> HermitianTwoNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("HermitianTwoNorm"))
    DistMatrix<F> B( A );
    DistMatrix<Base<F>,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, B, s );
    return InfinityNorm( s );
}

template<typename F>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CSE cse("SymmetricTwoNorm"))
    DistMatrix<F> B( A );
    DistMatrix<Base<F>,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVD( B, s );
    return MaxNorm( s );
}

#define PROTO(F) \
  template Base<F> TwoNorm( const Matrix<F>& A ); \
  template Base<F> TwoNorm( const AbstractDistMatrix<F>& A ); \
  template Base<F> HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianTwoNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A ); \
  template Base<F> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricTwoNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<F>& A );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
