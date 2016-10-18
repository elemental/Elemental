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
Base<F> TwoNorm( const Matrix<F>& A )
{
    DEBUG_CSE
    Matrix<Base<F>> s;
    SVD( A, s );
    return InfinityNorm( s );
}

template<typename F>
Base<F> HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_CSE
    Matrix<Base<F>> s;
    HermitianSVD( uplo, A, s );
    return InfinityNorm( s );
}

template<typename F>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A )
{
    DEBUG_CSE
    Matrix<F> B( A );
    Matrix<Base<F>> s;
    MakeSymmetric( uplo, B );
    SVDCtrl<Base<F>> ctrl;
    ctrl.overwrite = true;
    SVD( B, s, ctrl );
    return MaxNorm( s );
}

template<typename F> 
Base<F> TwoNorm( const ElementalMatrix<F>& A )
{
    DEBUG_CSE
    DistMatrix<Base<F>,VR,STAR> s( A.Grid() );
    SVD( A, s );
    return InfinityNorm( s );
}

template<typename F>
Base<F> HermitianTwoNorm( UpperOrLower uplo, const ElementalMatrix<F>& A )
{
    DEBUG_CSE
    DistMatrix<Base<F>,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, A, s );
    return InfinityNorm( s );
}

template<typename F>
Base<F> SymmetricTwoNorm( UpperOrLower uplo, const ElementalMatrix<F>& A )
{
    DEBUG_CSE
    DistMatrix<F> B( A );
    DistMatrix<Base<F>,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVDCtrl<Base<F>> ctrl;
    ctrl.overwrite = true;
    SVD( B, s, ctrl );
    return MaxNorm( s );
}

#define PROTO(F) \
  template Base<F> TwoNorm( const Matrix<F>& A ); \
  template Base<F> TwoNorm( const ElementalMatrix<F>& A ); \
  template Base<F> HermitianTwoNorm( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> HermitianTwoNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A ); \
  template Base<F> SymmetricTwoNorm( UpperOrLower uplo, const Matrix<F>& A ); \
  template Base<F> SymmetricTwoNorm \
  ( UpperOrLower uplo, const ElementalMatrix<F>& A );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
