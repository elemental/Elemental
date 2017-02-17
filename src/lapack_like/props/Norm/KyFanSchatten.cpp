/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
Base<Field> KyFanSchattenNorm( const Matrix<Field>& A, Int k, Base<Field> p )
{
    EL_DEBUG_CSE
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFanSchatten norm");

    typedef Base<Field> Real;
    Matrix<Real> s;
    SVD( A, s );

    Real sum = 0;
    for( Int j=k-1; j>=0; --j )
        sum += Pow( s(j), p );
    return Pow( sum, 1/p );
}

template<typename Field>
Base<Field> HermitianKyFanSchattenNorm
( UpperOrLower uplo, const Matrix<Field>& A, Int k, Base<Field> p )
{
    EL_DEBUG_CSE
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFanSchatten norm");

    typedef Base<Field> Real;
    Matrix<Real> s;
    HermitianSVD( uplo, A, s );

    Real sum = 0;
    for( Int j=k-1; j>=0; --j )
        sum += Pow( s(j), p );
    return Pow( sum, 1/p );
}

template<typename Field>
Base<Field> SymmetricKyFanSchattenNorm
( UpperOrLower uplo, const Matrix<Field>& A, Int k, Base<Field> p )
{
    EL_DEBUG_CSE
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFanSchatten norm");

    typedef Base<Field> Real;
    Matrix<Field> B( A );
    Matrix<Real> s;
    MakeSymmetric( uplo, B );
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    SVD( B, s, ctrl );

    Real sum = 0;
    for( Int j=k-1; j>=0; --j )
        sum += Pow( s(j), p );
    return Pow( sum, 1/p );
}

template<typename Field>
Base<Field>
KyFanSchattenNorm( const AbstractDistMatrix<Field>& A, Int k, Base<Field> p )
{
    EL_DEBUG_CSE
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFanSchatten norm");

    typedef Base<Field> Real;
    DistMatrix<Field> B( A );
    DistMatrix<Real,VR,STAR> s( A.Grid() );
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    SVD( B, s, ctrl );

    Real localSum = 0;
    auto sTop = s( IR(0,k), ALL );
    Matrix<Real>& sTopLoc = sTop.Matrix();
    const Int localHeight = sTop.LocalHeight();
    for( Int j=localHeight-1; j>=0; --j )
        localSum += Pow( sTopLoc(j), p );
    const Real sum = mpi::AllReduce( localSum, sTop.ColComm() );
    return Pow( sum, 1/p );
}

template<typename Field>
Base<Field> HermitianKyFanSchattenNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Int k, Base<Field> p )
{
    EL_DEBUG_CSE
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFanSchatten norm");

    typedef Base<Field> Real;
    DistMatrix<Real,VR,STAR> s( A.Grid() );
    HermitianSVD( uplo, A, s );

    Real localSum = 0;
    auto sTop = s( IR(0,k), ALL );
    Matrix<Real>& sTopLoc = sTop.Matrix();
    const Int localHeight = sTop.LocalHeight();
    for( Int j=localHeight-1; j>=0; --j )
        localSum += Pow( sTopLoc(j), p );
    const Real sum = mpi::AllReduce( localSum, sTop.ColComm() );
    return Pow( sum, 1/p );
}

template<typename Field>
Base<Field> SymmetricKyFanSchattenNorm
( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, Int k, Base<Field> p )
{
    EL_DEBUG_CSE
    if( k < 1 || k > Min(A.Height(),A.Width()) )
        LogicError("Invalid index of KyFanSchatten norm");

    typedef Base<Field> Real;
    DistMatrix<Field> B( A );
    DistMatrix<Real,VR,STAR> s( A.Grid() );
    MakeSymmetric( uplo, B );
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    SVD( B, s, ctrl );

    Real localSum = 0;
    auto sTop = s( IR(0,k), ALL );
    Matrix<Real>& sTopLoc = sTop.Matrix();
    const Int localHeight = sTop.LocalHeight();
    for( Int j=localHeight-1; j>=0; --j )
        localSum += Pow( sTopLoc(j), p );
    const Real sum = mpi::AllReduce( localSum, sTop.ColComm() );
    return Pow( sum, 1/p );
}

#define PROTO(Field) \
  template Base<Field> KyFanSchattenNorm \
  ( const Matrix<Field>& A, Int k, Base<Field> p ); \
  template Base<Field> KyFanSchattenNorm \
  ( const AbstractDistMatrix<Field>& A, Int k, Base<Field> p ); \
  template Base<Field> HermitianKyFanSchattenNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A, Int k, Base<Field> p ); \
  template Base<Field> HermitianKyFanSchattenNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, \
    Int k, Base<Field> p ); \
  template Base<Field> SymmetricKyFanSchattenNorm \
  ( UpperOrLower uplo, const Matrix<Field>& A, Int k, Base<Field> p ); \
  template Base<Field> SymmetricKyFanSchattenNorm \
  ( UpperOrLower uplo, const AbstractDistMatrix<Field>& A, \
    Int k, Base<Field> p );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
