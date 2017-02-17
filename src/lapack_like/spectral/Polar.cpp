/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Polar/QDWH.hpp"
#include "./Polar/SVD.hpp"

namespace El {

// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is
// Hermitian positive semi-definite. On exit, A is overwritten with Q.

template<typename Field>
PolarInfo Polar( Matrix<Field>& A, const PolarCtrl& ctrl )
{
    EL_DEBUG_CSE
    PolarInfo info;
    if( ctrl.qdwh )
        info.qdwhInfo = polar::QDWH( A, ctrl.qdwhCtrl );
    else
        polar::SVD( A );
    return info;
}

template<typename Field>
PolarInfo Polar( AbstractDistMatrix<Field>& A, const PolarCtrl& ctrl )
{
    EL_DEBUG_CSE
    PolarInfo info;
    if( ctrl.qdwh )
        info.qdwhInfo = polar::QDWH( A, ctrl.qdwhCtrl );
    else
        polar::SVD( A );
    return info;
}

template<typename Field>
PolarInfo Polar( Matrix<Field>& A, Matrix<Field>& P, const PolarCtrl& ctrl )
{
    EL_DEBUG_CSE
    PolarInfo info;
    if( ctrl.qdwh )
        info.qdwhInfo = polar::QDWH( A, P, ctrl.qdwhCtrl );
    else
        polar::SVD( A, P );
    return info;
}

template<typename Field>
PolarInfo Polar
( AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& P,
  const PolarCtrl& ctrl )
{
    EL_DEBUG_CSE
    PolarInfo info;
    if( ctrl.qdwh )
        info.qdwhInfo = polar::QDWH( A, P, ctrl.qdwhCtrl );
    else
        polar::SVD( A, P );
    return info;
}

template<typename Field>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  Matrix<Field>& A,
  const PolarCtrl& ctrl )
{
    EL_DEBUG_CSE
    PolarInfo info;
    if( ctrl.qdwh )
        info.qdwhInfo = herm_polar::QDWH( uplo, A, ctrl.qdwhCtrl );
    else
        HermitianSign( uplo, A );
    return info;
}

template<typename Field>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& A,
  const PolarCtrl& ctrl )
{
    EL_DEBUG_CSE
    PolarInfo info;
    if( ctrl.qdwh )
        info.qdwhInfo = herm_polar::QDWH( uplo, A, ctrl.qdwhCtrl );
    else
        HermitianSign( uplo, A );
    return info;
}

template<typename Field>
PolarInfo HermitianPolar
( UpperOrLower uplo, Matrix<Field>& A, Matrix<Field>& P, const PolarCtrl& ctrl )
{
    EL_DEBUG_CSE
    PolarInfo info;
    if( ctrl.qdwh )
        info.qdwhInfo = herm_polar::QDWH( uplo, A, P, ctrl.qdwhCtrl );
    else
        HermitianSign( uplo, A, P );
    return info;
}

template<typename Field>
PolarInfo HermitianPolar
( UpperOrLower uplo,
  AbstractDistMatrix<Field>& A,
  AbstractDistMatrix<Field>& P,
  const PolarCtrl& ctrl )
{
    EL_DEBUG_CSE
    PolarInfo info;
    if( ctrl.qdwh )
        info.qdwhInfo = herm_polar::QDWH( uplo, A, P, ctrl.qdwhCtrl );
    else
        HermitianSign( uplo, A, P );
    return info;
}

#define PROTO(Field) \
  template PolarInfo HermitianPolar \
  ( UpperOrLower uplo, \
    Matrix<Field>& A, \
    const PolarCtrl& ctrl ); \
  template PolarInfo HermitianPolar \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<Field>& A, \
    const PolarCtrl& ctrl ); \
  template PolarInfo HermitianPolar \
  ( UpperOrLower uplo, \
    Matrix<Field>& A, \
    Matrix<Field>& P, \
    const PolarCtrl& ctrl ); \
  template PolarInfo HermitianPolar \
  ( UpperOrLower uplo, \
    AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Field>& P, \
    const PolarCtrl& ctrl ); \
  template PolarInfo Polar \
  ( Matrix<Field>& A, \
    const PolarCtrl& ctrl ); \
  template PolarInfo Polar \
  ( AbstractDistMatrix<Field>& A, \
    const PolarCtrl& ctrl ); \
  template PolarInfo Polar \
  ( Matrix<Field>& A, \
    Matrix<Field>& P, \
    const PolarCtrl& ctrl ); \
  template PolarInfo Polar \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Field>& P, \
    const PolarCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
