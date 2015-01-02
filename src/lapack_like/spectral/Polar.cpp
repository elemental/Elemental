/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Polar/QDWH.hpp"
#include "./Polar/SVD.hpp"

namespace El {

// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is 
// Hermitian positive semi-definite. On exit, A is overwritten with Q.

template<typename F>
void Polar( Matrix<F>& A, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Polar"))
    if( ctrl.qdwh )
        ctrl.numIts = polar::QDWH( A, ctrl );
    else
        polar::SVD( A );
}

template<typename F>
void Polar( AbstractDistMatrix<F>& A, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Polar"))
    if( ctrl.qdwh )
        ctrl.numIts = polar::QDWH( A, ctrl );
    else
        polar::SVD( A );
}

template<typename F>
void Polar( Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Polar"))
    if( ctrl.qdwh )
        ctrl.numIts = polar::QDWH( A, P, ctrl );
    else
        polar::SVD( A, P );
}

template<typename F>
void Polar
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& P, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Polar"))
    if( ctrl.qdwh )
        ctrl.numIts = polar::QDWH( A, P, ctrl );
    else
        polar::SVD( A, P );
}

template<typename F>
void HermitianPolar( UpperOrLower uplo, Matrix<F>& A, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    if( ctrl.qdwh )
        ctrl.numIts = herm_polar::QDWH( uplo, A, ctrl );
    else
        HermitianSign( uplo, A );
}

template<typename F>
void HermitianPolar
( UpperOrLower uplo, AbstractDistMatrix<F>& A, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    if( ctrl.qdwh )
        ctrl.numIts = herm_polar::QDWH( uplo, A, ctrl );
    else
        HermitianSign( uplo, A );
}

template<typename F>
void HermitianPolar
( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    if( ctrl.qdwh )
        ctrl.numIts = herm_polar::QDWH( uplo, A, P, ctrl );
    else
        HermitianSign( uplo, A, P );
}

template<typename F>
void HermitianPolar
( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& P, 
  const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    if( ctrl.qdwh )
        ctrl.numIts = herm_polar::QDWH( uplo, A, P, ctrl );
    else
        HermitianSign( uplo, A, P );
}

#define PROTO(F) \
  template void Polar( Matrix<F>& A, const PolarCtrl& ctrl ); \
  template void Polar( AbstractDistMatrix<F>& A, const PolarCtrl& ctrl ); \
  template void Polar \
  ( Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl ); \
  template void Polar \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& P, \
    const PolarCtrl& ctrl ); \
  template void HermitianPolar \
  ( UpperOrLower uplo, Matrix<F>& A, const PolarCtrl& ctrl ); \
  template void HermitianPolar \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, const PolarCtrl& ctrl ); \
  template void HermitianPolar \
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl ); \
  template void HermitianPolar \
  ( UpperOrLower uplo, AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& P, \
    const PolarCtrl& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
