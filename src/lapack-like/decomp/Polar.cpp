/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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
void Polar( DistMatrix<F>& A, const PolarCtrl& ctrl )
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
void Polar( DistMatrix<F>& A, DistMatrix<F>& P, const PolarCtrl& ctrl )
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
( UpperOrLower uplo, DistMatrix<F>& A, const PolarCtrl& ctrl )
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
( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& P, const PolarCtrl& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    if( ctrl.qdwh )
        ctrl.numIts = herm_polar::QDWH( uplo, A, P, ctrl );
    else
        HermitianSign( uplo, A, P );
}

#define PROTO(F) \
  template void Polar( Matrix<F>& A, const PolarCtrl& ctrl ); \
  template void Polar( DistMatrix<F>& A, const PolarCtrl& ctrl ); \
  template void Polar \
  ( Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl ); \
  template void Polar \
  ( DistMatrix<F>& A, DistMatrix<F>& P, const PolarCtrl& ctrl ); \
  template void HermitianPolar \
  ( UpperOrLower uplo, Matrix<F>& A, const PolarCtrl& ctrl ); \
  template void HermitianPolar \
  ( UpperOrLower uplo, DistMatrix<F>& A, const PolarCtrl& ctrl ); \
  template void HermitianPolar \
  ( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P, const PolarCtrl& ctrl ); \
  template void HermitianPolar \
  ( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& P, const PolarCtrl& ctrl );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
