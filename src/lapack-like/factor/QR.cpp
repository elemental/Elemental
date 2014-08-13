/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./QR/ApplyQ.hpp"
#include "./QR/BusingerGolub.hpp"
#include "./QR/Cholesky.hpp"
#include "./QR/Householder.hpp"
#include "./QR/SolveAfter.hpp"
#include "./QR/Explicit.hpp"
#include "./QR/TS.hpp"

namespace El {

template<typename F> 
void QR( Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A );
}

template<typename F> 
void QR( AbstractDistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A );
}

template<typename F> 
void QR( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A, t, d );
}

template<typename F> 
void QR
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, 
  AbstractDistMatrix<Base<F>>& d )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A, t, d );
}

// Variants which perform (Businger-Golub) column-pivoting
// =======================================================

template<typename F> 
Int QR
( Matrix<F>& A, Matrix<Int>& p, const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    return qr::BusingerGolub( A, p, ctrl );
}

template<typename F> 
Int QR
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p, 
  const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    return qr::BusingerGolub( A, p, ctrl );
}

template<typename F> 
Int QR
( Matrix<F>& A, Matrix<F>& t, 
  Matrix<Base<F>>& d, Matrix<Int>& p, const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    return qr::BusingerGolub( A, t, d, p, ctrl );
}

template<typename F> 
Int QR
( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, 
  AbstractDistMatrix<Base<F>>& d, AbstractDistMatrix<Int>& p,
  const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    return qr::BusingerGolub( A, t, d, p, ctrl );
}

#define PROTO(F) \
  template void QR( Matrix<F>& A ); \
  template void QR( AbstractDistMatrix<F>& A ); \
  template void QR( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d ); \
  template void QR \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, \
    AbstractDistMatrix<Base<F>>& d ); \
  template Int QR \
  ( Matrix<F>& A, Matrix<Int>& p, \
    const QRCtrl<Base<F>> ctrl ); \
  template Int QR \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Int>& p, \
    const QRCtrl<Base<F>> ctrl ); \
  template Int QR \
  ( Matrix<F>& A, Matrix<F>& t, \
    Matrix<Base<F>>& d, Matrix<Int>& p, \
    const QRCtrl<Base<F>> ctrl ); \
  template Int QR \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& t, \
    AbstractDistMatrix<Base<F>>& d, AbstractDistMatrix<Int>& p, \
    const QRCtrl<Base<F>> ctrl ); \
  template void qr::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, Matrix<F>& B ); \
  template void qr::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& t, \
    const AbstractDistMatrix<Base<F>>& d, AbstractDistMatrix<F>& B ); \
  template void qr::Cholesky \
  ( Matrix<F>& A, Matrix<F>& R ); \
  template void qr::Cholesky \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R ); \
  template void qr::Explicit( Matrix<F>& A, bool colPiv ); \
  template void qr::Explicit( AbstractDistMatrix<F>& A, bool colPiv ); \
  template void qr::Explicit \
  ( Matrix<F>& A, Matrix<F>& R, bool colPiv ); \
  template void qr::Explicit \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R, bool colPiv ); \
  template void qr::Explicit \
  ( Matrix<F>& A, Matrix<F>& R, Matrix<Int>& p ); \
  template void qr::Explicit \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R, \
    AbstractDistMatrix<Int>& p ); \
  template void qr::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void qr::SolveAfter \
  ( Orientation orientation, \
    const AbstractDistMatrix<F      >& A, const AbstractDistMatrix<F>& t, \
    const AbstractDistMatrix<Base<F>>& d, const AbstractDistMatrix<F>& B, \
          AbstractDistMatrix<F      >& X ); \
  template qr::TreeData<F> qr::TS( const AbstractDistMatrix<F>& A ); \
  template void qr::ExplicitTS \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& R ); \
  template Matrix<F>& qr::ts::RootQR \
  ( const AbstractDistMatrix<F>& A, TreeData<F>& treeData ); \
  template const Matrix<F>& qr::ts::RootQR \
  ( const AbstractDistMatrix<F>& A, const TreeData<F>& treeData ); \
  template void qr::ts::Reduce \
  ( const AbstractDistMatrix<F>& A, TreeData<F>& treeData ); \
  template void qr::ts::Scatter \
  ( AbstractDistMatrix<F>& A, const TreeData<F>& treeData ); 

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
