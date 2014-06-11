/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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
void QR( DistMatrix<F>& A )
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
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    qr::Householder( A, t, d );
}

// Variants which perform (Businger-Golub) column-pivoting
// =======================================================

template<typename F> 
Int QR
( Matrix<F>& A, Matrix<Int>& pPerm, 
  const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    return qr::BusingerGolub( A, pPerm, ctrl );
}

template<typename F,Dist UPerm> 
Int QR
( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& pPerm, 
  const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    return qr::BusingerGolub( A, pPerm, ctrl );
}

template<typename F> 
Int QR
( Matrix<F>& A, Matrix<F>& t, 
  Matrix<Base<F>>& d, Matrix<Int>& pPerm, 
  const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    return qr::BusingerGolub( A, t, d, pPerm, ctrl );
}

template<typename F,Dist UPerm> 
Int QR
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, 
  DistMatrix<Base<F>,MD,STAR>& d, DistMatrix<Int,UPerm,STAR>& pPerm,
  const QRCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("QR"))
    return qr::BusingerGolub( A, t, d, pPerm, ctrl );
}

#define PROTO_DIST(F,U) \
   template qr::TreeData<F> qr::TS( const DistMatrix<F,U,STAR>& A ); \
  template void qr::ExplicitTS \
  ( DistMatrix<F,U,STAR>& A, DistMatrix<F,STAR,STAR>& R ); \
  template Matrix<F>& qr::ts::RootQR \
  ( const DistMatrix<F,U,STAR>& A, TreeData<F>& treeData ); \
  template const Matrix<F>& qr::ts::RootQR \
  ( const DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData ); \
  template void qr::ts::Reduce \
  ( const DistMatrix<F,U,STAR>& A, TreeData<F>& treeData ); \
  template void qr::ts::Scatter \
  ( DistMatrix<F,U,STAR>& A, const TreeData<F>& treeData ); 

#define PROTO(F) \
  template void QR( Matrix<F>& A ); \
  template void QR( DistMatrix<F>& A ); \
  template void QR( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d ); \
  template void QR \
  ( DistMatrix<F>& A, \
    DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d ); \
  template Int QR \
  ( Matrix<F>& A, Matrix<Int>& pPerm, \
    const QRCtrl<Base<F>> ctrl ); \
  template Int QR \
  ( DistMatrix<F>& A, DistMatrix<Int,VR,STAR>& pPerm, \
    const QRCtrl<Base<F>> ctrl ); \
  template Int QR \
  ( Matrix<F>& A, Matrix<F>& t, Matrix<Base<F>>& d, \
    Matrix<Int>& pPerm, const QRCtrl<Base<F>> ctrl ); \
  template Int QR \
  ( DistMatrix<F>& A, \
    DistMatrix<F,MD,STAR>& t, DistMatrix<Base<F>,MD,STAR>& d, \
    DistMatrix<Int,VR,STAR>& pPerm, const QRCtrl<Base<F>> ctrl ); \
  template void qr::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, Matrix<F>& B ); \
  template void qr::ApplyQ \
  ( LeftOrRight side, Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, \
    const DistMatrix<Base<F>,MD,STAR>& d, DistMatrix<F>& B ); \
  template void qr::Cholesky \
  ( Matrix<F>& A, Matrix<F>& R ); \
  template void qr::Cholesky \
  ( DistMatrix<F,VC,STAR>& A, DistMatrix<F,STAR,STAR>& R ); \
  template void qr::Explicit( Matrix<F>& A, bool colPiv ); \
  template void qr::Explicit( DistMatrix<F>& A, bool colPiv ); \
  template void qr::Explicit \
  ( Matrix<F>& A, Matrix<F>& R, bool colPiv ); \
  template void qr::Explicit \
  ( DistMatrix<F>& A, DistMatrix<F>& R, bool colPiv ); \
  template void qr::Explicit \
  ( Matrix<F>& A, Matrix<F>& R, Matrix<Int>& pPerm ); \
  template void qr::Explicit \
  ( DistMatrix<F>& A, DistMatrix<F>& R, DistMatrix<Int,VR,STAR>& pPerm ); \
  template void qr::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, const Matrix<F>& t, \
    const Matrix<Base<F>>& d, const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void qr::SolveAfter \
  ( Orientation orientation, \
    const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& t, \
    const DistMatrix<Base<F>,MD,STAR>& d, const DistMatrix<F>& B, \
          DistMatrix<F>& X ); \
  PROTO_DIST(F,MC  ) \
  PROTO_DIST(F,MD  ) \
  PROTO_DIST(F,MR  ) \
  PROTO_DIST(F,STAR) \
  PROTO_DIST(F,VC  ) \
  PROTO_DIST(F,VR  )

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
