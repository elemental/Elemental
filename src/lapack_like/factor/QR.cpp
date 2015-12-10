/*
   Copyright (c) 2009-2015, Jack Poulson
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
void QR
( Matrix<F>& A,
  Matrix<F>& t,
  Matrix<Base<F>>& d )
{
    DEBUG_ONLY(CSE cse("QR"))
    qr::Householder( A, t, d );
}

template<typename F> 
void QR
( ElementalMatrix<F>& A,
  ElementalMatrix<F>& t, 
  ElementalMatrix<Base<F>>& d )
{
    DEBUG_ONLY(CSE cse("QR"))
    qr::Householder( A, t, d );
}

template<typename F>
void QR
( DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<F,MR,STAR,BLOCK>& t )
{
    DEBUG_ONLY(CSE cse("QR"))
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    t.AlignWith( A );
    t.Resize( minDim, 1 ); 

    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );

    scalapack::QR( m, n, A.Buffer(), descA.data(), t.Buffer() );

    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
}

// Variants which perform (Businger-Golub) column-pivoting
// =======================================================

template<typename F> 
void QR
( Matrix<F>& A,
  Matrix<F>& t, 
  Matrix<Base<F>>& d,
  Permutation& Omega,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("QR"))
    qr::BusingerGolub( A, t, d, Omega, ctrl );
}

template<typename F> 
void QR
( ElementalMatrix<F>& A,
  ElementalMatrix<F>& t, 
  ElementalMatrix<Base<F>>& d,
  DistPermutation& Omega,
  const QRCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("QR"))
    qr::BusingerGolub( A, t, d, Omega, ctrl );
}

#define PROTO(F) \
  template void QR \
  ( Matrix<F>& A, \
    Matrix<F>& t, \
    Matrix<Base<F>>& d ); \
  template void QR \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& t, \
    ElementalMatrix<Base<F>>& d ); \
  template void QR \
  ( Matrix<F>& A, \
    Matrix<F>& t, \
    Matrix<Base<F>>& d, \
    Permutation& Omega, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void QR \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& t, \
    ElementalMatrix<Base<F>>& d, \
    DistPermutation& Omega, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void QR \
  ( DistMatrix<F,MC,MR,BLOCK>& A, \
    DistMatrix<F,MR,STAR,BLOCK>& t ); \
  template void qr::ExplicitTriang \
  ( Matrix<F>& A, const QRCtrl<Base<F>>& ctrl ); \
  template void qr::ExplicitTriang \
  ( ElementalMatrix<F>& A, const QRCtrl<Base<F>>& ctrl ); \
  template void qr::ExplicitUnitary \
  ( Matrix<F>& A, bool thinQR, const QRCtrl<Base<F>>& ctrl ); \
  template void qr::ExplicitUnitary \
  ( ElementalMatrix<F>& A, bool thinQR, const QRCtrl<Base<F>>& ctrl ); \
  template void qr::Explicit \
  ( Matrix<F>& A, \
    Matrix<F>& R, \
    bool thinQR, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::Explicit \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& R, \
    bool thinQR, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::Explicit \
  ( Matrix<F>& A, \
    Matrix<F>& R, \
    Matrix<Int>& Omega, \
    bool thinQR, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::Explicit \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& R, \
    ElementalMatrix<Int>& Omega, \
    bool thinQR, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::ApplyQ \
  ( LeftOrRight side, \
    Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& t, \
    const Matrix<Base<F>>& d, \
          Matrix<F>& B ); \
  template void qr::ApplyQ \
  ( LeftOrRight side, \
    Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& t, \
    const ElementalMatrix<Base<F>>& d, \
          ElementalMatrix<F>& B ); \
  template void qr::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& t, \
    const Matrix<Base<F>>& d, \
    const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void qr::SolveAfter \
  ( Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& t, \
    const ElementalMatrix<Base<F>>& d, \
    const ElementalMatrix<F>& B, \
          ElementalMatrix<F>& X ); \
  template void qr::Cholesky \
  ( Matrix<F>& A, \
    Matrix<F>& R ); \
  template void qr::Cholesky \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& R ); \
  template qr::TreeData<F> qr::TS( const ElementalMatrix<F>& A ); \
  template void qr::ExplicitTS \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<F>& R ); \
  template Matrix<F>& qr::ts::RootQR \
  ( const ElementalMatrix<F>& A, TreeData<F>& treeData ); \
  template const Matrix<F>& qr::ts::RootQR \
  ( const ElementalMatrix<F>& A, const TreeData<F>& treeData ); \
  template void qr::ts::Reduce \
  ( const ElementalMatrix<F>& A, TreeData<F>& treeData ); \
  template void qr::ts::Scatter \
  ( ElementalMatrix<F>& A, const TreeData<F>& treeData ); 

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
