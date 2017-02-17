/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./QR/ApplyQ.hpp"
#include "./QR/BusingerGolub.hpp"
#include "./QR/Cholesky.hpp"
#include "./QR/Householder.hpp"
#include "./QR/SolveAfter.hpp"
#include "./QR/Explicit.hpp"

#include "./QR/ColSwap.hpp"

#include "./QR/TS.hpp"

namespace El {

template<typename F>
void QR
( Matrix<F>& A,
  Matrix<F>& householderScalars,
  Matrix<Base<F>>& signature )
{
    EL_DEBUG_CSE
    qr::Householder( A, householderScalars, signature );
}

namespace qr {

// TODO(poulson): Provide external wrappers
template<typename F,typename=EnableIf<IsBlasScalar<F>>>
void ScaLAPACKHelper
( DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<F,MR,STAR,BLOCK>& householderScalars )
{
    EL_DEBUG_CSE
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    householderScalars.AlignWith( A );
    householderScalars.Resize( minDim, 1 );

    auto descA = FillDesc( A );
    scalapack::QR
    ( m, n, A.Buffer(), descA.data(), householderScalars.Buffer() );
#endif
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
void ScaLAPACKHelper
( DistMatrix<F,MC,MR,BLOCK>& A,
  DistMatrix<F,MR,STAR,BLOCK>& householderScalars )
{
    EL_DEBUG_CSE
    LogicError("ScaLAPACK does not support ",TypeName<F>());
}

} // namespace qr

template<typename F>
void QR
( AbstractDistMatrix<F>& A,
  AbstractDistMatrix<F>& householderScalars,
  AbstractDistMatrix<Base<F>>& signature )
{
    EL_DEBUG_CSE
    qr::Householder( A, householderScalars, signature );
}

// Variants which perform (Businger-Golub) column-pivoting
// =======================================================

template<typename F>
void QR
( Matrix<F>& A,
  Matrix<F>& householderScalars,
  Matrix<Base<F>>& signature,
  Permutation& Omega,
  const QRCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    qr::BusingerGolub( A, householderScalars, signature, Omega, ctrl );
}

template<typename F>
void QR
( AbstractDistMatrix<F>& A,
  AbstractDistMatrix<F>& householderScalars,
  AbstractDistMatrix<Base<F>>& signature,
  DistPermutation& Omega,
  const QRCtrl<Base<F>>& ctrl )
{
    EL_DEBUG_CSE
    qr::BusingerGolub( A, householderScalars, signature, Omega, ctrl );
}

#define PROTO(F) \
  template void QR \
  ( Matrix<F>& A, \
    Matrix<F>& householderScalars, \
    Matrix<Base<F>>& signature ); \
  template void QR \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& householderScalars, \
    AbstractDistMatrix<Base<F>>& signature ); \
  template void QR \
  ( Matrix<F>& A, \
    Matrix<F>& householderScalars, \
    Matrix<Base<F>>& signature, \
    Permutation& Omega, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void QR \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& householderScalars, \
    AbstractDistMatrix<Base<F>>& signature, \
    DistPermutation& Omega, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::ExplicitTriang \
  ( Matrix<F>& A, const QRCtrl<Base<F>>& ctrl ); \
  template void qr::ExplicitTriang \
  ( AbstractDistMatrix<F>& A, const QRCtrl<Base<F>>& ctrl ); \
  template void qr::ExplicitUnitary \
  ( Matrix<F>& A, bool thinQR, const QRCtrl<Base<F>>& ctrl ); \
  template void qr::ExplicitUnitary \
  ( AbstractDistMatrix<F>& A, bool thinQR, const QRCtrl<Base<F>>& ctrl ); \
  template void qr::Explicit \
  ( Matrix<F>& A, \
    Matrix<F>& R, \
    bool thinQR, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::Explicit \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& R, \
    bool thinQR, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::Explicit \
  ( Matrix<F>& A, \
    Matrix<F>& R, \
    Matrix<Int>& Omega, \
    bool thinQR, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::Explicit \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& R, \
    AbstractDistMatrix<Int>& Omega, \
    bool thinQR, \
    const QRCtrl<Base<F>>& ctrl ); \
  template void qr::NeighborColSwap \
  (       Matrix<F>& Q, \
          Matrix<F>& R, \
    Int j ); \
  template void qr::DisjointNeighborColSwaps \
  (       Matrix<F>& Q, \
          Matrix<F>& R, \
    const Matrix<Int>& colSwaps ); \
  template void qr::ApplyQ \
  ( LeftOrRight side, \
    Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
    const Matrix<Base<F>>& signature, \
          Matrix<F>& B ); \
  template void qr::ApplyQ \
  ( LeftOrRight side, \
    Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
    const AbstractDistMatrix<Base<F>>& signature, \
          AbstractDistMatrix<F>& B ); \
  template void qr::SolveAfter \
  ( Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& householderScalars, \
    const Matrix<Base<F>>& signature, \
    const Matrix<F>& B, \
          Matrix<F>& X ); \
  template void qr::SolveAfter \
  ( Orientation orientation, \
    const AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& householderScalars, \
    const AbstractDistMatrix<Base<F>>& signature, \
    const AbstractDistMatrix<F>& B, \
          AbstractDistMatrix<F>& X ); \
  template void qr::Cholesky \
  ( Matrix<F>& A, \
    Matrix<F>& R ); \
  template void qr::Cholesky \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& R ); \
  template qr::TreeData<F> qr::TS( const AbstractDistMatrix<F>& A ); \
  template void qr::ExplicitTS \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<F>& R ); \
  template Matrix<F>& qr::ts::RootQR \
  ( const AbstractDistMatrix<F>& A, TreeData<F>& treeData ); \
  template const Matrix<F>& qr::ts::RootQR \
  ( const AbstractDistMatrix<F>& A, const TreeData<F>& treeData ); \
  template void qr::ts::Reduce \
  ( const AbstractDistMatrix<F>& A, TreeData<F>& treeData ); \
  template void qr::ts::Scatter \
  ( AbstractDistMatrix<F>& A, const TreeData<F>& treeData );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
