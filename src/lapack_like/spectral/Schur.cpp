/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./Schur/CheckReal.hpp"
#include "./Schur/RealToComplex.hpp"
#include "./Schur/QuasiTriangEig.hpp"
#include "./Schur/QR.hpp"
#include "./Schur/SDC.hpp"
#include "./Schur/InverseFreeSDC.hpp"

namespace El {

template<typename F>
void Schur
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
    if( ctrl.useSDC )
    {
        if( fullTriangle )
        {
            Matrix<F> Q;
            schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
        }
        else
            schur::SDC( A, w, ctrl.sdcCtrl );
    }
    else
        schur::QR( A, w, fullTriangle, ctrl.time );
}

template<typename F>
void Schur
( Matrix<F>& A,
  Matrix<Complex<Base<F>>>& w,
  Matrix<F>& Q,
  bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
    if( ctrl.useSDC )
        schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
    else
        schur::QR( A, w, Q, fullTriangle, ctrl.time );
}

template<typename F>
void Schur
( ElementalMatrix<F>& A,
  ElementalMatrix<Complex<Base<F>>>& w, 
  bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
#ifdef EL_HAVE_SCALAPACK
    if( ctrl.useSDC )
    {
        if( fullTriangle )
        {
            DistMatrix<F> Q(A.Grid());
            schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
        }
        else
            schur::SDC( A, w, ctrl.sdcCtrl );
    }
    else
    {
        schur::QR( A, w, fullTriangle, ctrl.qrCtrl, ctrl.time );
    }
#else
    if( fullTriangle )
    {
        DistMatrix<F> Q(A.Grid());
        schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
    }
    else
        schur::SDC( A, w, ctrl.sdcCtrl );
#endif
}

template<typename F>
void Schur
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w, 
  bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
    schur::QR( A, w, fullTriangle, ctrl.qrCtrl, ctrl.time );
}

template<typename F>
void Schur
( ElementalMatrix<F>& A,
  ElementalMatrix<Complex<Base<F>>>& w, 
  ElementalMatrix<F>& Q,
  bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
#ifdef EL_HAVE_SCALAPACK
    if( ctrl.useSDC )
        schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
    else
        schur::QR( A, w, Q, fullTriangle, ctrl.qrCtrl, ctrl.time );
#else
    schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
#endif
}

template<typename F>
void Schur
( DistMatrix<F,MC,MR,BLOCK>& A,
  ElementalMatrix<Complex<Base<F>>>& w, 
  DistMatrix<F,MC,MR,BLOCK>& Q,
  bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_CSE
    schur::QR( A, w, Q, fullTriangle, ctrl.qrCtrl, ctrl.time );
}

#define PROTO(F) \
  template void Schur \
  ( Matrix<F>& A, \
    Matrix<Complex<Base<F>>>& w, \
    bool fullTriangle, \
    const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Complex<Base<F>>>& w, \
    bool fullTriangle, \
    const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( DistMatrix<F,MC,MR,BLOCK>& A, \
    ElementalMatrix<Complex<Base<F>>>& w, \
    bool fullTriangle, \
    const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( Matrix<F>& A, \
    Matrix<Complex<Base<F>>>& w, \
    Matrix<F>& Q, \
    bool fullTriangle, \
    const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Complex<Base<F>>>& w, \
    ElementalMatrix<F>& Q, \
    bool fullTriangle, \
    const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( DistMatrix<F,MC,MR,BLOCK>& A, \
    ElementalMatrix<Complex<Base<F>>>& w, \
    DistMatrix<F,MC,MR,BLOCK>& Q, \
    bool fullTriangle, \
    const SchurCtrl<Base<F>> ctrl ); \
  template void schur::CheckRealSchur \
  ( const Matrix<F>& U, bool standardForm ); \
  template void schur::CheckRealSchur \
  ( const ElementalMatrix<F>& U, bool standardForm ); \
  template void schur::QuasiTriangEig \
  ( const Matrix<F>& dMain, \
    const Matrix<F>& dSub, \
    const Matrix<F>& dSup, \
          Matrix<Complex<Base<F>>>& w ); \
  template void schur::QuasiTriangEig \
  ( const Matrix<F>& U, \
          Matrix<Complex<Base<F>>>& w ); \
  template Matrix<Complex<Base<F>>> schur::QuasiTriangEig \
  ( const Matrix<F>& U ); \
  template DistMatrix<Complex<Base<F>>,VR,STAR> \
  schur::QuasiTriangEig( const ElementalMatrix<F>& U ); \
  template void schur::QuasiTriangEig \
  ( const ElementalMatrix<F>& U, \
          ElementalMatrix<Complex<Base<F>>>& w );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void schur::RealToComplex \
  ( const Matrix<Real>& UQuasi, \
          Matrix<Complex<Real>>& U ); \
  template void schur::RealToComplex \
  ( const Matrix<Real>& UQuasi, \
    const Matrix<Real>& QQuasi, \
          Matrix<Complex<Real>>& U, \
          Matrix<Complex<Real>>& Q ); \
  template void schur::RealToComplex \
  ( const ElementalMatrix<Real>& UQuasi, \
          ElementalMatrix<Complex<Real>>& U ); \
  template void schur::RealToComplex \
  ( const ElementalMatrix<Real>& UQuasi, \
    const ElementalMatrix<Real>& QQuasi, \
          ElementalMatrix<Complex<Real>>& U, \
          ElementalMatrix<Complex<Real>>& Q );

#define EL_NO_INT_PROTO
#include <El/macros/Instantiate.h>

} // namespace El
