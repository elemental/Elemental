/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./Schur/CheckReal.hpp"
#include "./Schur/RealToComplex.hpp"
#include "./Schur/QuasiTriangEig.hpp"
#include "./Schur/QR.hpp"
#include "./Schur/SDC.hpp"
#include "./Schur/InverseFreeSDC.hpp"

namespace El {

template<typename F>
void Schur
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
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
        schur::QR( A, w, fullTriangle );
}

template<typename F>
void Schur
( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, Matrix<F>& Q, bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
    if( ctrl.useSDC )
        schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
    else
        schur::QR( A, w, Q, fullTriangle );
}

template<typename F>
void Schur
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, 
  bool fullTriangle, const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
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
        schur::QR( A, w, fullTriangle, ctrl.qrCtrl );
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
( BlockDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, 
  bool fullTriangle, const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef EL_HAVE_SCALAPACK
    schur::QR( A, w, fullTriangle, ctrl.qrCtrl );
#else
    LogicError("Block distributed Schur currently requires ScaLAPACK support");
#endif
}

template<typename F>
void Schur
( AbstractDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, 
  AbstractDistMatrix<F>& Q, bool fullTriangle, const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef EL_HAVE_SCALAPACK
    if( ctrl.useSDC )
        schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
    else
        schur::QR( A, w, Q, fullTriangle, ctrl.qrCtrl );
#else
    schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
#endif
}

template<typename F>
void Schur
( BlockDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, 
  BlockDistMatrix<F>& Q, bool fullTriangle, const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef EL_HAVE_SCALAPACK
    schur::QR( A, w, Q, fullTriangle, ctrl.qrCtrl );
#else
    LogicError("Block distributed Schur currently requires ScaLAPACK support");
#endif
}

#define PROTO(F) \
  template void Schur \
  ( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, \
    bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, \
    bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( BlockDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, \
    bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, \
    Matrix<F>& Q, bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( AbstractDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, \
    AbstractDistMatrix<F>& Q, bool fullTriangle, \
    const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( BlockDistMatrix<F>& A, AbstractDistMatrix<Complex<Base<F>>>& w, \
    BlockDistMatrix<F>& Q, bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void schur::CheckRealSchur \
  ( const Matrix<F>& U, bool standardForm ); \
  template void schur::CheckRealSchur \
  ( const AbstractDistMatrix<F>& U, bool standardForm ); \
  template void schur::QuasiTriangEig \
  ( const Matrix<F>& dMain, const Matrix<F>& dSub, const Matrix<F>& dSup, \
    Matrix<Complex<Base<F>>>& w ); \
  template void schur::QuasiTriangEig \
  ( const Matrix<F>& U, Matrix<Complex<Base<F>>>& w ); \
  template Matrix<Complex<Base<F>>> schur::QuasiTriangEig \
  ( const Matrix<F>& U ); \
  template DistMatrix<Complex<Base<F>>,VR,STAR> \
  schur::QuasiTriangEig( const AbstractDistMatrix<F>& U ); \
  template void schur::QuasiTriangEig \
  ( const AbstractDistMatrix<F>& U, AbstractDistMatrix<Complex<Base<F>>>& w );

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void schur::RealToComplex \
  ( const Matrix<Real>& UQuasi, Matrix<Complex<Real>>& U ); \
  template void schur::RealToComplex \
  ( const AbstractDistMatrix<Real>& UQuasi, \
          AbstractDistMatrix<Complex<Real>>& U );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
