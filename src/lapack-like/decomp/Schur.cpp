/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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
    if( ctrl.useSdc )
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
    if( ctrl.useSdc )
        schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
    else
        schur::QR( A, w, Q, fullTriangle );
}

template<typename F>
void Schur
( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, bool fullTriangle,
  const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef EL_HAVE_SCALAPACK
    if( ctrl.useSdc )
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
( BlockDistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, 
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
( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, DistMatrix<F>& Q,
  bool fullTriangle, const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef EL_HAVE_SCALAPACK
    if( ctrl.useSdc )
        schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
    else
        schur::QR( A, w, Q, fullTriangle, ctrl.qrCtrl );
#else
    schur::SDC( A, w, Q, fullTriangle, ctrl.sdcCtrl );
#endif
}

template<typename F>
void Schur
( BlockDistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, 
  BlockDistMatrix<F>& Q, bool fullTriangle, const SchurCtrl<Base<F>> ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("Schur"))
#ifdef EL_HAVE_SCALAPACK
    schur::QR( A, w, Q, fullTriangle, ctrl.qrCtrl );
#else
    LogicError("Block distributed Schur currently requires ScaLAPACK support");
#endif
}

#define PROTO_DIST(F,CDIST,RDIST) \
  template void schur::QuasiTriangEig \
  ( const DistMatrix<F>& U, DistMatrix<Complex<Base<F>>,CDIST,RDIST>& w );

#define PROTO(F) \
  template void Schur \
  ( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, \
    bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, \
    bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( BlockDistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, \
    bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( Matrix<F>& A, Matrix<Complex<Base<F>>>& w, \
    Matrix<F>& Q, bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( DistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, \
    DistMatrix<F>& Q, bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void Schur \
  ( BlockDistMatrix<F>& A, DistMatrix<Complex<Base<F>>,VR,STAR>& w, \
    BlockDistMatrix<F>& Q, bool fullTriangle, const SchurCtrl<Base<F>> ctrl ); \
  template void schur::CheckRealSchur \
  ( const Matrix<F>& U, bool standardForm ); \
  template void schur::CheckRealSchur \
  ( const DistMatrix<F>& U, bool standardForm ); \
  template void schur::QuasiTriangEig \
  ( const Matrix<F>& dMain, const Matrix<F>& dSub, const Matrix<F>& dSup, \
    Matrix<Complex<Base<F>>>& w ); \
  template void schur::QuasiTriangEig \
  ( const Matrix<F>& U, Matrix<Complex<Base<F>>>& w ); \
  template Matrix<Complex<Base<F>>> schur::QuasiTriangEig \
  ( const Matrix<F>& U ); \
  template DistMatrix<Complex<Base<F>>,VR,STAR> schur::QuasiTriangEig \
  ( const DistMatrix<F>& U ); \
  PROTO_DIST(F,CIRC,CIRC) \
  PROTO_DIST(F,MC,  MR  ) \
  PROTO_DIST(F,MC,  STAR) \
  PROTO_DIST(F,MD,  STAR) \
  PROTO_DIST(F,MR,  MC  ) \
  PROTO_DIST(F,MR,  STAR) \
  PROTO_DIST(F,STAR,MC  ) \
  PROTO_DIST(F,STAR,MD  ) \
  PROTO_DIST(F,STAR,MR  ) \
  PROTO_DIST(F,STAR,STAR) \
  PROTO_DIST(F,STAR,VC  ) \
  PROTO_DIST(F,STAR,VR  ) \
  PROTO_DIST(F,VC,  STAR) \
  PROTO_DIST(F,VR,  STAR)

#define PROTO_REAL(Real) \
  PROTO(Real) \
  template void schur::RealToComplex \
  ( const Matrix<Real>& UQuasi, Matrix<Complex<Real>>& U ); \
  template void schur::RealToComplex \
  ( const DistMatrix<Real>& UQuasi, DistMatrix<Complex<Real>>& U );

PROTO_REAL(float)
PROTO_REAL(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
