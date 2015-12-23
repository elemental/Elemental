/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./SVD/Chan.hpp"
#include "./SVD/Thresholded.hpp"

namespace El {

// Grab the full SVD of the general matrix A, A = U diag(s) V^H
// ============================================================
// NOTE: On exit, A is overwritten with U

template<typename F>
void SVD
( Matrix<F>& A,
  Matrix<Base<F>>& s,
  Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    if( ctrl.thresholded )
    {
        svd::Thresholded( A, s, V, ctrl.tol, ctrl.relative );
    }
    else
    {
        if( ctrl.seqQR )
            svd::QRSVD( A, s, V );
        else
            svd::DivideAndConquerSVD( A, s, V );
    }
}

template<typename F>
void SVD
( ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    if( ctrl.thresholded )
    {
        if( A.ColDist() == VC && A.RowDist() == STAR )
        {
            auto& ACast = static_cast<DistMatrix<F,VC,STAR>&>( A );
            svd::Thresholded( ACast, s, V, ctrl.tol, ctrl.relative );
        }
        else
            svd::Thresholded( A, s, V, ctrl.tol, ctrl.relative );
    }
    else
        svd::Chan( A, s, V, ctrl.fullChanRatio );
}

// Return the singular values
// ==========================

template<typename F>
void SVD( Matrix<F>& A, Matrix<Base<F>>& s )
{
    DEBUG_ONLY(CSE cse("SVD"))
    const Int m = A.Height();
    const Int n = A.Width();
    s.Resize( Min(m,n), 1 );
    lapack::SVD( m, n, A.Buffer(), A.LDim(), s.Buffer() );
}

template<typename F>
void SVD
( ElementalMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    // TODO: Add more options
    svd::Chan( A, s, ctrl.valChanRatio );
}

template<typename F>
void SVD( DistMatrix<F,MC,MR,BLOCK>& A, Matrix<Base<F>>& s )
{
    DEBUG_ONLY(CSE cse("SVD"))
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const int m = A.Height();
    const int n = A.Width();
    const int k = Min(m,n);

    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );

    s.Resize( k, 1 );
    scalapack::SingularValues( m, n, A.Buffer(), descA.data(), s.Buffer() ); 

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
}

template<typename F>
void SVD
( DistMatrix<F,MC,MR,BLOCK>& A,
  Matrix<Base<F>>& s,
  DistMatrix<F,MC,MR,BLOCK>& U,
  DistMatrix<F,MC,MR,BLOCK>& VH )
{
    DEBUG_ONLY(CSE cse("SVD"))
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    const int m = A.Height();
    const int n = A.Width();
    const int k = Min(m,n);
    Zeros( U, m, k );
    Zeros( VH, k, n );

    const int bHandle = blacs::Handle( A );
    const int context = blacs::GridInit( bHandle, A );
    auto descA = FillDesc( A, context );
    auto descU = FillDesc( U, context );
    auto descVH = FillDesc( VH, context );

    s.Resize( k, 1 );
    scalapack::SVD
    ( m, n,
      A.Buffer(), descA.data(),
      s.Buffer(),
      U.Buffer(), descU.data(),
      VH.Buffer(), descVH.data() ); 

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
}

namespace svd {

template<typename F>
void TSQR
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& sPre )
{
    DEBUG_ONLY(CSE cse("svd::TSQR"))

    DistMatrixReadWriteProxy<F,F,VC,STAR> AProx( APre );
    DistMatrixWriteProxy<Base<F>,Base<F>,CIRC,CIRC> sProx( sPre );
    auto& A = AProx.Get();
    auto& s = sProx.Get();
    const Int m = A.Height();
    const Int n = A.Width();
    if( m < n )
        LogicError("svd::TSQR assumes m >= n");
    s.Resize( n, n );

    const Int p = mpi::Size( A.ColComm() );
    if( p == 1 )
        SVD( A.Matrix(), s.Matrix() );

    qr::TreeData<F> treeData;
    treeData.QR0 = A.LockedMatrix();
    QR( treeData.QR0, treeData.t0, treeData.d0 );
    qr::ts::Reduce( A, treeData );
    if( A.ColRank() == 0 )
        SVD( qr::ts::RootQR(A,treeData), s.Matrix() );
    qr::ts::Scatter( A, treeData );
}

template<typename F>
void TSQR
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& sPre,
  ElementalMatrix<F>& VPre )
{
    DEBUG_ONLY(CSE cse("svd::TSQR"))

    DistMatrixReadWriteProxy<F,F,VC,STAR> AProx( APre );
    DistMatrixWriteProxy<Base<F>,Base<F>,CIRC,CIRC> sProx( sPre );
    DistMatrixWriteProxy<F,F,CIRC,CIRC> VProx( VPre );
    auto& A = AProx.Get();
    auto& s = sProx.Get();
    auto& V = VProx.Get();
    const Int m = A.Height();
    const Int n = A.Width();
    if( m < n )
        LogicError("svd::TSQR assumes m >= n");
    s.Resize( n, n );
    V.Resize( n, n );

    const Int p = mpi::Size( A.ColComm() );
    if( p == 1 )
        SVD( A.Matrix(), s.Matrix(), V.Matrix() );

    qr::TreeData<F> treeData;
    treeData.QR0 = A.LockedMatrix();
    QR( treeData.QR0, treeData.t0, treeData.d0 );
    qr::ts::Reduce( A, treeData );
    if( A.ColRank() == 0 )
        SVD( qr::ts::RootQR(A,treeData), s.Matrix(), V.Matrix() );
    qr::ts::Scatter( A, treeData );
}

} // namespace svd

#define PROTO(F) \
  template void SVD \
  ( Matrix<F>& A, \
    Matrix<Base<F>>& s ); \
  template void SVD \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( DistMatrix<F,MC,MR,BLOCK>& A, \
    Matrix<Base<F>>& s ); \
  template void SVD \
  ( Matrix<F>& A, \
    Matrix<Base<F>>& s, \
    Matrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& s, \
    ElementalMatrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( DistMatrix<F,MC,MR,BLOCK>& A, \
    Matrix<Base<F>>& s, \
    DistMatrix<F,MC,MR,BLOCK>& U, \
    DistMatrix<F,MC,MR,BLOCK>& VH ); \
  template void svd::TSQR \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& s ); \
  template void svd::TSQR \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& s, \
    ElementalMatrix<F>& V );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
