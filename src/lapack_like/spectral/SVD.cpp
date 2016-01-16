/*
   Copyright (c) 2009-2016, Jack Poulson
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

template<typename F>
void SVD
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    auto ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    SVD( ACopy, U, s, V, ctrlMod );
}

template<typename F>
void SVD
(     Matrix<F>& A,
      Matrix<F>& U,
      Matrix<Base<F>>& s,
      Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    if( !ctrl.overwrite )
    {
        auto ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        SVD( ACopy, U, s, V, ctrlMod );
        return;
    }

    if( ctrl.approach == THRESHOLDED_SVD )
    {
        U = A;
        svd::Thresholded( U, s, V, ctrl.tol, ctrl.relative );
    }
    else if( ctrl.approach == THIN_SVD )
    {
        U = A;
        if( ctrl.seqQR )
            svd::QRSVD( U, s, V );
        else
            svd::DivideAndConquerSVD( U, s, V );
    }
    else if( ctrl.approach == COMPACT_SVD )
    {
        // TODO: Extend THIN_SVD approach with a boolean?
        LogicError("This option not yet supported");
    }
    else // ctrl.approach == FULL_SVD
    {
        LogicError("This option not yet supported");
    }
}

template<typename F>
void SVD
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    DistMatrix<F> ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    SVD( ACopy, U, s, V, ctrlMod );
}

template<typename F>
void SVD
(       ElementalMatrix<F>& A,
        ElementalMatrix<F>& U,
        ElementalMatrix<Base<F>>& s, 
        ElementalMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    if( !ctrl.overwrite )
    {
        DistMatrix<F> ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        SVD( ACopy, U, s, V, ctrlMod );
        return;
    }

    if( ctrl.approach == THRESHOLDED_SVD )
    {
        Copy( A, U );
        if( U.ColDist() == VC && U.RowDist() == STAR )
        {
            auto& UCast = static_cast<DistMatrix<F,VC,STAR>&>( U );
            svd::Thresholded( UCast, s, V, ctrl.tol, ctrl.relative );
        }
        else
            svd::Thresholded( U, s, V, ctrl.tol, ctrl.relative );
    }
    else if( ctrl.approach == THIN_SVD )
    {
        Copy( A, U );
        svd::Chan( U, s, V, ctrl.fullChanRatio );
    }
    else if( ctrl.approach == COMPACT_SVD )
    {
        // TODO
        LogicError("This option is not yet supported");
    }
    else // ctrl.approach == FULL_SVD
    {
        // TODO
        LogicError("This option is not yet supported");
    }
}

// Return the singular values
// ==========================

template<typename F>
void SVD
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    auto ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    SVD( ACopy, s, ctrlMod );
}

template<typename F>
void SVD
(       Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    if( ctrl.approach != THIN_SVD )
        LogicError("Only the THIN_SVD option is currently supported");

    Matrix<F> AMod;
    if( ctrl.overwrite )
        View( AMod, A );
    else
        AMod = A;

    const Int m = AMod.Height();
    const Int n = AMod.Width();
    s.Resize( Min(m,n), 1 );
    lapack::SVD( m, n, AMod.Buffer(), AMod.LDim(), s.Buffer() );
}

template<typename F>
void SVD
( const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    if( ctrl.approach != THIN_SVD )
        LogicError("Only the THIN_SVD option is currently supported");

    DistMatrix<F> ACopy( A );
    svd::Chan( ACopy, s, ctrl.valChanRatio );
}

template<typename F>
void SVD
(       ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    if( ctrl.approach != THIN_SVD )
        LogicError("Only the THIN_SVD option is currently supported");

    svd::Chan( A, s, ctrl.valChanRatio );
}

template<typename F>
void SVD
( const DistMatrix<F,MC,MR,BLOCK>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    DistMatrix<F,MC,MR,BLOCK> ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    SVD( ACopy, s, ctrlMod );
}

template<typename F>
void SVD
(       DistMatrix<F,MC,MR,BLOCK>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    DistMatrix<F,MC,MR,BLOCK> AMod( A.Grid() );
    if( ctrl.overwrite )
        View( AMod, A );
    else
        AMod = A;

    const int m = AMod.Height();
    const int n = AMod.Width();
    const int k = Min(m,n);

    const int bHandle = blacs::Handle( AMod );
    const int context = blacs::GridInit( bHandle, AMod );
    auto descAMod = FillDesc( AMod, context );

    s.Resize( k, 1 );
    scalapack::SingularValues
    ( m, n, AMod.Buffer(), descAMod.data(), s.Buffer() ); 

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );
#endif
}

template<typename F>
void SVD
( const DistMatrix<F,MC,MR,BLOCK>& A,
        DistMatrix<F,MC,MR,BLOCK>& U,
        Matrix<Base<F>>& s,
        DistMatrix<F,MC,MR,BLOCK>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    DistMatrix<F,MC,MR,BLOCK> ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    SVD( ACopy, U, s, V, ctrlMod );
}

template<typename F>
void SVD
(       DistMatrix<F,MC,MR,BLOCK>& A,
        DistMatrix<F,MC,MR,BLOCK>& U,
        Matrix<Base<F>>& s,
        DistMatrix<F,MC,MR,BLOCK>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD"))
    AssertScaLAPACKSupport();
#ifdef EL_HAVE_SCALAPACK
    DistMatrix<F,MC,MR,BLOCK> AMod( A.Grid() );
    if( ctrl.overwrite )
        View( AMod, A );
    else
        AMod = A;

    const int m = AMod.Height();
    const int n = AMod.Width();
    const int k = Min(m,n);
    Zeros( U, m, k );
    DistMatrix<F,MC,MR,BLOCK> VH( AMod.Grid() );
    Zeros( VH, k, n );

    const int bHandle = blacs::Handle( AMod );
    const int context = blacs::GridInit( bHandle, AMod );
    auto descAMod = FillDesc( AMod, context );
    auto descU = FillDesc( U, context );
    auto descVH = FillDesc( VH, context );

    s.Resize( k, 1 );
    scalapack::SVD
    ( m, n,
      AMod.Buffer(), descAMod.data(),
      s.Buffer(),
      U.Buffer(), descU.data(),
      VH.Buffer(), descVH.data() ); 

    // TODO: Cache context, handle, and exit BLACS during El::Finalize()
    blacs::FreeGrid( context );
    blacs::FreeHandle( bHandle );

    Adjoint( VH, V );
#endif
}

namespace svd {

template<typename F>
void TSQR
( const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s )
{
    DEBUG_ONLY(CSE cse("svd::TSQR"))
    DistMatrix<F,VC,STAR> ACopy( A );
    TSQR( ACopy, s, true );
}

template<typename F>
void TSQR
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& sPre,
  bool overwrite )
{
    DEBUG_ONLY(CSE cse("svd::TSQR"))
    if( !overwrite )  
    {
        DistMatrix<F,VC,STAR> A( APre );
        TSQR( A, sPre, true );
        return;
    }

    DistMatrixReadProxy<F,F,VC,STAR> AProx( APre );
    DistMatrixWriteProxy<Base<F>,Base<F>,CIRC,CIRC> sProx( sPre );
    auto& A = AProx.Get();
    auto& s = sProx.Get();
    const Int m = A.Height();
    const Int n = A.Width();
    if( m < n )
        LogicError("svd::TSQR assumes m >= n");
    const Int minDim = Min(m,n);
    s.Resize( minDim, 1 );

    const Int p = mpi::Size( A.ColComm() );
    if( p == 1 )
    {
        SVD( A, s );
        return; 
    }

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
( const ElementalMatrix<F>& A,
        ElementalMatrix<F>& UPre,
        ElementalMatrix<Base<F>>& sPre,
        ElementalMatrix<F>& VPre )
{
    DEBUG_ONLY(CSE cse("svd::TSQR"))

    DistMatrixWriteProxy<F,F,VC,STAR> UProx( UPre );
    DistMatrixWriteProxy<Base<F>,Base<F>,CIRC,CIRC> sProx( sPre );
    DistMatrixWriteProxy<F,F,CIRC,CIRC> VProx( VPre );
    auto& U = UProx.Get();
    auto& s = sProx.Get();
    auto& V = VProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    if( m < n )
        LogicError("svd::TSQR assumes m >= n");
    const Int minDim = Min(m,n);
    s.Resize( minDim, 1 );
    V.Resize( minDim, minDim );

    const Int p = mpi::Size( A.ColComm() );
    if( p == 1 )
    {
        SVD( A, U, s, V );
        return;
    }

    Copy( A, U );
    qr::TreeData<F> treeData;
    treeData.QR0 = U.LockedMatrix();
    QR( treeData.QR0, treeData.t0, treeData.d0 );
    qr::ts::Reduce( U, treeData );
    if( U.ColRank() == 0 )
    {
        Matrix<F>& rootQR = qr::ts::RootQR(U,treeData);
        svd::QRSVD( rootQR, s.Matrix(), V.Matrix() );
    }
    qr::ts::Scatter( U, treeData );
}

} // namespace svd

#define PROTO(F) \
  template void SVD \
  (       Matrix<F>& A, \
          Matrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( const Matrix<F>& A, \
          Matrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  (       ElementalMatrix<F>& A, \
          ElementalMatrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( const ElementalMatrix<F>& A, \
          ElementalMatrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  (       DistMatrix<F,MC,MR,BLOCK>& A, \
          Matrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( const DistMatrix<F,MC,MR,BLOCK>& A, \
          Matrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  (       Matrix<F>& A, \
          Matrix<F>& U, \
          Matrix<Base<F>>& s, \
          Matrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( const Matrix<F>& A, \
          Matrix<F>& U, \
          Matrix<Base<F>>& s, \
          Matrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  (       ElementalMatrix<F>& A, \
          ElementalMatrix<F>& U, \
          ElementalMatrix<Base<F>>& s, \
          ElementalMatrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( const ElementalMatrix<F>& A, \
          ElementalMatrix<F>& U, \
          ElementalMatrix<Base<F>>& s, \
          ElementalMatrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  (       DistMatrix<F,MC,MR,BLOCK>& A, \
          DistMatrix<F,MC,MR,BLOCK>& U, \
          Matrix<Base<F>>& s, \
          DistMatrix<F,MC,MR,BLOCK>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void SVD \
  ( const DistMatrix<F,MC,MR,BLOCK>& A, \
          DistMatrix<F,MC,MR,BLOCK>& U, \
          Matrix<Base<F>>& s, \
          DistMatrix<F,MC,MR,BLOCK>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template void svd::TSQR \
  ( ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& s, \
    bool overwrite=false ); \
  template void svd::TSQR \
  ( const ElementalMatrix<F>& A, \
    ElementalMatrix<Base<F>>& s ); \
  template void svd::TSQR \
  ( const ElementalMatrix<F>& A, \
          ElementalMatrix<F>& U, \
          ElementalMatrix<Base<F>>& s, \
          ElementalMatrix<F>& V );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
