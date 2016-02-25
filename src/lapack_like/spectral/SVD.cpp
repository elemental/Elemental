/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

#include "./SVD/Chan.hpp"
#include "./SVD/Product.hpp"

namespace El {

// Grab the SVD of the general matrix A, A = U diag(s) V^H
// =======================================================

template<typename F>
void SVD
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD [const Matrix Decomp]"))
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
    DEBUG_ONLY(CSE cse("SVD [Matrix Decomp]"))
    typedef Base<F> Real;
    if( !ctrl.overwrite && ctrl.approach != PRODUCT_SVD )
    {
        auto ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        SVD( ACopy, U, s, V, ctrlMod );
        return;
    }
    if( ctrl.avoidComputingU && ctrl.avoidComputingV )
    {
        SVD( A, s, ctrl );
        return;
    }

    if( ctrl.approach == PRODUCT_SVD )
    {
        // TODO: switch to control structure
        svd::Product
        ( A, U, s, V,
          ctrl.tol, ctrl.relative,
          ctrl.avoidComputingU, ctrl.avoidComputingV );
    }
    else if( ctrl.approach == THIN_SVD ||
             ctrl.approach == FULL_SVD ||
             ctrl.approach == COMPACT_SVD )
    {
        const Int m = A.Height();
        const Int n = A.Width();
        const Int k = Min(m,n);
        const bool thin = ( ctrl.approach == THIN_SVD );
        const bool compact = ( ctrl.approach == COMPACT_SVD );
        const bool avoidU = ctrl.avoidComputingU;
        const bool avoidV = ctrl.avoidComputingV;
        s.Resize( k, 1 );
        Matrix<F> VAdj;

        if( ctrl.seqQR )
        {
            if( thin || compact )
            {
                if( !avoidU ) U.Resize( m, k );
                if( !avoidV ) VAdj.Resize( k, n );
            }
            else
            {
                if( !avoidU ) U.Resize( m, m );
                if( !avoidV ) VAdj.Resize( n, n );
            }
            lapack::QRSVD
            ( m, n,
              A.Buffer(), A.LDim(),
              s.Buffer(),
              U.Buffer(), U.LDim(),
              VAdj.Buffer(), VAdj.LDim(),
              (thin||compact), avoidU, avoidV );
        }
        else
        {
            if( thin || compact )
            {
                U.Resize( m, k );
                VAdj.Resize( k, n );
            }
            else
            {
                U.Resize( m, m );
                VAdj.Resize( n, n );
            }
            lapack::DivideAndConquerSVD
            ( m, n,
              A.Buffer(), A.LDim(),
              s.Buffer(),
              U.Buffer(), U.LDim(),
              VAdj.Buffer(), VAdj.LDim(),
              (thin||compact) );
        }

        if( compact )
        {
            const Real twoNorm = ( k==0 ? Real(0) : s.Get(0,0) );
            // Use Max(m,n)*twoNorm*eps unless a manual tolerance is specified
            Real thresh = Max(m,n)*twoNorm*limits::Epsilon<Real>();
            if( ctrl.tol != Real(0) )
            {
                if( ctrl.relative )
                    thresh = twoNorm*ctrl.tol;
                else
                    thresh = ctrl.tol;
            }
            Int rank = k;
            for( Int j=0; j<k; ++j )
            {
                if( s.Get(j,0) <= thresh )
                {
                    rank = j;
                    break;
                }
            }
            s.Resize( rank, 1 );
            if( !avoidU ) U.Resize( m, rank );
            if( !avoidV ) VAdj.Resize( rank, n );
        }
        if( !avoidV ) Adjoint( VAdj, V );
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
    DEBUG_ONLY(CSE cse("SVD [const ElementalMatrix Decomp]"))
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
    DEBUG_ONLY(CSE cse("SVD [ElementalMatrix Decomp]"))
    if( !ctrl.overwrite && ctrl.approach != PRODUCT_SVD )
    {
        DistMatrix<F> ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        SVD( ACopy, U, s, V, ctrlMod );
        return;
    }
    if( ctrl.avoidComputingU && ctrl.avoidComputingV )
    {
        SVD( A, s, ctrl );
        return;
    }

    if( ctrl.approach == PRODUCT_SVD )
    {
        // TODO: Switch to using control structure
        if( U.ColDist() == VC && U.RowDist() == STAR )
        {
            auto& UCast = static_cast<DistMatrix<F,VC,STAR>&>( U );
            svd::Product
            ( A, UCast, s, V,
              ctrl.tol, ctrl.relative,
              ctrl.avoidComputingU, ctrl.avoidComputingV );
        }
        else
            svd::Product
            ( A, U, s, V,
              ctrl.tol, ctrl.relative,
              ctrl.avoidComputingU, ctrl.avoidComputingV );
    }
    else
    {
        svd::Chan( A, U, s, V, ctrl );
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
    DEBUG_ONLY(CSE cse("SVD [const Matrix values]"))
    if( ctrl.approach == PRODUCT_SVD )
    {
        svd::Product( A, s, ctrl.tol, ctrl.relative );
    }
    else
    {
        auto ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        SVD( ACopy, s, ctrlMod );
    }
}

template<typename F>
void SVD
(       Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD [Matrix values]"))
    typedef Base<F> Real;

    Matrix<F> AMod;
    if( ctrl.overwrite )
        View( AMod, A );
    else
        AMod = A;

    if( ctrl.approach == THIN_SVD ||
        ctrl.approach == COMPACT_SVD ||
        ctrl.approach == FULL_SVD )
    {
        const Int m = AMod.Height();
        const Int n = AMod.Width();
        s.Resize( Min(m,n), 1 );
        lapack::SVD( m, n, AMod.Buffer(), AMod.LDim(), s.Buffer() );

        if( ctrl.approach == COMPACT_SVD )
        {
            const Real twoNorm = ( Min(m,n)==0 ? Real(0) : s.Get(0,0) );
            // Use Max(m,n)*twoNorm*eps unless a manual tolerance is specified
            Real thresh = Max(m,n)*twoNorm*limits::Epsilon<Real>();
            if( ctrl.tol != Real(0) )
            {
                if( ctrl.relative )
                    thresh = twoNorm*ctrl.tol;
                else
                    thresh = ctrl.tol;
            }
            Int rank = Min(m,n); 
            for( Int j=0; j<Min(m,n); ++j )
            {
                if( s.Get(j,0) <= thresh )
                {
                    rank = j;
                    break;
                }
            }
            s.Resize( rank, 1 );
        }
    }
    else
    {
        svd::Product( A, s, ctrl.tol, ctrl.relative );
    }
}

template<typename F>
void SVD
( const ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD [const ElementalMatrix values]"))

    if( ctrl.approach == THIN_SVD ||
        ctrl.approach == COMPACT_SVD ||
        ctrl.approach == FULL_SVD )
    {
        DistMatrix<F> ACopy( A );
        svd::Chan( ACopy, s, ctrl );
    }
    else
    {
        svd::Product( A, s, ctrl.tol, ctrl.relative );
    }
}

template<typename F>
void SVD
(       ElementalMatrix<F>& A,
        ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD [ElementalMatrix values]"))
    if( ctrl.approach == PRODUCT_SVD )
    {
        svd::Product( A, s, ctrl.tol, ctrl.relative );
        return;
    }

    if( !ctrl.overwrite )
    {
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        DistMatrix<F> ACopy( A ); 
        SVD( ACopy, s, ctrlMod );
        return;
    }
    svd::Chan( A, s, ctrl );
}

template<typename F>
void SVD
( const DistMatrix<F,MC,MR,BLOCK>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("SVD [const BlockMatrix values]"))
    if( ctrl.approach == PRODUCT_SVD )
        LogicError("Block product SVD not yet supported");

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
    DEBUG_ONLY(CSE cse("SVD [BlockMatrix values]"))
    typedef Base<F> Real;
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

    if( ctrl.approach == THIN_SVD ||
        ctrl.approach == COMPACT_SVD ||
        ctrl.approach == FULL_SVD )
    {
        const int bHandle = blacs::Handle( AMod );
        const int context = blacs::GridInit( bHandle, AMod );
        auto descAMod = FillDesc( AMod, context );

        const bool compact = ( ctrl.approach == COMPACT_SVD );
        s.Resize( k, 1 );
        scalapack::SingularValues
        ( m, n, AMod.Buffer(), descAMod.data(), s.Buffer() ); 
        if( compact )
        {
            const Real twoNorm = ( k==0 ? Real(0) : s.Get(0,0) );
            // Use Max(m,n)*twoNorm*eps unless a manual tolerance is specified
            Real thresh = Max(m,n)*twoNorm*limits::Epsilon<Real>();
            if( ctrl.tol != Real(0) )
            {
                if( ctrl.relative )
                    thresh = twoNorm*ctrl.tol;
                else
                    thresh = ctrl.tol;
            }
            Int rank = k;
            for( Int j=0; j<k; ++j )
            {
                if( s.Get(j,0) <= thresh )
                {
                    rank = j;
                    break;
                }
            }
            s.Resize( rank, 1 );
        }

        // TODO: Cache context, handle, and exit BLACS during El::Finalize()
        blacs::FreeGrid( context );
        blacs::FreeHandle( bHandle );
    }
    else
        LogicError("Block product SVD not yet supported");
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
    DEBUG_ONLY(CSE cse("SVD [const BlockMatrix Decomp]"))
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
    DEBUG_ONLY(CSE cse("SVD [BlockMatrix Decomp]"))
    typedef Base<F> Real;
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

    if( ctrl.approach == THIN_SVD || ctrl.approach == COMPACT_SVD )
    {
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

        const bool compact = ( ctrl.approach == COMPACT_SVD );
        if( compact )
        {
            const Real twoNorm = ( k==0 ? Real(0) : s.Get(0,0) );
            // Use Max(m,n)*twoNorm*eps unless a manual tolerance is specified
            Real thresh = Max(m,n)*twoNorm*limits::Epsilon<Real>();
            if( ctrl.tol != Real(0) )
            {
                if( ctrl.relative )
                    thresh = twoNorm*ctrl.tol;
                else
                    thresh = ctrl.tol;
            }
            Int rank = k;
            for( Int j=0; j<k; ++j )
            {
                if( s.Get(j,0) <= thresh )
                {
                    rank = j;
                    break;
                }
            }
            s.Resize( rank, 1 );
            U.Resize( m, rank );
            VH.Resize( rank, n );
        }

        // TODO: Cache context, handle, and exit BLACS during El::Finalize()
        blacs::FreeGrid( context );
        blacs::FreeHandle( bHandle );

        Adjoint( VH, V );
    }
    else
        LogicError
        ("Only Thin and Compact singular value options currently supported");
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
        const Int mRoot = rootQR.Height();
        const Int nRoot = rootQR.Width();
        const Int kRoot = Min(m,n);

        Matrix<F> URoot, VAdjRoot;
        URoot.Resize( mRoot, kRoot );
        VAdjRoot.Resize( kRoot, nRoot );
        lapack::QRSVD
        ( mRoot, nRoot,
          rootQR.Buffer(), rootQR.LDim(),
          s.Buffer(),
          URoot.Buffer(), URoot.LDim(),
          VAdjRoot.Buffer(), VAdjRoot.LDim() );

        rootQR = URoot; 
        Adjoint( VAdjRoot, V.Matrix() );
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
