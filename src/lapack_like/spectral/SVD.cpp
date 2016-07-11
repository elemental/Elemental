/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

#include "./SVD/Chan.hpp"
#include "./SVD/Product.hpp"

namespace El {

// Grab the SVD of the general matrix A, A = U diag(s) V^H
// =======================================================

namespace svd {

// The following exists primarily for benchmarking purposes
template<typename F,typename=EnableIf<IsBlasScalar<F>>>
SVDInfo ScaLAPACKHelper
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& sPre,
        AbstractDistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    AssertScaLAPACKSupport();
    SVDInfo info;
#ifdef EL_HAVE_SCALAPACK
    typedef Base<F> Real;
    DistMatrix<F,MC,MR,BLOCK> A( APre );
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx(sPre);
    DistMatrixWriteProxy<F,F,MC,MR,BLOCK> UProx(UPre);
    auto& s = sProx.Get();
    auto& U = UProx.Get();

    const int m = A.Height();
    const int n = A.Width();
    const int k = Min(m,n);

    auto approach = ctrl.bidiagSVDCtrl.approach;
    if( approach == THIN_SVD || approach == COMPACT_SVD )
    {
        Zeros( U, m, k );
        DistMatrix<F,MC,MR,BLOCK> VH( A.Grid() );
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

        const bool compact = ( approach == COMPACT_SVD );
        if( compact )
        {
            const Real twoNorm = ( k==0 ? Real(0) : s.Get(0,0) );
            const Real thresh =
              bidiag_svd::APosterioriThreshold
              ( m, n, twoNorm, ctrl.bidiagSVDCtrl );

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
    return info;
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
SVDInfo ScaLAPACKHelper
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& sPre,
        AbstractDistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    SVDInfo info;
    LogicError("ScaLAPACK does not support this datatype");
    return info;
}

// The following exists primarily for benchmarking purposes
template<typename F,typename=EnableIf<IsBlasScalar<F>>>
SVDInfo LAPACKHelper
(     Matrix<F>& A,
      Matrix<F>& U,
      Matrix<Base<F>>& s,
      Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    if( !ctrl.overwrite )
        LogicError("LAPACKHelper assumes ctrl.overwrite == true");
    auto approach = ctrl.bidiagSVDCtrl.approach;
    if( approach != THIN_SVD &&
        approach != FULL_SVD &&
        approach != COMPACT_SVD )
        LogicError("LAPACKHelper assumes THIN_SVD, FULL_SVD, or COMPACT_SVD");

    SVDInfo info;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min(m,n);
    const bool thin = ( approach == THIN_SVD );
    const bool compact = ( approach == COMPACT_SVD );
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    s.Resize( k, 1 );
    Matrix<F> VAdj;

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

    if( compact )
    {
        const Real twoNorm = ( k==0 ? Real(0) : s(0) );
        const Real thresh =
          bidiag_svd::APosterioriThreshold
          ( m, n, twoNorm, ctrl.bidiagSVDCtrl );

        Int rank = k;
        for( Int j=0; j<k; ++j )
        {
            if( s(j) <= thresh )
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

    return info;
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
SVDInfo LAPACKHelper
(     Matrix<F>& A,
      Matrix<F>& U,
      Matrix<Base<F>>& s,
      Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    SVDInfo info;
    LogicError("LAPACK does not support this datatype");
    return info;
}

} // namespace svd

template<typename F>
SVDInfo SVD
( const Matrix<F>& A,
        Matrix<F>& U,
        Matrix<Base<F>>& s,
        Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    auto ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    return SVD( ACopy, U, s, V, ctrlMod );
}

template<typename F>
SVDInfo SVD
(     Matrix<F>& A,
      Matrix<F>& U,
      Matrix<Base<F>>& s,
      Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    if( !ctrl.overwrite && ctrl.bidiagSVDCtrl.approach != PRODUCT_SVD )
    {
        auto ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        return SVD( ACopy, U, s, V, ctrlMod );
    }
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    if( avoidU && avoidV )
    {
        return SVD( A, s, ctrl );
    }
    if( ctrl.useLAPACK )
    {
        return svd::LAPACKHelper( A, U, s, V, ctrl );
    }

    SVDInfo info;
    auto approach = ctrl.bidiagSVDCtrl.approach;
    if( approach == PRODUCT_SVD )
    {
        auto tolType = ctrl.bidiagSVDCtrl.tolType;
        if( tolType == RELATIVE_TO_SELF_SING_VAL_TOL )
            LogicError("Product SVD's inherently require absolute SVD tol's");
        const bool relative = (tolType == RELATIVE_TO_MAX_SING_VAL_TOL);

        // TODO(poulson): switch to control structure
        info = svd::Product
        ( A, U, s, V,
          ctrl.bidiagSVDCtrl.tol, relative, avoidU, avoidV );
    }
    else if( approach == THIN_SVD ||
             approach == FULL_SVD ||
             approach == COMPACT_SVD )
    {
        info = svd::Chan( A, U, s, V, ctrl );
    }
    return info;
}

template<typename F>
SVDInfo SVD
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    DistMatrix<F> ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    return SVD( ACopy, U, s, V, ctrlMod );
}

template<typename F>
SVDInfo SVD
(       AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& U,
        AbstractDistMatrix<Base<F>>& s, 
        AbstractDistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    if( IsBlasScalar<F>::value && ctrl.useScaLAPACK )
    {
        return svd::ScaLAPACKHelper( A, U, s, V, ctrl );
    }
    auto approach = ctrl.bidiagSVDCtrl.approach;
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    if( !ctrl.overwrite && approach != PRODUCT_SVD )
    {
        DistMatrix<F> ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        return SVD( ACopy, U, s, V, ctrlMod );
    }
    if( avoidU && avoidV )
    {
        return SVD( A, s, ctrl );
    }

    SVDInfo info;
    if( approach == PRODUCT_SVD )
    {
        auto tolType = ctrl.bidiagSVDCtrl.tolType;
        if( tolType == RELATIVE_TO_SELF_SING_VAL_TOL )
            LogicError("Product SVD's inherently require absolute SVD tol's");
        const bool relative = (tolType == RELATIVE_TO_MAX_SING_VAL_TOL);

        // TODO: Switch to using control structure
        if( U.ColDist() == VC && U.RowDist() == STAR )
        {
            auto& UCast = static_cast<DistMatrix<F,VC,STAR>&>( U );
            info = svd::Product
            ( A, UCast, s, V,
              ctrl.bidiagSVDCtrl.tol, relative, avoidU, avoidV );
        }
        else
            info = svd::Product
            ( A, U, s, V,
              ctrl.bidiagSVDCtrl.tol, relative, avoidU, avoidV );
    }
    else
    {
        info = svd::Chan( A, U, s, V, ctrl );
    }
    return info;
}

// Return the singular values
// ==========================

namespace svd {

// The following exists primarily for benchmarking purposes
template<typename F,typename=EnableIf<IsBlasScalar<F>>>
SVDInfo LAPACKHelper
(       Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    if( !ctrl.overwrite )
        LogicError("LAPACKHelper assumes ctrl.overwrite == true");
    if( ctrl.bidiagSVDCtrl.approach != THIN_SVD &&
        ctrl.bidiagSVDCtrl.approach != FULL_SVD &&
        ctrl.bidiagSVDCtrl.approach != COMPACT_SVD )
        LogicError("LAPACKHelper assumes THIN_SVD, FULL_SVD, or COMPACT_SVD");

    SVDInfo info;
    const Int m = A.Height();
    const Int n = A.Width();
    s.Resize( Min(m,n), 1 );
    // This should map down to DQDS...
    lapack::SVD( m, n, A.Buffer(), A.LDim(), s.Buffer() );

    if( ctrl.bidiagSVDCtrl.approach == COMPACT_SVD )
    {
        const Real twoNorm = ( Min(m,n)==0 ? Real(0) : s(0) );
        const Real thresh =
          bidiag_svd::APosterioriThreshold
          ( m, n, twoNorm, ctrl.bidiagSVDCtrl );

        Int rank = Min(m,n); 
        for( Int j=0; j<Min(m,n); ++j )
        {
            if( s(j) <= thresh )
            {
                rank = j;
                break;
            }
        }
        s.Resize( rank, 1 );
    }

    return info;
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
SVDInfo LAPACKHelper
(       Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    SVDInfo info;
    LogicError("LAPACK does not support this datatype");
    return info;
}

// The following exists primarily for benchmarking purposes
template<typename F,typename=EnableIf<IsBlasScalar<F>>>
SVDInfo ScaLAPACKHelper
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<Base<F>>& sPre,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    AssertScaLAPACKSupport();
    SVDInfo info;
#ifdef EL_HAVE_SCALAPACK
    typedef Base<F> Real;
    DistMatrix<F,MC,MR,BLOCK> A( APre );
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx(sPre);
    auto& s = sProx.Get();
    const int m = A.Height();
    const int n = A.Width();
    const int k = Min(m,n);

    if( ctrl.bidiagSVDCtrl.approach == THIN_SVD ||
        ctrl.bidiagSVDCtrl.approach == COMPACT_SVD ||
        ctrl.bidiagSVDCtrl.approach == FULL_SVD )
    {
        const int bHandle = blacs::Handle( A );
        const int context = blacs::GridInit( bHandle, A );
        auto descA = FillDesc( A, context );

        const bool compact = ( ctrl.bidiagSVDCtrl.approach == COMPACT_SVD );
        s.Resize( k, 1 );
        scalapack::SingularValues
        ( m, n, A.Buffer(), descA.data(), s.Buffer() ); 
        if( compact )
        {
            const Real twoNorm = ( k==0 ? Real(0) : s.Get(0,0) );
            const Real thresh =
              bidiag_svd::APosterioriThreshold
              ( m, n, twoNorm, ctrl.bidiagSVDCtrl );

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
    return info;
}

template<typename F,typename=DisableIf<IsBlasScalar<F>>,typename=void>
SVDInfo ScaLAPACKHelper
( const AbstractDistMatrix<F>& APre,
        AbstractDistMatrix<Base<F>>& sPre,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    SVDInfo info;
    LogicError("ScaLAPACK does not support this datatype");
    return info;
}

} // namespace svd

template<typename F>
SVDInfo SVD
( const Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    if( ctrl.bidiagSVDCtrl.approach == PRODUCT_SVD )
    {
        auto tolType = ctrl.bidiagSVDCtrl.tolType;
        if( tolType == RELATIVE_TO_SELF_SING_VAL_TOL )
            LogicError("Product SVD's inherently require absolute SVD tol's");
        const bool relative = (tolType == RELATIVE_TO_MAX_SING_VAL_TOL);
        return svd::Product( A, s, ctrl.bidiagSVDCtrl.tol, relative );
    }
    else
    {
        auto ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        return SVD( ACopy, s, ctrlMod );
    }
}

template<typename F>
SVDInfo SVD
(       Matrix<F>& A,
        Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;

    Matrix<F> AMod;
    if( ctrl.overwrite )
        View( AMod, A );
    else
        AMod = A;

    if( IsBlasScalar<F>::value && ctrl.useLAPACK )
    {
        return svd::LAPACKHelper( A, s, ctrl );
    }

    SVDInfo info;
    if( ctrl.bidiagSVDCtrl.approach == THIN_SVD ||
        ctrl.bidiagSVDCtrl.approach == COMPACT_SVD ||
        ctrl.bidiagSVDCtrl.approach == FULL_SVD )
    {
        return svd::Chan( A, s, ctrl );
    }
    else
    {
        auto tolType = ctrl.bidiagSVDCtrl.tolType;
        if( tolType == RELATIVE_TO_SELF_SING_VAL_TOL )
            LogicError("Product SVD's inherently require absolute SVD tol's");
        const bool relative = (tolType == RELATIVE_TO_MAX_SING_VAL_TOL);

        info = svd::Product( A, s, ctrl.bidiagSVDCtrl.tol, relative );
    }
    return info;
}

template<typename F>
SVDInfo SVD
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    if( IsBlasScalar<F>::value && ctrl.useScaLAPACK )
    {
        return svd::ScaLAPACKHelper( A, s, ctrl );
    }
    if( ctrl.bidiagSVDCtrl.approach == THIN_SVD ||
        ctrl.bidiagSVDCtrl.approach == COMPACT_SVD ||
        ctrl.bidiagSVDCtrl.approach == FULL_SVD )
    {
        DistMatrix<F> ACopy( A );
        return svd::Chan( ACopy, s, ctrl );
    }
    else
    {
        auto tolType = ctrl.bidiagSVDCtrl.tolType;
        if( tolType == RELATIVE_TO_SELF_SING_VAL_TOL )
            LogicError("Product SVD's inherently require absolute SVD tol's");
        const bool relative = (tolType == RELATIVE_TO_MAX_SING_VAL_TOL);

        return svd::Product( A, s, ctrl.bidiagSVDCtrl.tol, relative );
    }
}

template<typename F>
SVDInfo SVD
(       AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    if( IsBlasScalar<F>::value && ctrl.useScaLAPACK )
    {
        return svd::ScaLAPACKHelper( A, s, ctrl );
    }
    if( ctrl.bidiagSVDCtrl.approach == PRODUCT_SVD )
    {
        auto tolType = ctrl.bidiagSVDCtrl.tolType;
        if( tolType == RELATIVE_TO_SELF_SING_VAL_TOL )
            LogicError("Product SVD's inherently require absolute SVD tol's");
        const bool relative = (tolType == RELATIVE_TO_MAX_SING_VAL_TOL);
        return svd::Product( A, s, ctrl.bidiagSVDCtrl.tol, relative );
    }

    if( !ctrl.overwrite )
    {
        auto ctrlMod( ctrl );
        ctrlMod.overwrite = true;
        DistMatrix<F> ACopy( A ); 
        return SVD( ACopy, s, ctrlMod );
    }
    return svd::Chan( A, s, ctrl );
}

namespace svd {

template<typename F>
SVDInfo TSQR
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<Base<F>>& s )
{
    DEBUG_CSE
    DistMatrix<F,VC,STAR> ACopy( A );
    return TSQR( ACopy, s, true );
}

template<typename F>
SVDInfo TSQR
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& sPre,
  bool overwrite )
{
    DEBUG_CSE
    if( !overwrite )  
    {
        DistMatrix<F,VC,STAR> A( APre );
        return TSQR( A, sPre, true );
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
        return SVD( A, s );
    }

    SVDInfo info;
    qr::TreeData<F> treeData;
    treeData.QR0 = A.LockedMatrix();
    QR( treeData.QR0, treeData.phase0, treeData.signature0 );
    qr::ts::Reduce( A, treeData );
    if( A.ColRank() == 0 )
        info = SVD( qr::ts::RootQR(A,treeData), s.Matrix() );
    // TODO(poulson): Broadcast info from root?
    qr::ts::Scatter( A, treeData );
    return info;
}

template<typename F>
SVDInfo TSQR
( const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<Base<F>>& sPre,
        AbstractDistMatrix<F>& VPre )
{
    DEBUG_CSE

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
        return SVD( A, U, s, V );
    }

    SVDInfo info;
    Copy( A, U );
    qr::TreeData<F> treeData;
    treeData.QR0 = U.LockedMatrix();
    QR( treeData.QR0, treeData.phase0, treeData.signature0 );
    qr::ts::Reduce( U, treeData );
    if( U.ColRank() == 0 )
    {
        Matrix<F>& rootQR = qr::ts::RootQR(U,treeData);
        const Int mRoot = rootQR.Height();
        const Int nRoot = rootQR.Width();
        const Int kRoot = Min(m,n);

        Matrix<F> URoot, VRoot;
        URoot.Resize( mRoot, kRoot );
        VRoot.Resize( nRoot, kRoot );
        SVD( rootQR, URoot, s.Matrix(), VRoot );

        rootQR = URoot; 
        V.Matrix() = VRoot;
    }
    qr::ts::Scatter( U, treeData );
    return info;
}

} // namespace svd

#define PROTO(F) \
  template SVDInfo SVD \
  (       Matrix<F>& A, \
          Matrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template SVDInfo SVD \
  ( const Matrix<F>& A, \
          Matrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template SVDInfo SVD \
  (       AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template SVDInfo SVD \
  ( const AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<Base<F>>& s, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template SVDInfo SVD \
  (       Matrix<F>& A, \
          Matrix<F>& U, \
          Matrix<Base<F>>& s, \
          Matrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template SVDInfo SVD \
  ( const Matrix<F>& A, \
          Matrix<F>& U, \
          Matrix<Base<F>>& s, \
          Matrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template SVDInfo SVD \
  (       AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<F>& U, \
          AbstractDistMatrix<Base<F>>& s, \
          AbstractDistMatrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template SVDInfo SVD \
  ( const AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<F>& U, \
          AbstractDistMatrix<Base<F>>& s, \
          AbstractDistMatrix<F>& V, \
    const SVDCtrl<Base<F>>& ctrl ); \
  template SVDInfo svd::TSQR \
  ( AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Base<F>>& s, \
    bool overwrite=false ); \
  template SVDInfo svd::TSQR \
  ( const AbstractDistMatrix<F>& A, \
    AbstractDistMatrix<Base<F>>& s ); \
  template SVDInfo svd::TSQR \
  ( const AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<F>& U, \
          AbstractDistMatrix<Base<F>>& s, \
          AbstractDistMatrix<F>& V );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
