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
template<typename Field,
         typename=EnableIf<IsBlasScalar<Field>>>
SVDInfo ScaLAPACKHelper
( const AbstractDistMatrix<Field>& APre,
        AbstractDistMatrix<Field>& UPre,
        AbstractDistMatrix<Base<Field>>& sPre,
        AbstractDistMatrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    AssertScaLAPACKSupport();
    SVDInfo info;
#ifdef EL_HAVE_SCALAPACK
    typedef Base<Field> Real;
    DistMatrix<Field,MC,MR,BLOCK> A( APre );
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx(sPre);
    DistMatrixWriteProxy<Field,Field,MC,MR,BLOCK> UProx(UPre);
    auto& s = sProx.Get();
    auto& U = UProx.Get();

    const int m = A.Height();
    const int n = A.Width();
    const int k = Min(m,n);

    auto approach = ctrl.bidiagSVDCtrl.approach;
    if( approach == THIN_SVD || approach == COMPACT_SVD )
    {
        Zeros( U, m, k );
        DistMatrix<Field,MC,MR,BLOCK> VH( A.Grid() );
        Zeros( VH, k, n );
        s.Resize( k, 1 );

        auto descA = FillDesc( A );
        auto descU = FillDesc( U );
        auto descVH = FillDesc( VH );
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

        Adjoint( VH, V );
    }
    else
        LogicError
        ("Only Thin and Compact singular value options currently supported");
#endif
    return info;
}

template<typename Field,
         typename=DisableIf<IsBlasScalar<Field>>,
         typename=void>
SVDInfo ScaLAPACKHelper
( const AbstractDistMatrix<Field>& APre,
        AbstractDistMatrix<Field>& UPre,
        AbstractDistMatrix<Base<Field>>& sPre,
        AbstractDistMatrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SVDInfo info;
    LogicError("ScaLAPACK does not support this datatype");
    return info;
}

// The following exists primarily for benchmarking purposes
template<typename Field,
         typename=EnableIf<IsBlasScalar<Field>>>
SVDInfo LAPACKHelper
(     Matrix<Field>& A,
      Matrix<Field>& U,
      Matrix<Base<Field>>& s,
      Matrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
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
    Matrix<Field> VAdj;

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

template<typename Field,
         typename=DisableIf<IsBlasScalar<Field>>,
         typename=void>
SVDInfo LAPACKHelper
(     Matrix<Field>& A,
      Matrix<Field>& U,
      Matrix<Base<Field>>& s,
      Matrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SVDInfo info;
    LogicError("LAPACK does not support this datatype");
    return info;
}

} // namespace svd

template<typename Field>
SVDInfo SVD
( const Matrix<Field>& A,
        Matrix<Field>& U,
        Matrix<Base<Field>>& s,
        Matrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    auto ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    return SVD( ACopy, U, s, V, ctrlMod );
}

template<typename Field>
SVDInfo SVD
(     Matrix<Field>& A,
      Matrix<Field>& U,
      Matrix<Base<Field>>& s,
      Matrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    const auto& bidiagSVDCtrl = ctrl.bidiagSVDCtrl;
    if( (bidiagSVDCtrl.wantU && bidiagSVDCtrl.accumulateU) ||
        (bidiagSVDCtrl.wantV && bidiagSVDCtrl.accumulateV) )
        LogicError("SVD does not support singular vector accumulation");

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

template<typename Field>
SVDInfo SVD
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrix<Field> ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.overwrite = true;
    return SVD( ACopy, U, s, V, ctrlMod );
}

template<typename Field>
SVDInfo SVD
(       AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& U,
        AbstractDistMatrix<Base<Field>>& s,
        AbstractDistMatrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    const auto& bidiagSVDCtrl = ctrl.bidiagSVDCtrl;
    if( (bidiagSVDCtrl.wantU && bidiagSVDCtrl.accumulateU) ||
        (bidiagSVDCtrl.wantV && bidiagSVDCtrl.accumulateV) )
        LogicError("SVD does not support singular vector accumulation");

    if( IsBlasScalar<Field>::value && ctrl.useScaLAPACK )
    {
        return svd::ScaLAPACKHelper( A, U, s, V, ctrl );
    }
    auto approach = ctrl.bidiagSVDCtrl.approach;
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    if( !ctrl.overwrite && approach != PRODUCT_SVD )
    {
        DistMatrix<Field> ACopy( A );
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
            auto& UCast = static_cast<DistMatrix<Field,VC,STAR>&>( U );
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
template<typename Field,
         typename=EnableIf<IsBlasScalar<Field>>>
SVDInfo LAPACKHelper
(       Matrix<Field>& A,
        Matrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
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

template<typename Field,
         typename=DisableIf<IsBlasScalar<Field>>,
         typename=void>
SVDInfo LAPACKHelper
(       Matrix<Field>& A,
        Matrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SVDInfo info;
    LogicError("LAPACK does not support this datatype");
    return info;
}

// The following exists primarily for benchmarking purposes
template<typename Field,
         typename=EnableIf<IsBlasScalar<Field>>>
SVDInfo ScaLAPACKHelper
( const AbstractDistMatrix<Field>& APre,
        AbstractDistMatrix<Base<Field>>& sPre,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    AssertScaLAPACKSupport();
    SVDInfo info;
#ifdef EL_HAVE_SCALAPACK
    typedef Base<Field> Real;
    DistMatrix<Field,MC,MR,BLOCK> A( APre );
    DistMatrixWriteProxy<Real,Real,STAR,STAR> sProx(sPre);
    auto& s = sProx.Get();
    const int m = A.Height();
    const int n = A.Width();
    const int k = Min(m,n);

    if( ctrl.bidiagSVDCtrl.approach == THIN_SVD ||
        ctrl.bidiagSVDCtrl.approach == COMPACT_SVD ||
        ctrl.bidiagSVDCtrl.approach == FULL_SVD )
    {
        const bool compact = ( ctrl.bidiagSVDCtrl.approach == COMPACT_SVD );
        s.Resize( k, 1 );

        auto descA = FillDesc( A );
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
    }
    else
        LogicError("Block product SVD not yet supported");
#endif
    return info;
}

template<typename Field,
         typename=DisableIf<IsBlasScalar<Field>>,
         typename=void>
SVDInfo ScaLAPACKHelper
( const AbstractDistMatrix<Field>& APre,
        AbstractDistMatrix<Base<Field>>& sPre,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SVDInfo info;
    LogicError("ScaLAPACK does not support this datatype");
    return info;
}

} // namespace svd

template<typename Field>
SVDInfo SVD
( const Matrix<Field>& A,
        Matrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
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

template<typename Field>
SVDInfo SVD
(       Matrix<Field>& A,
        Matrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE

    Matrix<Field> AMod;
    if( ctrl.overwrite )
        View( AMod, A );
    else
        AMod = A;

    if( IsBlasScalar<Field>::value && ctrl.useLAPACK )
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

template<typename Field>
SVDInfo SVD
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    if( IsBlasScalar<Field>::value && ctrl.useScaLAPACK )
    {
        return svd::ScaLAPACKHelper( A, s, ctrl );
    }
    if( ctrl.bidiagSVDCtrl.approach == THIN_SVD ||
        ctrl.bidiagSVDCtrl.approach == COMPACT_SVD ||
        ctrl.bidiagSVDCtrl.approach == FULL_SVD )
    {
        DistMatrix<Field> ACopy( A );
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

template<typename Field>
SVDInfo SVD
(       AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    if( IsBlasScalar<Field>::value && ctrl.useScaLAPACK )
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
        DistMatrix<Field> ACopy( A );
        return SVD( ACopy, s, ctrlMod );
    }
    return svd::Chan( A, s, ctrl );
}

namespace svd {

template<typename Field>
SVDInfo TSQR
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Base<Field>>& s )
{
    EL_DEBUG_CSE
    DistMatrix<Field,VC,STAR> ACopy( A );
    return TSQR( ACopy, s, true );
}

template<typename Field>
SVDInfo TSQR
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Base<Field>>& sPre,
  bool overwrite )
{
    EL_DEBUG_CSE
    if( !overwrite )
    {
        DistMatrix<Field,VC,STAR> A( APre );
        return TSQR( A, sPre, true );
    }

    DistMatrixReadProxy<Field,Field,VC,STAR> AProx( APre );
    DistMatrixWriteProxy<Base<Field>,Base<Field>,CIRC,CIRC> sProx( sPre );
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
    qr::TreeData<Field> treeData;
    treeData.QR0 = A.LockedMatrix();
    QR( treeData.QR0, treeData.householderScalars0, treeData.signature0 );
    qr::ts::Reduce( A, treeData );
    if( A.ColRank() == 0 )
        info = SVD( qr::ts::RootQR(A,treeData), s.Matrix() );
    // TODO(poulson): Broadcast info from root?
    qr::ts::Scatter( A, treeData );
    return info;
}

template<typename Field>
SVDInfo TSQR
( const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& UPre,
        AbstractDistMatrix<Base<Field>>& sPre,
        AbstractDistMatrix<Field>& VPre )
{
    EL_DEBUG_CSE

    DistMatrixWriteProxy<Field,Field,VC,STAR> UProx( UPre );
    DistMatrixWriteProxy<Base<Field>,Base<Field>,CIRC,CIRC> sProx( sPre );
    DistMatrixWriteProxy<Field,Field,CIRC,CIRC> VProx( VPre );
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
    qr::TreeData<Field> treeData;
    treeData.QR0 = U.LockedMatrix();
    QR( treeData.QR0, treeData.householderScalars0, treeData.signature0 );
    qr::ts::Reduce( U, treeData );
    if( U.ColRank() == 0 )
    {
        Matrix<Field>& rootQR = qr::ts::RootQR(U,treeData);
        const Int mRoot = rootQR.Height();
        const Int nRoot = rootQR.Width();
        const Int kRoot = Min(m,n);

        Matrix<Field> URoot, VRoot;
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

#define PROTO(Field) \
  template SVDInfo SVD \
  (       Matrix<Field>& A, \
          Matrix<Base<Field>>& s, \
    const SVDCtrl<Base<Field>>& ctrl ); \
  template SVDInfo SVD \
  ( const Matrix<Field>& A, \
          Matrix<Base<Field>>& s, \
    const SVDCtrl<Base<Field>>& ctrl ); \
  template SVDInfo SVD \
  (       AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Base<Field>>& s, \
    const SVDCtrl<Base<Field>>& ctrl ); \
  template SVDInfo SVD \
  ( const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Base<Field>>& s, \
    const SVDCtrl<Base<Field>>& ctrl ); \
  template SVDInfo SVD \
  (       Matrix<Field>& A, \
          Matrix<Field>& U, \
          Matrix<Base<Field>>& s, \
          Matrix<Field>& V, \
    const SVDCtrl<Base<Field>>& ctrl ); \
  template SVDInfo SVD \
  ( const Matrix<Field>& A, \
          Matrix<Field>& U, \
          Matrix<Base<Field>>& s, \
          Matrix<Field>& V, \
    const SVDCtrl<Base<Field>>& ctrl ); \
  template SVDInfo SVD \
  (       AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& U, \
          AbstractDistMatrix<Base<Field>>& s, \
          AbstractDistMatrix<Field>& V, \
    const SVDCtrl<Base<Field>>& ctrl ); \
  template SVDInfo SVD \
  ( const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& U, \
          AbstractDistMatrix<Base<Field>>& s, \
          AbstractDistMatrix<Field>& V, \
    const SVDCtrl<Base<Field>>& ctrl ); \
  template SVDInfo svd::TSQR \
  ( AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Base<Field>>& s, \
    bool overwrite ); \
  template SVDInfo svd::TSQR \
  ( const AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Base<Field>>& s ); \
  template SVDInfo svd::TSQR \
  ( const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& U, \
          AbstractDistMatrix<Base<Field>>& s, \
          AbstractDistMatrix<Field>& V );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
