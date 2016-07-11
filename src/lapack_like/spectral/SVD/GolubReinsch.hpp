/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVD_GOLUBREINSCH_HPP
#define EL_SVD_GOLUBREINSCH_HPP

#include "./Util.hpp"

namespace El {

namespace bidiag_svd {

template<typename Real>
bidiag_svd::QRInfo
QRAlg
( Matrix<Real>& mainDiag,
  Matrix<Real>& superDiag,
  const BidiagSVDCtrl<Real>& ctrl );

template<typename F>
bidiag_svd::QRInfo
QRAlg
( Matrix<Base<F>>& mainDiag,
  Matrix<Base<F>>& superDiag,
  Matrix<F>& U,
  Matrix<F>& V,
  const BidiagSVDCtrl<Base<F>>& ctrl );

} // namespace bidiag_svd

namespace svd {

template<typename F>
SVDInfo GolubReinsch
( Matrix<F>& A,
  Matrix<F>& U,
  Matrix<Base<F>>& s, 
  Matrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    if( avoidU && avoidV )
    {
        return SVD( A, s, ctrl );
    }
    if( uplo == LOWER )
        LogicError("Lower bidiagonal not yet supported");
    SVDInfo info;

    // Bidiagonalize A
    Timer timer;
    Matrix<F> phaseP, phaseQ;
    if( ctrl.time )
        timer.Start();
    Bidiag( A, phaseP, phaseQ );
    if( ctrl.time )
        Output("Reduction to bidiagonal: ",timer.Stop()," seconds");

    // Grab copies of the diagonal and sub/super-diagonal of A
    // (Force the buffer holding e to be of length at least n)
    // NOTE: lapack::BidiagSVDQRAlg expects e to be of length k
    s = GetRealPartOfDiagonal(A);
    Matrix<Real> ePlus(k,1);
    auto e = ePlus(IR(0,k-1),ALL);
    e = GetRealPartOfDiagonal(A,offdiagonal);

    // TODO: If compact SVD, identify the rank with DQDS first?

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and V
    if( ctrl.time )
        timer.Start();
    bidiag_svd::QRAlg( s, e, U, V, ctrl.bidiagSVDCtrl );
    if( ctrl.time )
        Output("Bidiag SVD: ",timer.Stop()," seconds");

    Int rank = k;
    const bool compact = ( ctrl.bidiagSVDCtrl.approach == COMPACT_SVD );
    if( compact )
    {
        const Real twoNorm = ( k==0 ? Real(0) : s(0) );
        const Real thresh =
          bidiag_svd::APosterioriThreshold
          ( m, n, twoNorm, ctrl.bidiagSVDCtrl );

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
        if( !avoidV ) V.Resize( n, rank );
    }

    // Backtransform U and V
    if( ctrl.time )
        timer.Start();
    if( !avoidU ) bidiag::ApplyQ( LEFT, NORMAL, A, phaseQ, U );
    if( !avoidV ) bidiag::ApplyP( LEFT, NORMAL, A, phaseP, V );
    if( ctrl.time )
        Output("GolubReinsch backtransformation: ",timer.Stop()," seconds");

    return info;
}

template<typename F>
SVDInfo GolubReinsch
( DistMatrix<F>& A,
  DistMatrix<F>& U,
  AbstractDistMatrix<Base<F>>& s, 
  DistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    const Grid& g = A.Grid();
    if( avoidU && avoidV )
    {
        return SVD( A, s, ctrl );
    }
    if( uplo == LOWER )
        LogicError("Lower bidiagonal not yet supported");
    SVDInfo info;

    // Bidiagonalize A
    Timer timer;
    DistMatrix<F,STAR,STAR> phaseP(g), phaseQ(g);
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    Bidiag( A, phaseP, phaseQ );
    if( ctrl.time && g.Rank() == 0 )
        Output("Reduction to bidiagonal: ",timer.Stop()," seconds");

    // Grab copies of the diagonal and sub/super-diagonal of A
    auto d_MD_STAR = GetRealPartOfDiagonal(A);
    auto e_MD_STAR = GetRealPartOfDiagonal(A,offdiagonal);

    // NOTE: lapack::BidiagSVDQRAlg expects e to be of length k
    typedef Base<F> Real;
    DistMatrix<Real,STAR,STAR>
      d_STAR_STAR( d_MD_STAR ),
      eHat_STAR_STAR( k, 1, g );
    auto e_STAR_STAR = eHat_STAR_STAR( IR(0,k-1), ALL );
    e_STAR_STAR = e_MD_STAR;

    // Initialize U and V to the appropriate identity matrices
    auto bidiagSVDCtrlMod = ctrl.bidiagSVDCtrl;
    bidiagSVDCtrlMod.accumulateU = true;
    bidiagSVDCtrlMod.accumulateV = true;
    DistMatrix<F,VC,STAR> U_VC_STAR( g );
    if( bidiagSVDCtrlMod.wantU )
    {
        U_VC_STAR.AlignWith( A );
        Identity( U_VC_STAR, m, k );
    }
    DistMatrix<F,VC,STAR> V_VC_STAR( g );
    if( bidiagSVDCtrlMod.wantV )
    {
        V_VC_STAR.AlignWith( V );
        Identity( V_VC_STAR, n, k );
    }

    // TODO: If compact SVD, identify the rank with DQDS first?

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and V
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    bidiag_svd::QRAlg
    ( d_STAR_STAR.Matrix(), e_STAR_STAR.Matrix(), 
      U_VC_STAR.Matrix(), V_VC_STAR.Matrix(), bidiagSVDCtrlMod );
    if( ctrl.time )
    {
        mpi::Barrier( g.Comm() );
        if( g.Rank() == 0 )
            Output("Bidiag SVD: ",timer.Stop()," seconds");
    }

    Int rank = k;
    const bool compact = ( ctrl.bidiagSVDCtrl.approach == COMPACT_SVD );
    if( compact )
    {
        const Real twoNorm = ( k==0 ? Real(0) : d_STAR_STAR.Get(0,0) );
        const Real thresh =
          bidiag_svd::APosterioriThreshold
          ( m, n, twoNorm, ctrl.bidiagSVDCtrl );

        for( Int j=0; j<k; ++j ) 
        {
            if( d_STAR_STAR.Get(j,0) <= thresh )
            {
                rank = j;
                break;
            }
        }

        d_STAR_STAR.Resize( rank, 1 );
        if( !avoidU ) U_VC_STAR.Resize( m, rank );
        if( !avoidV ) V_VC_STAR.Resize( n, rank );
    }
    // Copy out the appropriate subset of the singular values
    Copy( d_STAR_STAR, s );

    if( m >= n )
    {
        if( !avoidU )
        {
            U.Resize( m, rank );
            auto UT = U( IR(0,n  ), ALL );
            auto UB = U( IR(n,END), ALL );
            auto UT_VC_STAR = U_VC_STAR( IR(0,n), ALL );
            UT = UT_VC_STAR;
            Zero( UB );
        }
        if( !avoidV ) V = V_VC_STAR;
    }
    else
    {
        if( !avoidU ) U = U_VC_STAR;
        if( !avoidV )
        {
            auto VT_VC_STAR = V_VC_STAR( IR(0,m), IR(0,rank) );
            V.Resize( n, rank );
            auto VT = V( IR(0,m  ), ALL );
            auto VB = V( IR(m,END), ALL );
            VT = VT_VC_STAR;
            Zero( VB );
        }
    }

    // Backtransform U and V
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    if( !avoidU ) bidiag::ApplyQ( LEFT, NORMAL, A, phaseQ, U );
    if( !avoidV ) bidiag::ApplyP( LEFT, NORMAL, A, phaseP, V );
    if( ctrl.time && g.Rank() == 0 )
        Output("GolubReinsch backtransformation: ",timer.Stop()," seconds");

    return info;
}

template<typename F>
void GolubReinsch
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& UPre,
  AbstractDistMatrix<Base<F>>& s, 
  AbstractDistMatrix<F>& VPre,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return GolubReinsch( A, U, s, V, ctrl );
}

template<typename F>
SVDInfo GolubReinsch
( Matrix<F>& A,
  Matrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    SVDInfo info;

    // Bidiagonalize A
    Timer timer;
    Matrix<F> phaseP, phaseQ;
    if( ctrl.time )
        timer.Start();
    Bidiag( A, phaseP, phaseQ );
    if( ctrl.time )
        Output("Reduction to bidiagonal: ",timer.Stop()," seconds");

    // Grab copies of the diagonal and sub/super-diagonal of A
    // (Force the buffer holding e to be of length at least n)
    // NOTE: lapack::BidiagDQDS expects e to be of length k
    s = GetRealPartOfDiagonal(A);
    Matrix<Real> ePlus(k,1);
    auto e = ePlus(IR(0,k-1),ALL);
    e = GetRealPartOfDiagonal(A,offdiagonal);

    if( ctrl.time )
        timer.Start();
    bidiag_svd::QRAlg( s, e, ctrl.bidiagSVDCtrl );
    if( ctrl.time )
        Output("Bidiag SVD: ",timer.Stop()," seconds");

    const bool compact = ( ctrl.bidiagSVDCtrl.approach == COMPACT_SVD );
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
    }

    return info;
}

template<typename F>
SVDInfo GolubReinsch
( DistMatrix<F>& A,
  AbstractDistMatrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    const Grid& g = A.Grid();
    SVDInfo info;

    // Bidiagonalize A
    Timer timer;
    DistMatrix<F,STAR,STAR> phaseP(g), phaseQ(g);
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    Bidiag( A, phaseP, phaseQ );
    if( ctrl.time && g.Rank() == 0 )
        Output("Reduction to bidiagonal: ",timer.Stop()," seconds");

    // Grab copies of the diagonal and sub/super-diagonal of A
    auto d_MD_STAR = GetRealPartOfDiagonal(A);
    auto e_MD_STAR = GetRealPartOfDiagonal(A,offdiagonal);

    // In order to use serial DQDS kernels, we need the full bidiagonal matrix
    // on each process
    //
    // NOTE: lapack::BidiagDQDS expects e to be of length k
    typedef Base<F> Real;
    DistMatrix<Real,STAR,STAR>
      d_STAR_STAR( d_MD_STAR ),
      eHat_STAR_STAR( k, 1, g );
    auto e_STAR_STAR = eHat_STAR_STAR( IR(0,k-1), ALL );
    e_STAR_STAR = e_MD_STAR;

    // Compute the singular values of the bidiagonal matrix
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    bidiag_svd::QRAlg
    ( d_STAR_STAR.Matrix(), e_STAR_STAR.Matrix(), ctrl.bidiagSVDCtrl );
    if( ctrl.time && g.Rank() == 0 )
        Output("Bidiag SVD: ",timer.Stop()," seconds");

    const bool compact = ( ctrl.bidiagSVDCtrl.approach == COMPACT_SVD );
    if( compact )
    {
        const Real twoNorm = ( k==0 ? Real(0) : d_STAR_STAR.Get(0,0) );
        const Real thresh =
          bidiag_svd::APosterioriThreshold
          ( m, n, twoNorm, ctrl.bidiagSVDCtrl );

        Int rank = k;
        for( Int j=0; j<k; ++j )
        {
            if( d_STAR_STAR.Get(j,0) <= thresh )
            {
                rank = j;
                break;
            }
        }
        d_STAR_STAR.Resize( rank, 1 );
    }

    // Copy out the appropriate subset of the singular values
    Copy( d_STAR_STAR, s );

    return info;
}

template<typename F>
SVDInfo GolubReinsch
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    return GolubReinsch( A, s, ctrl );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_GOLUBREINSCH_HPP
