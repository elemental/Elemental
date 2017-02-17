/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HESS_SCHUR_AED_NIBBLE_HPP
#define EL_HESS_SCHUR_AED_NIBBLE_HPP

#include "./SpikeDeflation.hpp"

namespace El {
namespace hess_schur {
namespace aed {

// The spike value will be overwritten
template<typename Real>
AEDInfo NibbleHelper
( Matrix<Real>& H,
  Real& spikeValue,
  Matrix<Complex<Real>>& w,
  Matrix<Real>& V,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    AEDInfo info;

    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    Zeros( V, 0, 0 );
    if( n == 1 )
    {
        w(0) = H(0,0);
        if( Abs(spikeValue) <= Max( smallNum, ulp*Abs(w(0).real()) ) )
        {
            // The offdiagonal entry was small enough to deflate
            info.numDeflated = 1;
            spikeValue = zero;
        }
        else
        {
            // The offdiagonal entry was too large to deflate
            info.numShiftCandidates = 1;
        }
        return info;
    }

    // NOTE(poulson): We could only copy the upper-Hessenberg portion of H
    auto T( H ); // TODO(poulson): Reuse this matrix?
    Identity( V, n, n );
    auto ctrlSub( ctrl );
    ctrlSub.winBeg = 0;
    ctrlSub.winEnd = n;
    ctrlSub.fullTriangle = true;
    ctrlSub.wantSchurVecs = true;
    ctrlSub.demandConverged = false;
    ctrlSub.alg = ( ctrl.recursiveAED ? HESSENBERG_SCHUR_AED
                                      : HESSENBERG_SCHUR_MULTIBULGE );
    auto infoSub = HessenbergSchur( T, w, V, ctrlSub );
    EL_DEBUG_ONLY(
      if( infoSub.numUnconverged != 0 )
          Output(infoSub.numUnconverged," eigenvalues did not converge");
    )

    vector<Real> work(2*n);
    info = SpikeDeflation( T, V, spikeValue, infoSub.numUnconverged, work );
    if( ctrl.progress )
    {
        if( info.numUnconverged > 0 )
            Output
            ("  ",info.numUnconverged," AED eigenvalues did not converge");
        Output("  ",info.numDeflated," of ",n," AED eigenvalues deflated");
    }
    const Int spikeSize = info.numUnconverged + info.numShiftCandidates;
    if( spikeSize == 0 )
    {
        // The entire spike has deflated
        spikeValue = zero;
    }
    if( info.numDeflated > 0 )
    {
        // TODO(poulson):
        //
        // Follow LAPACK's lead and sort the diagonal blocks of T from largest
        // to smallest magnitude to improve the accuracy of the algorithm on
        // graded matrices.
        //
        // Perhaps this should be its own subroutine that accepts a SortType,
        // which can take the values {ASCENDING, DESCENDING, UNSORTED}.
    }

    // Reform the eigenvalues and shift candidates by looping over the converged
    // eigenvalues from last to first
    for( Int i=n-1; i>=info.numUnconverged; )
    {
        if( i == info.numUnconverged || T(i,i-1) == zero )
        {
            // 1x1 block
            w(i) = T(i,i);
            i -= 1;
        }
        else
        {
            // 2x2 block
            Real alpha00 = T(i-1,i-1);
            Real alpha10 = T(i,  i-1);
            Real alpha01 = T(i-1,i  );
            Real alpha11 = T(i,  i  );
            schur::TwoByTwo
            ( alpha00, alpha01,
              alpha10, alpha11,
              w(i-1), w(i) );
            i -= 2;
        }
    }
    if( spikeSize >= n && spikeValue != zero )
    {
        // There were no deflations and we have a coupling
        V.Resize( 0, 0 );
        return info;
    }

    // Either we deflated at least one eigenvalue or we can simply
    // rotate the deflation window into Schur form
    auto spikeInd = IR(0,spikeSize);
    auto TTL = T(spikeInd,spikeInd);
    auto TTR = T(spikeInd,IR(spikeSize,END));
    auto VL = V(ALL,spikeInd);

    Matrix<Real> householderScalarsT;
    if( spikeSize > 1 && spikeValue != zero )
    {
        // The spike needs to be reduced to length one while maintaining
        // the Hessenberg form of the deflation window
        for( Int i=0; i<spikeSize; ++i )
            work[i] = V(0,i);

        // Compute a Householder reflector for condensing the spike
        Real beta = work[0];
        Real tau = lapack::Reflector( spikeSize, beta, &work[1], 1 );
        work[0] = Real(1);

        // Force T to be upper Hessenberg
        MakeTrapezoidal( UPPER, T, -1 );

        lapack::ApplyReflector
        ( true, spikeSize, n,
          &work[0], 1, tau,
          T.Buffer(), T.LDim(),
          &work[n] );

        lapack::ApplyReflector
        ( false, spikeSize, spikeSize,
          &work[0], 1, tau,
          TTL.Buffer(), TTL.LDim(),
          &work[n] );

        lapack::ApplyReflector
        ( false, n, spikeSize,
          &work[0], 1, tau,
          VL.Buffer(), VL.LDim(),
          &work[n] );

        Hessenberg( UPPER, TTL, householderScalarsT );
        hessenberg::ApplyQ
        ( LEFT, UPPER, ADJOINT, TTL, householderScalarsT, TTR );
    }

    spikeValue *= V(0,0);
    // NOTE(poulson): We could copy only the upper Hessenberg part of T
    H = T;
    MakeTrapezoidal( UPPER, H, -1 );

    if( spikeSize > 1 && spikeValue != zero )
    {
        hessenberg::ApplyQ
        ( RIGHT, UPPER, NORMAL, TTL, householderScalarsT, VL );
    }

    return info;
}

template<typename Real>
AEDInfo NibbleHelper
( Matrix<Complex<Real>>& H,
  Complex<Real>& spikeValue,
  Matrix<Complex<Real>>& w,
  Matrix<Complex<Real>>& V,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    typedef Complex<Real> Field;
    const Int n = H.Height();
    AEDInfo info;

    const Real zero(0);
    const Real ulp = limits::Precision<Real>();
    const Real safeMin = limits::SafeMin<Real>();
    const Real smallNum = safeMin*(Real(n)/ulp);

    Zeros( V, 0, 0 );
    if( n == 1 )
    {
        w(0) = H(0,0);
        if( OneAbs(spikeValue) <= Max( smallNum, ulp*OneAbs(w(0)) ) )
        {
            // The offdiagonal entry was small enough to deflate
            info.numDeflated = 1;
            spikeValue = zero;
        }
        else
        {
            // The offdiagonal entry was too large to deflate
            info.numShiftCandidates = 1;
        }
        return info;
    }

    // NOTE(poulson): We could only copy the upper-Hessenberg portion of H
    auto T( H ); // TODO(poulson): Reuse this matrix?
    Identity( V, n, n );
    auto ctrlSub( ctrl );
    ctrlSub.winBeg = 0;
    ctrlSub.winEnd = n;
    ctrlSub.fullTriangle = true;
    ctrlSub.wantSchurVecs = true;
    ctrlSub.demandConverged = false;
    ctrlSub.alg = ( ctrl.recursiveAED ? HESSENBERG_SCHUR_AED
                                      : HESSENBERG_SCHUR_MULTIBULGE );
    auto infoSub = HessenbergSchur( T, w, V, ctrlSub );
    EL_DEBUG_ONLY(
      if( infoSub.numUnconverged != 0 )
          Output(infoSub.numUnconverged," eigenvalues did not converge");
    )

    vector<Field> work(2*n);
    info = SpikeDeflation( T, V, spikeValue, infoSub.numUnconverged, work );
    if( ctrl.progress )
    {
        if( info.numUnconverged > 0 )
            Output
            ("  ",info.numUnconverged," AED eigenvalues did not converge");
        Output("  ",info.numDeflated," of ",n," AED eigenvalues deflated");
    }
    const Int spikeSize = info.numUnconverged + info.numShiftCandidates;
    if( spikeSize == 0 )
    {
        // The entire spike has deflated
        spikeValue = zero;
    }
    if( info.numDeflated > 0 )
    {
        // TODO(poulson):
        //
        // Follow LAPACK's lead and sort the diagonal of T from largest
        // to smallest magnitude to improve the accuracy of the algorithm on
        // graded matrices.
        //
        // Perhaps this should be its own subroutine that accepts a SortType,
        // which can take the values {ASCENDING, DESCENDING, UNSORTED}.
    }

    // Reform the eigenvalues and shift candidates by looping over the converged
    // eigenvalues from last to first
    for( Int i=n-1; i>=info.numUnconverged; --i )
        w(i) = T(i,i);

    if( spikeSize >= n && spikeValue != zero )
    {
        V.Resize( 0, 0 );
        return info;
    }

    // Either we deflated at least one eigenvalue or we can simply
    // rotate the deflation window into Schur form
    auto spikeInd = IR(0,spikeSize);
    auto TTL = T(spikeInd,spikeInd);
    auto TTR = T(spikeInd,IR(spikeSize,END));
    auto VL = V(ALL,spikeInd);
    Matrix<Field> householderScalarsT;

    if( spikeSize > 1 && spikeValue != zero )
    {
        // The spike needs to be reduced to length one while maintaining
        // the Hessenberg form of the deflation window
        for( Int i=0; i<spikeSize; ++i )
            work[i] = Conj(V(0,i));

        // Compute a Householder reflector for condensing the spike
        Field beta = work[0];
        Field tau = lapack::Reflector( spikeSize, beta, &work[1], 1 );
        work[0] = Real(1);

        // Force T to be upper Hessenberg
        MakeTrapezoidal( UPPER, T, -1 );

        lapack::ApplyReflector
        ( true, spikeSize, n,
          &work[0], 1, tau,
          T.Buffer(), T.LDim(),
          &work[n] );

        lapack::ApplyReflector
        ( false, spikeSize, spikeSize,
          &work[0], 1, Conj(tau),
          TTL.Buffer(), TTL.LDim(),
          &work[n] );

        lapack::ApplyReflector
        ( false, n, spikeSize,
          &work[0], 1, Conj(tau),
          VL.Buffer(), VL.LDim(),
          &work[n] );

        Hessenberg( UPPER, TTL, householderScalarsT );
        hessenberg::ApplyQ
        ( LEFT, UPPER, ADJOINT, TTL, householderScalarsT, TTR );
    }

    spikeValue *= Conj(V(0,0));
    // NOTE(poulson): We could copy only the upper Hessenberg part of T
    H = T;
    MakeTrapezoidal( UPPER, H, -1 );

    if( spikeSize > 1 && spikeValue != zero )
    {
        hessenberg::ApplyQ
        ( RIGHT, UPPER, NORMAL, TTL, householderScalarsT, VL );
    }

    return info;
}

template<typename Field>
AEDInfo Nibble
( Matrix<Field>& H,
  Int deflationSize,
  Matrix<Complex<Base<Field>>>& w,
  Matrix<Field>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );
    AEDInfo info;

    if( winBeg > winEnd )
        return info;
    if( deflationSize < 1 )
        return info;

    const Int blockSize = Min( deflationSize, winEnd-winBeg );
    const Int deflateBeg = winEnd-blockSize;
    auto deflateInd = IR(deflateBeg,winEnd);

    // If the deflation window touches the beginning of the full window,
    // then there is no spike
    Field spikeValue =
      ( deflateBeg==winBeg ? Field(0) : H(deflateBeg,deflateBeg-1) );

    Matrix<Field> V;
    auto HDefl = H( deflateInd, deflateInd );
    auto wDefl = w( deflateInd, ALL );
    info = NibbleHelper( HDefl, spikeValue, wDefl, V, ctrl );
    if( deflateBeg > winBeg )
        H(deflateBeg,deflateBeg-1) = spikeValue;
    if( V.Height() == 0 && V.Width() == 0 )
    {
        // It was signalled that no transformation was required
        return info;
    }

    if( ctrl.fullTriangle )
    {
        auto H12 = H( deflateInd, IR(winEnd,END) );
        multibulge::TransformRows( V, H12 );
    }

    const Int applyBeg = ( ctrl.fullTriangle ? 0 : winBeg );
    auto H01 = H( IR(applyBeg,deflateBeg), deflateInd );
    multibulge::TransformColumns( V, H01 );

    if( ctrl.wantSchurVecs )
    {
        auto Z1 = Z(ALL,deflateInd);
        multibulge::TransformColumns( V, Z1 );
    }

    return info;
}

template<typename Field>
AEDInfo Nibble
( DistMatrix<Field,MC,MR,BLOCK>& H,
  Int deflationSize,
  DistMatrix<Complex<Base<Field>>,STAR,STAR>& w,
  DistMatrix<Field,MC,MR,BLOCK>& Z,
  const HessenbergSchurCtrl& ctrl )
{
    EL_DEBUG_CSE
    const Int n = H.Height();
    const Grid& grid = H.Grid();
    AEDInfo info;

    const Int winBeg = ( ctrl.winBeg==END ? n : ctrl.winBeg );
    const Int winEnd = ( ctrl.winEnd==END ? n : ctrl.winEnd );

    if( winBeg > winEnd )
        return info;
    if( deflationSize < 1 )
        return info;

    const Int blockSize = Min( deflationSize, winEnd-winBeg );
    const Int deflateBeg = winEnd-blockSize;
    auto deflateInd = IR(deflateBeg,winEnd);
    auto HDefl = H( deflateInd, deflateInd );
    auto wDefl = w( deflateInd, ALL );

    const int owner = HDefl.Owner(0,0);
    DistMatrix<Field,CIRC,CIRC> HDefl_CIRC_CIRC( grid, owner );
    HDefl_CIRC_CIRC = HDefl;
    Field spikeValue =
      ( deflateBeg==winBeg ? Field(0) : H.Get(deflateBeg,deflateBeg-1) );
    Int VSize = 0;
    Matrix<Field> V;
    if( HDefl_CIRC_CIRC.CrossRank() == HDefl_CIRC_CIRC.Root() )
    {
        info =
          NibbleHelper
          ( HDefl_CIRC_CIRC.Matrix(), spikeValue, wDefl.Matrix(), V, ctrl );
        VSize = V.Height();
    }
    El::Broadcast( wDefl, HDefl_CIRC_CIRC.CrossComm(), HDefl_CIRC_CIRC.Root() );

    if( deflateBeg > winBeg )
    {
        mpi::Broadcast
        ( spikeValue, HDefl_CIRC_CIRC.Root(), HDefl_CIRC_CIRC.CrossComm() );
        H.Set( deflateBeg, deflateBeg-1, spikeValue );
    }

    // TODO(poulson): Combine these broadcasts?
    mpi::Broadcast
    ( info.numDeflated, HDefl_CIRC_CIRC.Root(), HDefl_CIRC_CIRC.CrossComm() );
    mpi::Broadcast
    ( info.numShiftCandidates,
      HDefl_CIRC_CIRC.Root(), HDefl_CIRC_CIRC.CrossComm() );
    mpi::Broadcast
    ( VSize, HDefl_CIRC_CIRC.Root(), HDefl_CIRC_CIRC.CrossComm() );
    if( VSize == 0 )
    {
        return info;
    }
    V.Resize( VSize, VSize );
    El::Broadcast( V, HDefl_CIRC_CIRC.CrossComm(), HDefl_CIRC_CIRC.Root() );
    HDefl = HDefl_CIRC_CIRC;

    if( ctrl.fullTriangle )
    {
        auto H12 = H( deflateInd, IR(winEnd,END) );
        multibulge::TransformRows( V, H12 );
    }

    const Int applyBeg = ( ctrl.fullTriangle ? 0 : winBeg );
    auto H01 = H( IR(applyBeg,deflateBeg), deflateInd );
    multibulge::TransformColumns( V, H01 );

    if( ctrl.wantSchurVecs )
    {
        auto Z1 = Z(ALL,deflateInd);
        multibulge::TransformColumns( V, Z1 );
    }

    return info;
}

} // namespace aed
} // namespace hess_schur
} // namespace El

#endif // ifndef EL_HESS_SCHUR_AED_NIBBLE_HPP
