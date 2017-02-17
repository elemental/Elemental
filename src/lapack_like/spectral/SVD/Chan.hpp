/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVD_CHAN_HPP
#define EL_SVD_CHAN_HPP

#include "./GolubReinsch.hpp"

namespace El {
namespace svd {

template<typename Field>
SVDInfo ChanUpper
( Matrix<Field>& A,
  Matrix<Field>& U,
  Matrix<Base<Field>>& s,
  Matrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( ctrl.fullChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )
    typedef Base<Field> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const double heightRatio = ctrl.fullChanRatio;
    const SVDApproach approach = ctrl.bidiagSVDCtrl.approach;
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    if( avoidU && avoidV )
    {
        return SVD( A, s, ctrl );
    }

    SVDInfo info;
    Timer timer;
    if( avoidU )
    {
        if( m > heightRatio*n )
        {
            Matrix<Field> householderScalars;
            Matrix<Real> signature;
            if( ctrl.time )
                timer.Start();
            QR( A, householderScalars, signature );
            if( ctrl.time )
                Output("Chan QR reduction: ",timer.Stop()," seconds");

            Matrix<Field> R;
            auto AT = A( IR(0,n), IR(0,n) );
            R = AT;
            MakeTrapezoidal( UPPER, R );

            // NOTE: This is only appropriate because U is not formed
            info = svd::GolubReinsch( R, U, s, V, ctrl );
        }
        else
        {
            // NOTE: This is only appropriate because U is not formed
            info = svd::GolubReinsch( A, U, s, V, ctrl );
        }
    }
    else
    {
        // This branch handles (avoidU,avoidV) in {(false,false),(false,true)}
        if( m > heightRatio*n )
        {
            Matrix<Field> householderScalars;
            Matrix<Real> signature;
            if( ctrl.time )
                timer.Start();
            QR( A, householderScalars, signature );
            if( ctrl.time )
                Output("Chan QR reduction: ",timer.Stop()," seconds");

            Matrix<Field> R;
            auto AT = A( IR(0,n), IR(0,n) );
            R = AT;
            MakeTrapezoidal( UPPER, R );

            if( approach == FULL_SVD )
            {
                Identity( U, m, m );
                auto UTL = U( IR(0,n), IR(0,n) );
                info = svd::GolubReinsch( R, UTL, s, V, ctrl );
                if( ctrl.time )
                    timer.Start();
                qr::ApplyQ( LEFT, NORMAL, A, householderScalars, signature, U );
                if( ctrl.time )
                    Output("Chan backtransformation: ",timer.Stop()," seconds");
            }
            else
            {
                Zeros( U, m, n );
                auto UT = U( IR(0,n), IR(0,n) );
                info = svd::GolubReinsch( R, UT, s, V, ctrl );
                const Int rank = UT.Width();
                U.Resize( m, rank );
                // (U,s,V) holds an SVD of the R from the QR fact. of original A
                if( ctrl.time )
                    timer.Start();
                qr::ApplyQ( LEFT, NORMAL, A, householderScalars, signature, U );
                if( ctrl.time )
                    Output("Chan backtransformation: ",timer.Stop()," seconds");
            }
        }
        else
        {
            if( approach == FULL_SVD )
            {
                Identity( U, m, m );
                auto UL = U( IR(0,m), IR(0,n) );
                info = svd::GolubReinsch( A, UL, s, V, ctrl );
            }
            else
            {
                info = svd::GolubReinsch( A, U, s, V, ctrl );
            }
        }
    }
    return info;
}

template<typename Field>
SVDInfo ChanUpper
( DistMatrix<Field>& A,
  DistMatrix<Field>& U,
  AbstractDistMatrix<Base<Field>>& s,
  DistMatrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, U, s, V );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( ctrl.fullChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )
    typedef Base<Field> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const double heightRatio = ctrl.fullChanRatio;
    const SVDApproach approach = ctrl.bidiagSVDCtrl.approach;
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    if( avoidU && avoidV )
    {
        return SVD( A, s, ctrl );
    }

    SVDInfo info;
    Timer timer;
    if( avoidU )
    {
        if( m > heightRatio*n )
        {
            DistMatrix<Field,MD,STAR> householderScalars(g);
            DistMatrix<Real,MD,STAR> signature(g);
            if( ctrl.time && g.Rank() == 0 )
                timer.Start();
            QR( A, householderScalars, signature );
            if( ctrl.time && g.Rank() == 0 )
                Output("Chan QR reduction: ",timer.Stop()," seconds");

            DistMatrix<Field> R(g);
            auto AT = A( IR(0,n), IR(0,n) );
            R = AT;
            MakeTrapezoidal( UPPER, R );

            // NOTE: This is only appropriate because U is not formed
            info = svd::GolubReinsch( R, U, s, V, ctrl );
        }
        else
        {
            // NOTE: This is only appropriate because U is not formed
            info = svd::GolubReinsch( A, U, s, V, ctrl );
        }
    }
    else
    {
        // This branch handles (avoidU,avoidV) in {(false,false),(false,true)}
        if( m > heightRatio*n )
        {
            DistMatrix<Field,MD,STAR> householderScalars(g);
            DistMatrix<Real,MD,STAR> signature(g);
            if( ctrl.time && g.Rank() == 0 )
                timer.Start();
            QR( A, householderScalars, signature );
            if( ctrl.time && g.Rank() == 0 )
                Output("Chan QR reduction: ",timer.Stop()," seconds");

            DistMatrix<Field> R(g);
            auto AT = A( IR(0,n), IR(0,n) );
            R = AT;
            MakeTrapezoidal( UPPER, R );

            if( approach == FULL_SVD )
            {
                Identity( U, m, m );
                auto UTL = U( IR(0,n), IR(0,n) );
                info = svd::GolubReinsch( R, UTL, s, V, ctrl );
                if( ctrl.time && g.Rank() == 0 )
                    timer.Start();
                qr::ApplyQ( LEFT, NORMAL, A, householderScalars, signature, U );
                if( ctrl.time && g.Rank() == 0 )
                    Output("Chan backtransformation: ",timer.Stop()," seconds");
            }
            else
            {
                Zeros( U, m, n );
                auto UT = U( IR(0,n), IR(0,n) );
                info = svd::GolubReinsch( R, UT, s, V, ctrl );
                const Int rank = UT.Width();
                U.Resize( m, rank );
                // (U,s,V) holds an SVD of the R from the QR fact. of original A
                if( ctrl.time && g.Rank() == 0 )
                    timer.Start();
                qr::ApplyQ( LEFT, NORMAL, A, householderScalars, signature, U );
                if( ctrl.time && g.Rank() == 0 )
                    Output("Chan backtransformation: ",timer.Stop()," seconds");
            }
        }
        else
        {
            if( approach == FULL_SVD )
            {
                Identity( U, m, m );
                auto UL = U( IR(0,m), IR(0,n) );
                info = svd::GolubReinsch( A, UL, s, V, ctrl );
            }
            else
            {
                info = svd::GolubReinsch( A, U, s, V, ctrl );
            }
        }
    }
    return info;
}

template<typename Field>
SVDInfo ChanUpper
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& UPre,
  AbstractDistMatrix<Base<Field>>& s,
  AbstractDistMatrix<Field>& VPre,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    DistMatrixWriteProxy<Field,Field,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<Field,Field,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return ChanUpper( A, U, s, V, ctrl );
}

template<typename Field>
SVDInfo ChanUpper
( Matrix<Field>& A,
  Matrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( ctrl.valChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint");
    )
    if( A.Height() >= ctrl.valChanRatio*A.Width() )
    {
        qr::ExplicitTriang( A );
        return GolubReinsch( A, s, ctrl );
    }
    else
    {
        return GolubReinsch( A, s, ctrl );
    }
}

template<typename Field>
SVDInfo ChanUpper
( DistMatrix<Field>& A,
  AbstractDistMatrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, s );
      if( ctrl.valChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint");
    )
    if( A.Height() >= ctrl.valChanRatio*A.Width() )
    {
        qr::ExplicitTriang( A );
        return GolubReinsch( A, s, ctrl );
    }
    else
    {
        return GolubReinsch( A, s, ctrl );
    }
}

template<typename Field>
SVDInfo ChanUpper
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    return ChanUpper( A, s, ctrl );
}

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H using Chan's  //
// algorithm. On exit, A is overwritten with U.                               //
//----------------------------------------------------------------------------//

template<typename Field>
SVDInfo Chan
( Matrix<Field>& A,
  Matrix<Field>& U,
  Matrix<Base<Field>>& s,
  Matrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( ctrl.fullChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )
    SVDInfo info;

    // Check if we need to rescale the matrix, and do so if necessary
    // TODO(poulson): Switch to SafeScale
    Base<Field> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        info = svd::ChanUpper( A, U, s, V, ctrl );
    }
    else
    {
        // TODO: Avoid the explicit copy by explicitly forming the Q from LQ
        Matrix<Field> AAdj;
        Adjoint( A, AAdj );
        auto ctrlMod( ctrl );
        ctrlMod.bidiagSVDCtrl.wantU = ctrl.bidiagSVDCtrl.wantV;
        ctrlMod.bidiagSVDCtrl.wantV = ctrl.bidiagSVDCtrl.wantU;
        ctrlMod.bidiagSVDCtrl.accumulateU = ctrl.bidiagSVDCtrl.accumulateV;
        ctrlMod.bidiagSVDCtrl.accumulateV = ctrl.bidiagSVDCtrl.accumulateU;
        info = svd::ChanUpper( AAdj, V, s, U, ctrlMod );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;

    return info;
}

template<typename Field>
SVDInfo Chan
( DistMatrix<Field>& A,
  DistMatrix<Field>& U,
  AbstractDistMatrix<Base<Field>>& s,
  DistMatrix<Field>& V,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, U, s, V );
      if( ctrl.fullChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )
    SVDInfo info;

    // Check if we need to rescale the matrix, and do so if necessary
    // TODO(poulson): Switch to SafeScale
    Base<Field> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        info = svd::ChanUpper( A, U, s, V, ctrl );
    }
    else
    {
        // TODO: Avoid the explicit copy by explicitly forming the Q from LQ
        DistMatrix<Field> AAdj(A.Grid());
        Adjoint( A, AAdj );
        auto ctrlMod( ctrl );
        ctrlMod.bidiagSVDCtrl.wantU = ctrl.bidiagSVDCtrl.wantV;
        ctrlMod.bidiagSVDCtrl.wantV = ctrl.bidiagSVDCtrl.wantU;
        ctrlMod.bidiagSVDCtrl.accumulateU = ctrl.bidiagSVDCtrl.accumulateV;
        ctrlMod.bidiagSVDCtrl.accumulateV = ctrl.bidiagSVDCtrl.accumulateU;
        info = svd::ChanUpper( AAdj, V, s, U, ctrlMod );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;

    return info;
}

template<typename Field>
SVDInfo Chan
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& UPre,
  AbstractDistMatrix<Base<Field>>& s,
  AbstractDistMatrix<Field>& VPre,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    DistMatrixWriteProxy<Field,Field,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<Field,Field,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    return Chan( A, U, s, V, ctrl );
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A.                          //
//----------------------------------------------------------------------------//

template<typename Field>
SVDInfo Chan
( Matrix<Field>& A,
  Matrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SVDInfo info;

    // Check if we need to rescale the matrix, and do so if necessary
    // TODO(poulson): Switch to SafeScale
    Base<Field> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        info = svd::ChanUpper( A, s, ctrl );
    }
    else
    {
        // Explicit formation of the Q from an LQ factorization is not yet
        // optimized
        Matrix<Field> AAdj;
        Adjoint( A, AAdj );
        info = svd::ChanUpper( AAdj, s, ctrl );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;

    return info;
}

template<typename Field>
SVDInfo Chan
( DistMatrix<Field>& A,
  AbstractDistMatrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SVDInfo info;

    // Check if we need to rescale the matrix, and do so if necessary
    // TODO(poulson): Switch to SafeScale
    Base<Field> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        info = svd::ChanUpper( A, s, ctrl );
    }
    else
    {
        // Explicit formation of the Q from an LQ factorization is not yet
        // optimized
        DistMatrix<Field> AAdj( A.Grid() );
        Adjoint( A, AAdj );
        info = svd::ChanUpper( AAdj, s, ctrl );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;

    return info;
}

template<typename Field>
SVDInfo Chan
( AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Base<Field>>& s,
  const SVDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    return Chan( A, s, ctrl );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_CHAN_HPP
