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

template<typename F>
inline void
ChanUpper
( DistMatrix<F>& A,
  DistMatrix<F>& U,
  ElementalMatrix<Base<F>>& s, 
  DistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("svd::ChanUpper [DistMatrix Decomp]");
      AssertSameGrids( A, U, s, V );
      if( A.Height() < A.Width() )
          LogicError("A must be at least as tall as it is wide");
      if( ctrl.fullChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int m = A.Height();
    const Int n = A.Width();
    const double heightRatio = ctrl.fullChanRatio;
    const SVDApproach approach = ctrl.approach;
    const bool avoidU = ctrl.avoidComputingU;
    const bool avoidV = ctrl.avoidComputingV;
    if( avoidU && avoidV )
    {
        SVD( A, s, ctrl );
        return;
    }

    Timer timer;
    if( avoidU )
    {
        if( m > heightRatio*n )
        {
            DistMatrix<F,MD,STAR> t(g);
            DistMatrix<Real,MD,STAR> d(g);
            if( ctrl.time && g.Rank() == 0 )
                timer.Start();
            QR( A, t, d );
            if( ctrl.time && g.Rank() == 0 )
                Output("Chan QR reduction: ",timer.Stop()," seconds");

            DistMatrix<F> R(g);
            auto AT = A( IR(0,n), IR(0,n) );
            R = AT;
            MakeTrapezoidal( UPPER, R );

            // NOTE: This is only appropriate because U is not formed
            svd::GolubReinsch( R, U, s, V, ctrl );
        }
        else
        {
            // NOTE: This is only appropriate because U is not formed
            svd::GolubReinsch( A, U, s, V, ctrl );
        }
    }
    else
    {
        // This branch handles (avoidU,avoidV) in {(false,false),(false,true)}
        if( m > heightRatio*n )
        {
            DistMatrix<F,MD,STAR> t(g);
            DistMatrix<Real,MD,STAR> d(g);
            if( ctrl.time && g.Rank() == 0 )
                timer.Start();
            QR( A, t, d );
            if( ctrl.time && g.Rank() == 0 )
                Output("Chan QR reduction: ",timer.Stop()," seconds");

            DistMatrix<F> R(g);
            auto AT = A( IR(0,n), IR(0,n) );
            R = AT;
            MakeTrapezoidal( UPPER, R );

            if( approach == FULL_SVD )
            {
                Identity( U, m, m );
                auto UTL = U( IR(0,n), IR(0,n) );
                svd::GolubReinsch( R, UTL, s, V, ctrl );
                if( ctrl.time && g.Rank() == 0 )
                    timer.Start();
                qr::ApplyQ( LEFT, NORMAL, A, t, d, U );
                if( ctrl.time && g.Rank() == 0 )
                    Output("Chan backtransformation: ",timer.Stop()," seconds");
            }
            else
            {
                Zeros( U, m, n );
                auto UT = U( IR(0,n), IR(0,n) );
                svd::GolubReinsch( R, UT, s, V, ctrl );
                const Int rank = UT.Width();
                U.Resize( m, rank );
                // (U,s,V) holds an SVD of the R from the QR fact. of original A
                if( ctrl.time && g.Rank() == 0 )
                    timer.Start();
                qr::ApplyQ( LEFT, NORMAL, A, t, d, U );
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
                svd::GolubReinsch( A, UL, s, V, ctrl );
            }
            else
            {
                svd::GolubReinsch( A, U, s, V, ctrl );
            }
        }
    }
}

template<typename F>
inline void
ChanUpper
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& UPre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::ChanUpper [ElementalMatrix Decomp]"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    ChanUpper( A, U, s, V, ctrl );
}

template<typename F>
inline void
ChanUpper
( DistMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("svd::ChanUpper [DistMatrix values]");
      AssertSameGrids( A, s );
      if( ctrl.valChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint");
    )
    if( A.Height() >= ctrl.valChanRatio*A.Width() )
    {
        qr::ExplicitTriang( A );
        GolubReinsch( A, s, ctrl );
    }
    else
    {
        GolubReinsch( A, s, ctrl );
    }
}

template<typename F>
inline void
ChanUpper
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::ChanUpper [ElementalMatrix values]"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    ChanUpper( A, s, ctrl );
}

//----------------------------------------------------------------------------//
// Grab the full SVD of the general matrix A, A = U diag(s) V^H using Chan's  //
// algorithm. On exit, A is overwritten with U.                               //
//----------------------------------------------------------------------------//

template<typename F>
inline void
Chan
( DistMatrix<F>& A,
  DistMatrix<F>& U,
  ElementalMatrix<Base<F>>& s, 
  DistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("svd::Chan [DistMatrix Decomp]");
      AssertSameGrids( A, U, s, V );
      if( ctrl.fullChanRatio <= 1.0 )
          LogicError("Nonsensical switchpoint for SVD");
    )

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::ChanUpper( A, U, s, V, ctrl );
    }
    else
    {
        // TODO: Avoid the explicit copy by explicitly forming the Q from LQ
        DistMatrix<F> AAdj(A.Grid());
        Adjoint( A, AAdj );
        auto ctrlMod( ctrl );
        ctrlMod.avoidComputingU = ctrl.avoidComputingV;
        ctrlMod.avoidComputingV = ctrl.avoidComputingU;
        svd::ChanUpper( AAdj, V, s, U, ctrlMod );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;
}

template<typename F>
inline void
Chan
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& UPre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::Chan [ElementalMatrix Decomp]"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    Chan( A, U, s, V, ctrl );
}

//----------------------------------------------------------------------------//
// Grab the singular values of the general matrix A.                          //
//----------------------------------------------------------------------------//

template<typename F>
inline void
Chan
( DistMatrix<F>& A,
  ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::Chan [DistMatrix values]"))

    // Check if we need to rescale the matrix, and do so if necessary
    Base<F> scale;
    bool needRescaling = svd::CheckScale( A, scale );
    if( needRescaling )
        A *= scale;

    // TODO: Switch between different algorithms. For instance, starting 
    //       with a QR decomposition of tall-skinny matrices.
    if( A.Height() >= A.Width() )
    {
        svd::ChanUpper( A, s, ctrl );
    }
    else
    {
        // Explicit formation of the Q from an LQ factorization is not yet
        // optimized
        DistMatrix<F> AAdj( A.Grid() );
        Adjoint( A, AAdj );
        svd::ChanUpper( AAdj, s, ctrl );
    }

    // Rescale the singular values if necessary
    if( needRescaling )
        s *= 1/scale;
}

template<typename F>
inline void
Chan
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s, 
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::Chan [ElementalMatrix values]"))
    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    Chan( A, s, ctrl );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_CHAN_HPP
