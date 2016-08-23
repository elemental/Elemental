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
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    if( avoidU && avoidV )
    {
        return SVD( A, s, ctrl );
    }
    SVDInfo info;

    // Bidiagonalize A
    Timer timer;
    Matrix<F> phaseP, phaseQ;
    if( ctrl.time )
        timer.Start();
    Bidiag( A, phaseP, phaseQ );
    if( ctrl.time )
        Output("Reduction to bidiagonal: ",timer.Stop()," seconds");

    // Compute the SVD of the bidiagonal matrix.
    // (We can guarantee that accumulation was not requested.)
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    const UpperOrLower uplo = ( m>=n ? UPPER : LOWER );
    auto mainDiag = GetRealPartOfDiagonal( A );
    auto offDiag = GetRealPartOfDiagonal( A, offdiagonal );
    if( ctrl.time )
        timer.Start();
    info.bidiagSVDInfo =
      BidiagSVD( uplo, mainDiag, offDiag, U, s, V, ctrl.bidiagSVDCtrl );
    if( ctrl.time )
        Output("Bidiag SVD: ",timer.Stop()," seconds");

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
    const bool avoidU = !ctrl.bidiagSVDCtrl.wantU;
    const bool avoidV = !ctrl.bidiagSVDCtrl.wantV;
    const Grid& g = A.Grid();
    if( avoidU && avoidV )
    {
        return SVD( A, s, ctrl );
    }
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
    const UpperOrLower uplo = ( m>=n ? UPPER : LOWER );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    auto mainDiag = GetRealPartOfDiagonal(A);
    auto offDiag = GetRealPartOfDiagonal(A,offdiagonal);

    // Run the bidiagonal SVD
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    info.bidiagSVDInfo =
      BidiagSVD( uplo, mainDiag, offDiag, U, s, V, ctrl.bidiagSVDCtrl );
    if( ctrl.time )
    {
        mpi::Barrier( g.Comm() );
        if( g.Rank() == 0 )
            Output("Bidiag SVD: ",timer.Stop()," seconds");
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
    SVDInfo info;

    // Bidiagonalize A
    Timer timer;
    Matrix<F> phaseP, phaseQ;
    if( ctrl.time )
        timer.Start();
    Bidiag( A, phaseP, phaseQ );
    if( ctrl.time )
        Output("Reduction to bidiagonal: ",timer.Stop()," seconds");

    // Compute the singular values of the bidiagonal matrix
    const UpperOrLower uplo = ( m>=n ? UPPER : LOWER );
    const Int offdiagonal = ( uplo==UPPER ? 1 : -1 );
    auto mainDiag = GetRealPartOfDiagonal( A );
    auto offDiag = GetRealPartOfDiagonal( A, offdiagonal );
    if( ctrl.time )
        timer.Start();
    info.bidiagSVDInfo =
      BidiagSVD( uplo, mainDiag, offDiag, s, ctrl.bidiagSVDCtrl );
    if( ctrl.time )
        Output("Bidiag SVD: ",timer.Stop()," seconds");

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
    const UpperOrLower uplo = ( m>=n ? UPPER : LOWER );
    const Int offdiagonal = ( uplo==UPPER ? 1 : -1 );
    auto mainDiag = GetRealPartOfDiagonal(A);
    auto offDiag = GetRealPartOfDiagonal(A,offdiagonal);
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    info.bidiagSVDInfo =
      BidiagSVD( uplo, mainDiag, offDiag, s, ctrl.bidiagSVDCtrl );
    if( ctrl.time && g.Rank() == 0 )
        Output("Bidiag SVD: ",timer.Stop()," seconds");

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
