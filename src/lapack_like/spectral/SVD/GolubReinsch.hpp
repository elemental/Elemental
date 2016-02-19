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
inline void
GolubReinsch
( DistMatrix<F>& A,
  DistMatrix<F>& U,
  ElementalMatrix<Base<F>>& s, 
  DistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch [DistMatrix Decomp]"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    const char uplo = ( m>=n ? 'U' : 'L' );
    const bool avoidU = ctrl.avoidComputingU;
    const bool avoidV = ctrl.avoidComputingV;
    const Grid& g = A.Grid();
    if( avoidU && avoidV )
    {
        SVD( A, s, ctrl );
        return;
    }

    // Bidiagonalize A
    Timer timer;
    DistMatrix<F,STAR,STAR> tP(g), tQ(g);
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    Bidiag( A, tP, tQ );
    if( ctrl.time && g.Rank() == 0 )
        Output("Reduction to bidiagonal: ",timer.Stop()," seconds");

    // Grab copies of the diagonal and sub/super-diagonal of A
    auto d_MD_STAR = GetRealPartOfDiagonal(A);
    auto e_MD_STAR = GetRealPartOfDiagonal(A,offdiagonal);

    // NOTE: lapack::BidiagQRAlg expects e to be of length k
    typedef Base<F> Real;
    DistMatrix<Real,STAR,STAR>
      d_STAR_STAR( d_MD_STAR ),
      eHat_STAR_STAR( k, 1, g );
    auto e_STAR_STAR = eHat_STAR_STAR( IR(0,k-1), ALL );
    e_STAR_STAR = e_MD_STAR;

    // Initialize U and VAdj to the appropriate identity matrices
    DistMatrix<F,VC,STAR> U_VC_STAR( g );
    if( !avoidU )
    {
        U_VC_STAR.AlignWith( A );
        Identity( U_VC_STAR, m, k );
    }

    DistMatrix<F,STAR,VC> VAdj_STAR_VC( g );
    if( !avoidV )
    {
        VAdj_STAR_VC.AlignWith( V );
        Identity( VAdj_STAR_VC, k, n );
    }

    // TODO: If compact SVD, identify the rank with DQDS first?

    // Compute the SVD of the bidiagonal matrix and accumulate the Givens
    // rotations into our local portion of U and VAdj
    Matrix<F>& ULoc = U_VC_STAR.Matrix();
    Matrix<F>& VAdjLoc = VAdj_STAR_VC.Matrix();
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    lapack::BidiagQRAlg
    ( uplo, k, VAdjLoc.Width(), ULoc.Height(),
      d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(), 
      VAdjLoc.Buffer(), VAdjLoc.LDim(), 
      ULoc.Buffer(), ULoc.LDim() );
    if( ctrl.time )
    {
        mpi::Barrier( g.Comm() );
        if( g.Rank() == 0 )
            Output("BidiagQRAlg: ",timer.Stop()," seconds");
    }

    Int rank = k;
    const bool compact = ( ctrl.approach == COMPACT_SVD );
    if( compact )
    {
        const Real twoNorm = ( k==0 ? Real(0) : d_STAR_STAR.Get(0,0) );
        // Use Max(m,n)*twoNorm*eps unless a manual tolerance is specified
        Real thresh = Max(m,n)*twoNorm*limits::Epsilon<Real>();
        if( ctrl.tol != Real(0) )
        {
            if( ctrl.relative )
                thresh = twoNorm*ctrl.tol;
            else
                thresh = ctrl.tol;
        }
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
        if( !avoidV ) VAdj_STAR_VC.Resize( rank, n );
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
        if( !avoidV ) Adjoint( VAdj_STAR_VC, V );
    }
    else
    {
        if( !avoidU ) U = U_VC_STAR;
        if( !avoidV )
        {
            auto VAdjL_STAR_VC = VAdj_STAR_VC( IR(0,rank), IR(0,m) );
            auto VT = V( IR(0,m  ), ALL );
            auto VB = V( IR(m,END), ALL );
            Adjoint( VAdjL_STAR_VC, VT );
            Zero( VB );
        }
    }

    // Backtransform U and V
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    if( !avoidU ) bidiag::ApplyQ( LEFT, NORMAL, A, tQ, U );
    if( !avoidV ) bidiag::ApplyP( LEFT, NORMAL, A, tP, V );
    if( ctrl.time && g.Rank() == 0 )
        Output("GolubReinsch backtransformation: ",timer.Stop()," seconds");
}

template<typename F>
inline void
GolubReinsch
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& UPre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch [ElementalMatrix Decomp]"))
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    GolubReinsch( A, U, s, V, ctrl );
}

#ifdef EL_HAVE_FLA_BSVD
template<typename F>
inline void
GolubReinschFlame
( DistMatrix<F>& A,
  DistMatrix<F>& U,
  ElementalMatrix<Base<F>>& s, 
  DistMatrix<F>& V,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinschFlame [DistMatrix Decomp]"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    const bool avoidU = ctrl.avoidComputingU;
    const bool avoidV = ctrl.avoidComputingV;
    const Grid& g = A.Grid();
    if( avoidU && avoidV )
    {
        SVD( A, s, ctrl );
        return;
    }

    // Bidiagonalize A
    Timer timer;
    DistMatrix<F,STAR,STAR> tP(g), tQ(g);
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    Bidiag( A, tP, tQ );
    if( ctrl.time && g.Rank() == 0 )
        Output("Reduction to bidiagonal: ",timer.Stop()," seconds");

    // Grab copies of the diagonal and sub/super-diagonal of A
    auto d_MD_STAR = GetRealPartOfDiagonal(A);
    auto e_MD_STAR = GetRealPartOfDiagonal(A,offdiagonal);

    // In order to use serial QR kernels, we need the full bidiagonal matrix
    // on each process
    typedef Base<F> Real;
    DistMatrix<Real,STAR,STAR>
      d_STAR_STAR( d_MD_STAR ),
      e_STAR_STAR( e_MD_STAR );

    // Initialize U and V to the appropriate identity matrices
    DistMatrix<F,VC,STAR> U_VC_STAR(g), V_VC_STAR(g);
    if( !avoidU )
    {
        U_VC_STAR.AlignWith( A );
        Identity( U_VC_STAR, m, k );
    }
    if( !avoidV )
    {
        V_VC_STAR.AlignWith( V );
        Identity( V_VC_STAR, n, k );
    }

    // Since libFLAME, to the best of my current knowledge, only supports the
    // upper-bidiagonal case, we may instead work with the adjoint in the 
    // lower-bidiagonal case.
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    if( m >= n )
    {
        flame::BidiagSVD
        ( k, U_VC_STAR.LocalHeight(), V_VC_STAR.LocalHeight(),
          d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          U_VC_STAR.Buffer(), U_VC_STAR.LDim(),
          V_VC_STAR.Buffer(), V_VC_STAR.LDim() );
    }
    else
    {
        flame::BidiagSVD
        ( k, V_VC_STAR.LocalHeight(), U_VC_STAR.LocalHeight(),
          d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer(),
          V_VC_STAR.Buffer(), V_VC_STAR.LDim(),
          U_VC_STAR.Buffer(), U_VC_STAR.LDim() );
    }
    if( ctrl.time )
    {
        mpi::Barrier( g.Comm() );
        if( g.Rank() == 0 )
            Output("BidiagQRAlg: ",timer.Stop()," seconds");
    }

    Int rank = k;
    const bool compact = ( ctrl.approach == COMPACT_SVD );
    if( compact )
    {
        const Real twoNorm = ( k==0 ? Real(0) : d_STAR_STAR.Get(0,0) );
        // Use Max(m,n)*twoNorm*eps unless a manual tolerance is specified
        Real thresh = Max(m,n)*twoNorm*limits::Epsilon<Real>();
        if( ctrl.tol != Real(0) )
        {
            if( ctrl.relative )
                thresh = twoNorm*ctrl.tol;
            else
                thresh = ctrl.tol;
        }
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
            auto UT_VC_STAR = U_VC_STAR( IR(0,n), IR(0,rank) );
            auto UT = U( IR(0,n), ALL );
            auto UB = U( IR(n,END), ALL );
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
            auto VT = V( IR(0,m  ), ALL );
            auto VB = V( IR(m,END), ALL ); 
            VT = VT_VC_STAR;
            Zero( VB );
        }
    }

    // Backtransform U and V
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    if( !avoidU ) bidiag::ApplyQ( LEFT, NORMAL, A, tQ, U );
    if( !avoidV ) bidiag::ApplyP( LEFT, NORMAL, A, tP, V );
    if( ctrl.time && g.Rank() == 0 )
        Output("GolubReinsch backtransformation: ",timer.Stop()," seconds");
}

template<typename F>
inline void
GolubReinschFlame
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& UPre,
  ElementalMatrix<Base<F>>& s, 
  ElementalMatrix<F>& VPre,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinschFlame [ElementalMatrix Decomp]"))
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixWriteProxy<F,F,MC,MR> VProx( VPre );
    auto& A = AProx.Get();
    auto& U = UProx.Get();
    auto& V = VProx.Get();
    GolubReinschFlame( A, U, s, V, ctrl );
}

template<>
inline void
GolubReinsch
( ElementalMatrix<double>& A,
  ElementalMatrix<double>& U,
  ElementalMatrix<double>& s, 
  ElementalMatrix<double>& V,
  const SVDCtrl<double>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch<double> [ElementalMatrix Decomp]"))
    if( ctrl.avoidLibflame )
        GolubReinsch( A, U, s, V, ctrl );
    else
        GolubReinschFlame( A, U, s, V, ctrl );
}

template<>
inline void
GolubReinsch
( ElementalMatrix<Complex<double>>& A,
  ElementalMatrix<Complex<double>>& U,
  ElementalMatrix<double>& s, 
  ElementalMatrix<Complex<double>>& V,
  const SVDCtrl<double>& ctrl )
{
    DEBUG_ONLY(
      CSE cse("svd::GolubReinsch<Complex<double>> [ElementalMatrix Decomp]")
    )
    if( ctrl.avoidLibflame )
        GolubReinsch( A, U, s, V, ctrl );
    else
        GolubReinschFlame( A, U, s, V, ctrl );
}
#endif // EL_HAVE_FLA_BSVD

template<typename F>
inline void
GolubReinsch
( DistMatrix<F>& A,
  ElementalMatrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch [DistMatrix values]"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = Min( m, n );
    const Int offdiagonal = ( m>=n ? 1 : -1 );
    const Grid& g = A.Grid();

    // Bidiagonalize A
    Timer timer;
    DistMatrix<F,STAR,STAR> tP(g), tQ(g);
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    Bidiag( A, tP, tQ );
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

    // Compute the singular values of the bidiagonal matrix via DQDS
    if( ctrl.time && g.Rank() == 0 )
        timer.Start();
    lapack::BidiagDQDS( k, d_STAR_STAR.Buffer(), e_STAR_STAR.Buffer() );
    if( ctrl.time && g.Rank() == 0 )
        Output("DQDS: ",timer.Stop()," seconds");
    const bool compact = ( ctrl.approach == COMPACT_SVD );
    if( compact )
    {
        const Real twoNorm = ( k==0 ? Real(0) : d_STAR_STAR.Get(0,0) );
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
}

template<typename F>
inline void
GolubReinsch
( ElementalMatrix<F>& APre,
  ElementalMatrix<Base<F>>& s,
  const SVDCtrl<Base<F>>& ctrl )
{
    DEBUG_ONLY(CSE cse("svd::GolubReinsch [ElementalMatrix values]"))
    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    GolubReinsch( A, s, ctrl );
}

} // namespace svd
} // namespace El

#endif // ifndef EL_SVD_GOLUBREINSCH_HPP
