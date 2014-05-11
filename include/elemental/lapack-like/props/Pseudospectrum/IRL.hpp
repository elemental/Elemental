/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_IRL_HPP
#define ELEM_PSEUDOSPECTRUM_IRL_HPP

#include "./Lanczos.hpp"

namespace elem {
namespace pspec {

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real>>& HDiagList,
  const std::vector<std::vector<Real>>& HSubdiagList,
  const Matrix<Int>& activeConverged,
  Matrix<Real>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    const Real normCap = NormCap<Real>();
    const Int numShifts = activeEsts.Height();
    if( numShifts == 0 )
        return;
    const Int basisSize = HDiagList[0].size();
    std::vector<Real> HDiag, HSubdiag, w(basisSize);
    for( Int j=0; j<numShifts; ++j )
    {
        HDiag = HDiagList[j];
        HSubdiag = HSubdiagList[j];
        if( !activeConverged.Get(j,0) )
        {
            if( !HasNan(HDiag) && !HasNan(HSubdiag) )
            {
                lapack::SymmetricTridiagEig
                ( basisSize, HDiag.data(), HSubdiag.data(), w.data(),
                  basisSize-1, basisSize-1 );
                const Real est = Sqrt(w[0]);
                activeEsts.Set( j, 0, Min(est,normCap) );
            }
            else
               activeEsts.Set( j, 0, normCap );
        }
    }
}

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real>>& HDiagList,
  const std::vector<std::vector<Real>>& HSubdiagList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
        DistMatrix<Real,MR,STAR>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    ComputeNewEstimates
    ( HDiagList, HSubdiagList, activeConverged.LockedMatrix(), 
      activeEsts.Matrix() );
}

template<typename Real>
inline void
Restart
( const std::vector<std::vector<Real>>& HDiagList,
  const std::vector<std::vector<Real>>& HSubdiagList,
  const Matrix<Int>& activeConverged,
  std::vector<Matrix<Complex<Real>>>& VList )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Restart"))
    const Int n = VList[0].Height();
    const Int numShifts = VList[0].Width();
    if( numShifts == 0 )
        return;
    const Int basisSize = HDiagList[0].size();
    std::vector<Real> HDiag, HSubdiag, w(basisSize);
    Matrix<Real> q(basisSize,1);
    Matrix<Complex<Real>> u(n,1);
    for( Int j=0; j<numShifts; ++j )
    {
        HDiag = HDiagList[j];
        HSubdiag = HSubdiagList[j];

        if( !activeConverged.Get(j,0) )
        {
            if( !HasNan(HDiag) && !HasNan(HSubdiag) )
            {
                lapack::SymmetricTridiagEig
                ( basisSize, HDiag.data(), HSubdiag.data(), w.data(), 
                  q.Buffer(), basisSize, basisSize-1, basisSize-1 );

                Zeros( u, n, 1 );
                for( Int k=0; k<basisSize; ++k )
                {
                    const Matrix<Complex<Real>>& V = VList[k];
                    auto v = LockedView( V, 0, j, n, 1 ); 
                    Axpy( q.Get(k,0), v, u );
                }
                Matrix<Complex<Real>>& V = VList[0];
                auto v = View( V, 0, j, n, 1 );
                v = u;
            }
        }
    }
}

template<typename Real>
inline void
Restart
( const std::vector<std::vector<Real>>& HDiagList,
  const std::vector<std::vector<Real>>& HSubdiagList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
  std::vector<DistMatrix<Complex<Real>>>& VList )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Restart"))
    const Int basisSize = HDiagList[0].size();
    std::vector<Matrix<Complex<Real>>> VLocList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        VLocList[j] = View( VList[j].Matrix() );
    Restart
    ( HDiagList, HSubdiagList, activeConverged.LockedMatrix(), VLocList );
}

template<typename Real>
inline Matrix<Int>
IRL
( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::IRL"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();

    const Int maxIts = psCtrl.maxIts;
    const Int basisSize = psCtrl.basisSize;
    const bool reorthog = psCtrl.reorthog;
    const bool deflate = psCtrl.deflate;
    const bool progress = psCtrl.progress;

    // Keep track of the number of iterations per shift
    Matrix<Int> itCounts;
    Ones( itCounts, numShifts, 1 );

    // Keep track of the pivoting history if deflation is requested
    Matrix<Int> preimage;
    Matrix<C> pivShifts( shifts );
    if( deflate )
    {
        preimage.Resize( numShifts, 1 );
        for( Int j=0; j<numShifts; ++j )
            preimage.Set( j, 0, j );
    }

    // The Hessenberg algorithm currently needs explicit access to the adjoint
    Matrix<C> UAdj, activeShiftsConj;
    if( !psCtrl.schur )
        Adjoint( U, UAdj );

    // Simultaneously run IRL for different shifts
    std::vector<Matrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        Zeros( VList[j], n, numShifts );
    Gaussian( VList[0], n, numShifts );
    std::vector<std::vector<Real>> 
        HDiagList(numShifts), HSubdiagList(numShifts);
    std::vector<Real> realComponents;
    std::vector<Complex<Real>> components;

    Matrix<Int> activeConverged;
    Zeros( activeConverged, numShifts, 1 );

    psCtrl.snapCtrl.ResetCounts();

    Timer timer, subtimer;
    Int numIts=0, numDone=0;
    Matrix<Real> estimates(numShifts,1);
    Zeros( estimates, numShifts, 1 );
    Matrix<Real> lastActiveEsts;
    Matrix<Int> activePreimage;
    while( true )
    {
        const Int numActive = ( deflate ? numShifts-numDone : numShifts );
        auto activeShifts = View( pivShifts, 0, 0, numActive, 1 );
        auto activeEsts = View( estimates, 0, 0, numActive, 1 );
        auto activeItCounts = View( itCounts, 0, 0, numActive, 1 );
        for( Int j=0; j<basisSize+1; ++j )
            activeVList[j] = View( VList[j], 0, 0, n, numActive ); 
        if( deflate )
        {
            activePreimage = View( preimage, 0, 0, numActive, 1 );
            Zeros( activeConverged, numActive, 1 );
        }

        // Reset the Rayleigh quotients
        for( Int j=0; j<numActive; ++j )
        {
            HDiagList[j].resize(0);
            HDiagList[j].reserve(basisSize);
            HSubdiagList[j].resize(0);
            HSubdiagList[j].reserve(basisSize);
        }

        if( progress )
            timer.Start();
        Matrix<Real> colNorms;
        ColumnNorms( activeVList[0], colNorms );
        InvBetaScale( colNorms, activeVList[0] );
        for( Int j=0; j<basisSize; ++j )
        {
            lastActiveEsts = activeEsts;
            activeVList[j+1] = activeVList[j];
            if( psCtrl.schur )
            {
                if( progress )
                    subtimer.Start();
                MultiShiftTrsm
                ( LEFT, UPPER, NORMAL, 
                  C(1), U, activeShifts, activeVList[j+1] );
                MultiShiftTrsm
                ( LEFT, UPPER, ADJOINT, 
                  C(1), U, activeShifts, activeVList[j+1] );
                if( progress )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops = 
                        (8.*n*n*numActiveShifts)/(msTime*1.e9);
                    std::cout << "  MultiShiftTrsm's: " << msTime 
                              << " seconds, " << gflops << " GFlops" 
                              << std::endl;
                }
            }
            else
            {
                if( progress )
                    subtimer.Start();
                Conjugate( activeShifts, activeShiftsConj );
                MultiShiftHessSolve
                ( UPPER, NORMAL, 
                  C(1), U, activeShifts, activeVList[j+1] );
                MultiShiftHessSolve
                ( LOWER, NORMAL, 
                  C(1), UAdj, activeShiftsConj, activeVList[j+1] );
                if( progress )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops = 
                        (32.*n*n*numActiveShifts)/(msTime*1.e9);
                    std::cout << "  MultiShiftHessSolve's: " << msTime
                              << " seconds, " << gflops << " GFlops" 
                              << std::endl;
                }
            }

            // Orthogonalize with respect to the old iterate
            if( j > 0 )
            {
                ExtractList( HSubdiagList, realComponents, j-1 );
                ColumnSubtractions
                ( realComponents, activeVList[j-1], activeVList[j+1] );
            }

            // Orthogonalize with respect to the last iterate
            InnerProducts( activeVList[j], activeVList[j+1], realComponents );
            PushBackList( HDiagList, realComponents );
            ColumnSubtractions
            ( realComponents, activeVList[j], activeVList[j+1] );

            if( reorthog )
            {
                // Explicitly (re)orthogonalize against all previous vectors
                for( Int i=0; i<j; ++i )
                {
                    InnerProducts
                    ( activeVList[i], activeVList[j+1], components );
                    ColumnSubtractions
                    ( components, activeVList[i], activeVList[j+1] );
                }
            }

            // Compute the norm of what is left
            ColumnNorms( activeVList[j+1], realComponents );
            PushBackList( HSubdiagList, realComponents );

            // TODO: Handle lucky breakdowns
            InvBetaScale( realComponents, activeVList[j+1] );

            ComputeNewEstimates
            ( HDiagList, HSubdiagList, activeConverged, activeEsts );
            // We will have the same estimate two iterations in a row when
            // restarting
            if( j != 0 ) 
                activeConverged =
                    FindConverged
                    ( lastActiveEsts, activeEsts, activeItCounts, psCtrl.tol );

            psCtrl.snapCtrl.Iterate();
        }
        if( progress )
            subtimer.Start();
        Restart
        ( HDiagList, HSubdiagList, activeVList, activeConverged, activeEsts );
        if( progress )
            std::cout << "IRL restart: " << subtimer.Stop()
                      << " seconds" << std::endl;

        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        numIts += basisSize;
        if( progress )
        {
            const double iterTime = timer.Stop();
            std::cout << "iteration " << numIts << ": " << iterTime
                      << " seconds, " << numDone << " of " << numShifts
                      << " converged" << std::endl;
        }
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate && numActiveDone != 0 )
        {
            Deflate
            ( activeShifts, activePreimage, activeVList[0], activeEsts, 
              activeConverged, activeItCounts, progress );
            lastActiveEsts = activeEsts;
        }

        // Save snapshots of the estimates at the requested rate
        Snapshot
        ( preimage, estimates, itCounts, numIts, deflate, psCtrl.snapCtrl );
    } 

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );
    FinalSnapshot( invNorms, itCounts, psCtrl.snapCtrl );

    return itCounts;
}

template<typename Real>
inline DistMatrix<Int,VR,STAR>
IRL
( const DistMatrix<Complex<Real>        >& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::IRL"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();
    const Grid& g = U.Grid();

    const Int maxIts = psCtrl.maxIts;
    const Int basisSize = psCtrl.basisSize;
    const bool reorthog = psCtrl.reorthog;
    const bool deflate = psCtrl.deflate;
    const bool progress = psCtrl.progress;

    // Keep track of the number of iterations per shift
    DistMatrix<Int,VR,STAR> itCounts(g);
    Ones( itCounts, numShifts, 1 );

    // Keep track of the pivoting history if deflation is requested
    DistMatrix<Int,VR,STAR> preimage(g);
    DistMatrix<C,  VR,STAR> pivShifts( shifts );
    if( deflate )
    {
        preimage.AlignWith( shifts );
        preimage.Resize( numShifts, 1 );
        const Int numLocShifts = preimage.LocalHeight();
        for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
        {
            const Int i = preimage.GlobalRow(iLoc);
            preimage.SetLocal( iLoc, 0, i );
        }
    }

    // The Hessenberg algorithm currently requires explicit adjoint access
    DistMatrix<C,VC,STAR> U_VC_STAR(g), UAdj_VC_STAR(g);
    DistMatrix<C,VR,STAR> activeShiftsConj(g);
    DistMatrix<C,STAR,VR> activeV_STAR_VR(g);
    if( !psCtrl.schur )
    {
        U_VC_STAR = U;
        Adjoint( U, UAdj_VC_STAR );
    }

    // Simultaneously run IRL for different shifts
    std::vector<DistMatrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
    {
        VList[j].SetGrid( g );
        Zeros( VList[j], n, numShifts );
    }
    Gaussian( VList[0], n, numShifts );
    const Int numMRShifts = VList[0].LocalWidth();
    std::vector<std::vector<Real>> HDiagList(numMRShifts), 
                                   HSubdiagList(numMRShifts);
    std::vector<Real> realComponents;
    std::vector<Complex<Real>> components;

    DistMatrix<Int,MR,STAR> activeConverged(g);
    Zeros( activeConverged, numShifts, 1 );

    psCtrl.snapCtrl.ResetCounts();

    Timer timer, subtimer;
    Int numIts=0, numDone=0;
    DistMatrix<Real,MR,STAR> estimates(g), lastActiveEsts(g);
    estimates.AlignWith( shifts );
    Zeros( estimates, numShifts, 1 );
    DistMatrix<Int,VR,STAR> activePreimage(g);
    while( true )
    {
        const Int numActive = ( deflate ? numShifts-numDone : numShifts );
        auto activeShifts = View( pivShifts, 0, 0, numActive, 1 );
        auto activeEsts = View( estimates, 0, 0, numActive, 1 );
        auto activeItCounts = View( itCounts, 0, 0, numActive, 1 );
        for( Int j=0; j<basisSize+1; ++j )
            activeVList[j] = View( VList[j], 0, 0, n, numActive ); 
        if( deflate )
        {
            activePreimage = View( preimage, 0, 0, numActive, 1 );
            Zeros( activeConverged, numActive, 1 );
        }

        // Reset the Rayleigh quotients
        const Int numActiveMR = estimates.LocalHeight();
        for( Int jLoc=0; jLoc<numActiveMR; ++jLoc )
        {
            HDiagList[jLoc].resize(0);
            HDiagList[jLoc].reserve(basisSize);
            HSubdiagList[jLoc].resize(0);
            HSubdiagList[jLoc].reserve(basisSize);
        }

        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                timer.Start();
        }
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnNorms( activeVList[0], colNorms );
        InvBetaScale( colNorms, activeVList[0] );
        for( Int j=0; j<basisSize; ++j )
        {
            lastActiveEsts = activeEsts;
            activeVList[j+1] = activeVList[j];
            if( psCtrl.schur )
            {
                if( progress )
                { 
                    mpi::Barrier( g.Comm() );
                    if( g.Rank() == 0 )
                        subtimer.Start();
                }
                MultiShiftTrsm
                ( LEFT, UPPER, NORMAL, 
                  C(1), U, activeShifts, activeVList[j+1] );
                MultiShiftTrsm
                ( LEFT, UPPER, ADJOINT, 
                  C(1), U, activeShifts, activeVList[j+1] );
                if( progress )
                {
                    mpi::Barrier( g.Comm() );
                    if( g.Rank() == 0 )
                    {
                        const double msTime = subtimer.Stop();
                        const Int numActiveShifts = activeShifts.Height();
                        const double gflops = 
                            (8.*n*n*numActiveShifts)/(msTime*1.e9);
                        std::cout << "  MultiShiftTrsm's: " << msTime 
                                  << " seconds, " << gflops << " GFlops" 
                                  << std::endl;
                    }
                }
            }
            else
            {
                if( progress )
                {
                    mpi::Barrier( g.Comm() );
                    if( g.Rank() == 0 )
                        subtimer.Start();
                }
                // NOTE: This redistribution sequence might not be necessary
                activeV_STAR_VR = activeVList[j+1];
                Conjugate( activeShifts, activeShiftsConj );
                MultiShiftHessSolve
                ( UPPER, NORMAL, C(1), U_VC_STAR, activeShifts,
                  activeV_STAR_VR );
                MultiShiftHessSolve
                ( LOWER, NORMAL, C(1), UAdj_VC_STAR, activeShiftsConj,
                  activeV_STAR_VR );
                activeVList[j+1] = activeV_STAR_VR;
                if( progress )
                {
                    mpi::Barrier( g.Comm() );
                    if( g.Rank() == 0 )
                    {
                        const double msTime = subtimer.Stop();
                        const Int numActiveShifts = activeShifts.Height();
                        const double gflops =
                            (32.*n*n*numActiveShifts)/(msTime*1.e9);
                        std::cout << "  MultiShiftHessSolve's: " << msTime
                                  << " seconds, " << gflops << " GFlops"
                                  << std::endl;
                    }
                }
            }

            // Orthogonalize with respect to the old iterate
            if( j > 0 )
            {
                ExtractList( HSubdiagList, realComponents, j-1 );
                ColumnSubtractions
                ( realComponents, activeVList[j-1], activeVList[j+1] );
            }

            // Orthogonalize with respect to the last iterate
            InnerProducts( activeVList[j], activeVList[j+1], realComponents );
            PushBackList( HDiagList, realComponents );
            ColumnSubtractions
            ( realComponents, activeVList[j], activeVList[j+1] );

            if( reorthog )
            {
                // Explicitly (re)orthogonalize against all previous vectors
                for( Int i=0; i<j; ++i )
                {
                    InnerProducts
                    ( activeVList[i], activeVList[j+1], components );
                    ColumnSubtractions
                    ( components, activeVList[i], activeVList[j+1] );
                }
            }

            // Compute the norm of what is left
            ColumnNorms( activeVList[j+1], realComponents );
            PushBackList( HSubdiagList, realComponents );

            // TODO: Handle lucky breakdowns
            InvBetaScale( realComponents, activeVList[j+1] );

            ComputeNewEstimates
            ( HDiagList, HSubdiagList, activeConverged, activeEsts );
            // We will have the same estimate two iterations in a row when
            // restarting
            if( j != 0 ) 
                activeConverged =
                    FindConverged
                    ( lastActiveEsts, activeEsts, activeItCounts, psCtrl.tol );

            psCtrl.snapCtrl.Iterate();
        }
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                subtimer.Start();
        }
        Restart( HDiagList, HSubdiagList, activeConverged, activeVList );
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                std::cout << "IRL computations: " << subtimer.Stop()
                          << " seconds" << std::endl;
        }

        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        numIts += basisSize;
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
            {
                const double iterTime = timer.Stop();
                std::cout << "iteration " << numIts << ": " << iterTime
                          << " seconds, " << numDone << " of " << numShifts
                          << " converged" << std::endl;
            }
        }
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate && numActiveDone != 0 )
        {
            Deflate
            ( activeShifts, activePreimage, activeVList[0], activeEsts, 
              activeConverged, activeItCounts, progress );
            lastActiveEsts = activeEsts;
        }

        // Save snapshots of the estimates at the requested rate
        Snapshot
        ( preimage, estimates, itCounts, numIts, deflate, psCtrl.snapCtrl );
    } 

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );
    FinalSnapshot( invNorms, itCounts, psCtrl.snapCtrl );

    return itCounts;
}

} // namespace pspec
} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_IRL_HPP
