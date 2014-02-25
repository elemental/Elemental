/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_KRYLOVSPECTRAL_HPP
#define ELEM_PSEUDOSPECTRUM_KRYLOVSPECTRAL_HPP

#include "./Lanczos.hpp"

namespace elem {
namespace pspec {

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real> >& HDiagList,
  const std::vector<std::vector<Real> >& HSubdiagList,
  const Matrix<Int>& activeConverged,
  Matrix<Real>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    const Real normCap = NormCap<Real>();
    const Int numShifts = activeEsts.Height();
    if( numShifts == 0 )
        return;
    const Int krylovSize = HDiagList[0].size();
    std::vector<Real> HDiag, HSubdiag, w(krylovSize);
    for( Int j=0; j<numShifts; ++j )
    {
        HDiag = HDiagList[j];
        HSubdiag = HSubdiagList[j];
        if( !activeConverged.Get(j,0) )
        {
            if( !HasNan( HDiag, HSubdiag ) )
            {
                lapack::SymmetricTridiagEig
                ( 'N', 'I', krylovSize, HDiag.data(), HSubdiag.data(), 0, 0,
                  krylovSize, krylovSize, 0, w.data(), 0, 1 );
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
( const std::vector<std::vector<Real> >& HDiagList,
  const std::vector<std::vector<Real> >& HSubdiagList,
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
( const std::vector<std::vector<Real> >& HDiagList,
  const std::vector<std::vector<Real> >& HSubdiagList,
  const Matrix<Int>& activeConverged,
  std::vector<Matrix<Complex<Real> > >& VList )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Restart"))
    const Int n = VList[0].Height();
    const Int numShifts = VList[0].Width();
    if( numShifts == 0 )
        return;
    const Int krylovSize = HDiagList[0].size();
    std::vector<Real> HDiag, HSubdiag, w(krylovSize);
    Matrix<Real> q(krylovSize,1);
    Matrix<Complex<Real>> u(n,1);
    for( Int j=0; j<numShifts; ++j )
    {
        HDiag = HDiagList[j];
        HSubdiag = HSubdiagList[j];

        if( !activeConverged.Get(j,0) )
        {
            if( !HasNan( HDiag, HSubdiag ) )
            {
                lapack::SymmetricTridiagEig
                ( 'V', 'I', krylovSize, HDiag.data(), HSubdiag.data(), 
                  0, 0, krylovSize, krylovSize, 0, w.data(), 
                  q.Buffer(), krylovSize );

                Zeros( u, n, 1 );
                for( Int k=0; k<krylovSize; ++k )
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
( const std::vector<std::vector<Real> >& HDiagList,
  const std::vector<std::vector<Real> >& HSubdiagList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
  std::vector<DistMatrix<Complex<Real> > >& VList )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Restart"))
    const Int krylovSize = HDiagList[0].size();
    std::vector<Matrix<Complex<Real>>> VLocList(krylovSize+1);
    for( Int j=0; j<krylovSize+1; ++j )
        VLocList[j] = View( VList[j].Matrix() );
    Restart
    ( HDiagList, HSubdiagList, activeConverged.LockedMatrix(), VLocList );
}

template<typename Real>
inline Matrix<Int>
TriangularKrylovSpectral
( const Matrix<Complex<Real> >& U, const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, const Int krylovSize=10, bool reorthog=true,
  bool deflate=true, Int maxIts=1000, Real tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::TriangularKrylovSpectral"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();

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

    // Simultaneously run a Krylov-spectral method for various shifts
    std::vector<Matrix<C>> VList(krylovSize+1), activeVList(krylovSize+1);
    for( Int j=0; j<krylovSize+1; ++j )
        Zeros( VList[j], n, numShifts );
    Gaussian( VList[0], n, numShifts );
    std::vector<std::vector<Real>> 
        HDiagList(numShifts), HSubdiagList(numShifts);
    std::vector<Complex<Real>> components;

    Matrix<Int> activeConverged;
    Zeros( activeConverged, numShifts, 1 );

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
        for( Int j=0; j<krylovSize+1; ++j )
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
            HDiagList[j].reserve(krylovSize);
            HSubdiagList[j].resize(0);
            HSubdiagList[j].reserve(krylovSize);
        }

        if( progress )
            timer.Start();
        Matrix<Real> colNorms;
        ColumnNorms( activeVList[0], colNorms );
        InvBetaScale( colNorms, activeVList[0] );
        for( Int j=0; j<krylovSize; ++j )
        {
            lastActiveEsts = activeEsts;
            activeVList[j+1] = activeVList[j];
            if( progress )
                subtimer.Start();
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeVList[j+1] );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeVList[j+1] );
            if( progress )
            {
                const double msTime = subtimer.Stop();
                const double gflops = (4.*n*n*numShifts)/(msTime*1.e9);
                std::cout << "  MultiShiftTrsm's: " << msTime << " seconds, "
                          << gflops << " GFlops" << std::endl;
            }
            if( j != 0 )
            {
                ColumnSubtractions
                ( HSubdiagList, activeVList[j-1], activeVList[j+1] );
            }
            InnerProducts( activeVList[j], activeVList[j+1], HDiagList );
            ColumnSubtractions( HDiagList, activeVList[j], activeVList[j+1] );
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
            ColumnNorms( activeVList[j+1], HSubdiagList );
            // TODO: Handle lucky breakdowns
            InvBetaScale( HSubdiagList, activeVList[j+1] );

            ComputeNewEstimates
            ( HDiagList, HSubdiagList, activeConverged, activeEsts );
            // We will have the same estimate two iterations in a row when
            // restarting
            if( j != 0 ) 
                activeConverged =
                    FindConverged
                    ( lastActiveEsts, activeEsts, activeItCounts, tol );
        }
        if( progress )
            subtimer.Start();
        Restart
        ( HDiagList, HSubdiagList, activeVList, activeConverged, activeEsts );
        if( progress )
            std::cout << "Krylov-spectral contraction: " << subtimer.Stop()
                      << " seconds" << std::endl;

        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        numIts += krylovSize;
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
    } 
    if( numDone != numShifts )
        RuntimeError("Two-norm estimates did not converge in time");

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );

    return itCounts;
}

template<typename Real>
inline DistMatrix<Int,VR,STAR>
TriangularKrylovSpectral
( const DistMatrix<Complex<Real>        >& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
        Int krylovSize=10, bool reorthog=true,
  bool deflate=true, Int maxIts=1000, Real tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::TriangularKrylovSpectral"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();
    const Grid& g = U.Grid();
    if( deflate && U.Grid().Rank() == 0 ) 
        std::cerr << "WARNING: Deflation swaps not yet optimized!" << std::endl;

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
        for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        {
            const Int j = preimage.ColShift() + jLoc*preimage.ColStride();
            preimage.SetLocal( jLoc, 0, j );
        }
    }

    // Simultaneously run a Krylov-spectral method for various shifts
    std::vector<DistMatrix<C>> VList(krylovSize+1), activeVList(krylovSize+1);
    for( Int j=0; j<krylovSize+1; ++j )
    {
        VList[j].SetGrid( U.Grid() );
        Zeros( VList[j], n, numShifts );
    }
    Gaussian( VList[0], n, numShifts );
    const Int numMRShifts = VList[0].LocalWidth();
    std::vector<std::vector<Real>> HDiagList(numMRShifts), 
                                   HSubdiagList(numMRShifts);
    std::vector<Complex<Real>> components;

    DistMatrix<Int,MR,STAR> activeConverged(g);
    Zeros( activeConverged, numShifts, 1 );

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
        for( Int j=0; j<krylovSize+1; ++j )
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
            HDiagList[jLoc].reserve(krylovSize);
            HSubdiagList[jLoc].resize(0);
            HSubdiagList[jLoc].reserve(krylovSize);
        }

        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                timer.Start();
        }
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnNorms( activeVList[0], colNorms );
        InvBetaScale( colNorms, activeVList[0] );
        for( Int j=0; j<krylovSize; ++j )
        {
            lastActiveEsts = activeEsts;
            activeVList[j+1] = activeVList[j];
            if( progress )
            { 
                mpi::Barrier( U.Grid().Comm() );
                if( U.Grid().Rank() == 0 )
                    subtimer.Start();
            }
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeVList[j+1] );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeVList[j+1] );
            if( progress )
            {
                mpi::Barrier( U.Grid().Comm() );
                if( U.Grid().Rank() == 0 )
                {
                    const double msTime = subtimer.Stop();
                    const double gflops = (4.*n*n*numShifts)/(msTime*1.e9);
                    std::cout << "  MultiShiftTrsm's: " << msTime 
                              << " seconds, " << gflops << " GFlops" 
                              << std::endl;
                }
            }
            if( j != 0 )
            {
                ColumnSubtractions
                ( HSubdiagList, activeVList[j-1], activeVList[j+1] );
            }
            InnerProducts( activeVList[j], activeVList[j+1], HDiagList );
            ColumnSubtractions( HDiagList, activeVList[j], activeVList[j+1] );
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
            ColumnNorms( activeVList[j+1], HSubdiagList );
            // TODO: Handle lucky breakdowns
            InvBetaScale( HSubdiagList, activeVList[j+1] );

            ComputeNewEstimates
            ( HDiagList, HSubdiagList, activeConverged, activeEsts );
            // We will have the same estimate two iterations in a row when
            // restarting
            if( j != 0 ) 
                activeConverged =
                    FindConverged
                    ( lastActiveEsts, activeEsts, activeItCounts, tol );
        }
        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                subtimer.Start();
        }
        Restart( HDiagList, HSubdiagList, activeConverged, activeVList );
        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                std::cout << "Krylov-spectral computations: " << subtimer.Stop()
                          << " seconds" << std::endl;
        }

        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        numIts += krylovSize;
        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
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
    } 
    if( numDone != numShifts )
        RuntimeError("Two-norm estimates did not converge in time");

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );

    return itCounts;
}

} // namespace pspec
} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_KRYLOVSPECTRAL_HPP
