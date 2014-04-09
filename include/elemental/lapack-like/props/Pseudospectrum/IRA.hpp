/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_IRA_HPP
#define ELEM_PSEUDOSPECTRUM_IRA_HPP

#include "./Lanczos.hpp"

namespace elem {
namespace pspec {

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<Matrix<Complex<Real> > >& HList,
  const Matrix<Int>& activeConverged,
  Matrix<Real>& activeEsts,
  Int n )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    const Real normCap = NormCap<Real>();
    const Int numShifts = activeEsts.Height();
    if( numShifts == 0 )
        return;
    Matrix<Complex<Real>> H, HTL;
    Matrix<Complex<Real>> w(n,1);
    for( Int j=0; j<numShifts; ++j )
    {
        H = HList[j];
        HTL = View( H, 0, 0, n, n );
        if( !activeConverged.Get(j,0) )
        {
            if( !HasNan(HTL) )
            {
                lapack::HessenbergEig
                ( n, HTL.Buffer(), HTL.LDim(), w.Buffer() );
                Real estSquared=0;
                for( Int k=0; k<n; ++k )
                    if( w.GetRealPart(k,0) > estSquared )
                        estSquared = w.GetRealPart(k,0);
                activeEsts.Set( j, 0, Min(Sqrt(estSquared),normCap) );
            }
            else
               activeEsts.Set( j, 0, normCap );
        }
    }
}

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<Matrix<Complex<Real> > >& HList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
        DistMatrix<Real,MR,STAR>& activeEsts,
        Int n )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    ComputeNewEstimates
    ( HList, activeConverged.LockedMatrix(), activeEsts.Matrix(), n );
}

template<typename Real>
inline void
Restart
( const std::vector<Matrix<Complex<Real> > >& HList,
  const Matrix<Int>& activeConverged,
  std::vector<Matrix<Complex<Real> > >& VList )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Restart"))
    const Int n = VList[0].Height();
    const Int numShifts = VList[0].Width();
    if( numShifts == 0 )
        return;
    const Int basisSize = HList[0].Width();
    Matrix<Complex<Real>> H, HTL, Q(basisSize,basisSize);
    Matrix<Complex<Real>> w(basisSize,1), u(n,1);
    for( Int j=0; j<numShifts; ++j )
    {
        H = HList[j];
        HTL = View( H, 0, 0, basisSize, basisSize );

        if( !activeConverged.Get(j,0) )
        {
            if( !HasNan(HTL) )
            {
                // TODO: Switch to lapack::HessenbergEig
                lapack::Eig
                ( basisSize, HTL.Buffer(), HTL.LDim(), w.Buffer(), 
                  Q.Buffer(), Q.LDim() );

                Real maxReal=0;
                Int maxIdx=0;
                for( Int k=0; k<basisSize; ++k )
                {
                    if( w.GetRealPart(k,0) > maxReal )
                    {
                        maxReal = w.GetRealPart(k,0);
                        maxIdx = k;
                    }
                }

                Zeros( u, n, 1 );
                for( Int k=0; k<basisSize; ++k )
                {
                    const Matrix<Complex<Real>>& V = VList[k];
                    auto v = LockedView( V, 0, j, n, 1 ); 
                    Axpy( Q.Get(k,maxIdx), v, u );
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
( const std::vector<Matrix<Complex<Real> > >& HList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
  std::vector<DistMatrix<Complex<Real> > >& VList )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Restart"))
    const Int basisSize = HList[0].Width();
    std::vector<Matrix<Complex<Real>>> VLocList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        VLocList[j] = View( VList[j].Matrix() );
    Restart( HList, activeConverged.LockedMatrix(), VLocList );
}

template<typename Real>
inline void
Deflate
( std::vector<Matrix<Complex<Real> > >& HList,
  Matrix<Complex<Real> >& activeShifts,
  Matrix<Int           >& activePreimage,
  Matrix<Complex<Real> >& activeXOld,
  Matrix<Complex<Real> >& activeX,
  Matrix<Real          >& activeEsts,
  Matrix<Int           >& activeConverged,
  Matrix<Int           >& activeItCounts,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Deflate"))
    Timer timer;
    if( progress )
        timer.Start();
    const Int numActive = activeX.Width();
    Int swapTo = numActive-1;
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( activeConverged.Get(swapFrom,0) )
        {
            if( swapTo != swapFrom )
            {
                std::swap( HList[swapFrom], HList[swapTo] );
                RowSwap( activeShifts, swapFrom, swapTo );
                RowSwap( activePreimage, swapFrom, swapTo );
                RowSwap( activeEsts, swapFrom, swapTo );
                RowSwap( activeItCounts, swapFrom, swapTo );
                ColumnSwap( activeXOld, swapFrom, swapTo );
                ColumnSwap( activeX,    swapFrom, swapTo );
            }
            --swapTo;
        }
    }
    if( progress )
        std::cout << "Deflation took " << timer.Stop() << " seconds"
                  << std::endl;
}

template<typename Real>
inline void
Deflate
( std::vector<Matrix<Complex<Real> > >& HList,
  DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
  DistMatrix<Complex<Real>        >& activeXOld,
  DistMatrix<Complex<Real>        >& activeX,
  DistMatrix<Real,         MR,STAR>& activeEsts,
  DistMatrix<Int,          MR,STAR>& activeConverged,
  DistMatrix<Int,          VR,STAR>& activeItCounts,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Deflate"))
    Timer timer;
    if( progress && activeShifts.Grid().Rank() == 0 )
        timer.Start();
    const Int numActive = activeX.Width();
    Int swapTo = numActive-1;

    DistMatrix<Complex<Real>,STAR,STAR> shiftsCopy( activeShifts );
    DistMatrix<Int,STAR,STAR> preimageCopy( activePreimage );
    DistMatrix<Real,STAR,STAR> estimatesCopy( activeEsts );
    DistMatrix<Int, STAR,STAR> itCountsCopy( activeItCounts );
    DistMatrix<Int, STAR,STAR> convergedCopy( activeConverged );
    DistMatrix<Complex<Real>,VC,STAR> XOldCopy( activeXOld ), XCopy( activeX );

    const Int n = ( activeX.LocalWidth()>0 ? HList[0].Height() : 0 );
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( convergedCopy.Get(swapFrom,0) )
        {
            if( swapTo != swapFrom )
            {
                // TODO: Avoid this large latency penalty
                if( activeX.IsLocalCol(swapFrom) &&
                    activeX.IsLocalCol(swapTo) )
                {
                    const Int localFrom = activeX.LocalCol(swapFrom);
                    const Int localTo = activeX.LocalCol(swapTo);
                    std::swap( HList[localFrom], HList[localTo] );
                }
                else if( activeX.IsLocalCol(swapFrom) )
                {
                    const Int localFrom = activeX.LocalCol(swapFrom);
                    const Int partner = activeX.ColOwner(swapTo);
                    DEBUG_ONLY(
                        if( HList[localFrom].LDim() != n )
                            LogicError("Leading dimension was incorrect");
                    )
                    mpi::TaggedSendRecv
                    ( HList[localFrom].Buffer(), n*n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                }
                else if( activeX.IsLocalCol(swapTo) )
                {
                    const Int localTo = activeX.LocalCol(swapTo);
                    DEBUG_ONLY(
                        if( HList[localTo].LDim() != n )
                            LogicError("Leading dimension was incorrect");
                    )
                    const Int partner = activeX.ColOwner(swapFrom);
                    mpi::TaggedSendRecv
                    ( HList[localTo].Buffer(), n*n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                }

                RowSwap( shiftsCopy, swapFrom, swapTo );
                RowSwap( preimageCopy, swapFrom, swapTo );
                RowSwap( estimatesCopy, swapFrom, swapTo );
                RowSwap( itCountsCopy, swapFrom, swapTo );
                ColumnSwap( XOldCopy, swapFrom, swapTo );
                ColumnSwap( XCopy,    swapFrom, swapTo );
            }
            --swapTo;
        }
    }

    activeShifts   = shiftsCopy;
    activePreimage = preimageCopy;
    activeEsts     = estimatesCopy;
    activeItCounts = itCountsCopy;
    activeXOld     = XOldCopy;
    activeX        = XCopy;

    if( progress )
    {
        mpi::Barrier( activeShifts.Grid().Comm() );
        if( activeShifts.Grid().Rank() == 0 )
            std::cout << "Deflation took " << timer.Stop() << " seconds"
                      << std::endl;
    }
}

template<typename Real>
inline Matrix<Int>
TriangularIRA
( const Matrix<Complex<Real> >& U, const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, const Int basisSize=10, 
  Int maxIts=1000, Real tol=1e-6, bool progress=false, bool deflate=true,
  SnapshotCtrl snapCtrl=SnapshotCtrl() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::TriangularIRA"))
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

    // Simultaneously run IRA for different shifts
    std::vector<Matrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        Zeros( VList[j], n, numShifts );
    Gaussian( VList[0], n, numShifts );
    std::vector<Matrix<Complex<Real>>> HList(numShifts);
    std::vector<Real> realComponents;
    std::vector<Complex<Real>> components;

    Matrix<Int> activeConverged;
    Zeros( activeConverged, numShifts, 1 );

    snapCtrl.numSaveCount = 0;
    snapCtrl.imgSaveCount = 0;

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
        HList.resize( numActive );
        for( Int j=0; j<numActive; ++j )
            Zeros( HList[j], basisSize+1, basisSize );

        if( progress )
            timer.Start();
        ColumnNorms( activeVList[0], realComponents );
        InvBetaScale( realComponents, activeVList[0] );
        for( Int j=0; j<basisSize; ++j )
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
                const Int numActiveShifts = activeShifts.Height();
                const double gflops = (8.*n*n*numActiveShifts)/(msTime*1.e9);
                std::cout << "  MultiShiftTrsm's: " << msTime << " seconds, "
                          << gflops << " GFlops" << std::endl;
            }

            // Orthogonalize with respect to the old iterate
            if( j > 0 )
            {
                ExtractList( HList, components, j, j-1 );
                // TODO: Conjugate components?
                PlaceList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVList[j-1], activeVList[j+1] );
            }

            // Orthogonalize with respect to the last iterate
            InnerProducts( activeVList[j], activeVList[j+1], components );
            PlaceList( HList, components, j, j );
            ColumnSubtractions
            ( components, activeVList[j], activeVList[j+1] );

            // Explicitly (re)orthogonalize against all previous vectors
            for( Int i=0; i<j-1; ++i )
            {
                InnerProducts
                ( activeVList[i], activeVList[j+1], components );
                PlaceList( HList, components, i, j );
                ColumnSubtractions
                ( components, activeVList[i], activeVList[j+1] );
            }
            if( j > 0 )
            {
                InnerProducts
                ( activeVList[j-1], activeVList[j+1], components );
                UpdateList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVList[j-1], activeVList[j+1] );
            }

            // Compute the norm of what is left
            ColumnNorms( activeVList[j+1], realComponents );
            PlaceList( HList, realComponents, j+1, j );

            // TODO: Handle lucky breakdowns
            InvBetaScale( realComponents, activeVList[j+1] );

            ComputeNewEstimates( HList, activeConverged, activeEsts, j+1 );
            // We will have the same estimate two iterations in a row when
            // restarting
            if( j != 0 ) 
                activeConverged =
                    FindConverged
                    ( lastActiveEsts, activeEsts, activeItCounts, tol );

            if( snapCtrl.numFreq > 0 )
                ++snapCtrl.numSaveCount;
            if( snapCtrl.imgFreq > 0 )
                ++snapCtrl.imgSaveCount;
        }
        if( progress )
            subtimer.Start();
        Restart( HList, activeVList, activeConverged, activeEsts );
        if( progress )
            std::cout << "IRA restart: " << subtimer.Stop()
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
        Snapshot( estimates, preimage, numIts, deflate, snapCtrl );
    } 

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );

    return itCounts;
}

template<typename Real>
inline Matrix<Int>
HessenbergIRA
( const Matrix<Complex<Real> >& H, const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, const Int basisSize=10, 
  Int maxIts=1000, Real tol=1e-6, bool progress=false, bool deflate=true,
  SnapshotCtrl snapCtrl=SnapshotCtrl() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::HessenbergIRA"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = H.Height();
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

    Matrix<C> HAdj;
    Adjoint( H, HAdj );
    Matrix<C> activeShiftsConj;

    // Simultaneously run IRA for different shifts
    std::vector<Matrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        Zeros( VList[j], n, numShifts );
    Gaussian( VList[0], n, numShifts );
    std::vector<Matrix<Complex<Real>>> HList(numShifts);
    std::vector<Real> realComponents;
    std::vector<Complex<Real>> components;

    Matrix<Int> activeConverged;
    Zeros( activeConverged, numShifts, 1 );

    snapCtrl.numSaveCount = 0;
    snapCtrl.imgSaveCount = 0;

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
        HList.resize( numActive );
        for( Int j=0; j<numActive; ++j )
            Zeros( HList[j], basisSize+1, basisSize );

        if( progress )
            timer.Start();
        ColumnNorms( activeVList[0], realComponents );
        InvBetaScale( realComponents, activeVList[0] );
        for( Int j=0; j<basisSize; ++j )
        {
            lastActiveEsts = activeEsts;
            activeVList[j+1] = activeVList[j];
            if( progress )
                subtimer.Start();
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), H, activeShifts, activeVList[j+1] );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), HAdj, activeShiftsConj, activeVList[j+1] );
            if( progress )
            {
                const double msTime = subtimer.Stop();
                const Int numActiveShifts = activeShifts.Height();
                const double gflops = (32.*n*n*numActiveShifts)/(msTime*1.e9);
                std::cout << "  MultiShiftHessSolve's: " << msTime 
                          << " seconds, " << gflops << " GFlops" << std::endl;
            }

            // Orthogonalize with respect to the old iterate
            if( j > 0 )
            {
                ExtractList( HList, components, j, j-1 );
                // TODO: Conjugate components?
                PlaceList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVList[j-1], activeVList[j+1] );
            }

            // Orthogonalize with respect to the last iterate
            InnerProducts( activeVList[j], activeVList[j+1], components );
            PlaceList( HList, components, j, j );
            ColumnSubtractions
            ( components, activeVList[j], activeVList[j+1] );

            // Explicitly (re)orthogonalize against all previous vectors
            for( Int i=0; i<j-1; ++i )
            {
                InnerProducts
                ( activeVList[i], activeVList[j+1], components );
                PlaceList( HList, components, i, j );
                ColumnSubtractions
                ( components, activeVList[i], activeVList[j+1] );
            }
            if( j > 0 )
            {
                InnerProducts
                ( activeVList[j-1], activeVList[j+1], components );
                UpdateList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVList[j-1], activeVList[j+1] );
            }

            // Compute the norm of what is left
            ColumnNorms( activeVList[j+1], realComponents );
            PlaceList( HList, realComponents, j+1, j );

            // TODO: Handle lucky breakdowns
            InvBetaScale( realComponents, activeVList[j+1] );

            ComputeNewEstimates( HList, activeConverged, activeEsts, j+1 );
            // We will have the same estimate two iterations in a row when
            // restarting
            if( j != 0 ) 
                activeConverged =
                    FindConverged
                    ( lastActiveEsts, activeEsts, activeItCounts, tol );    

            if( snapCtrl.numFreq > 0 )
                ++snapCtrl.numSaveCount;
            if( snapCtrl.imgFreq > 0 )
                ++snapCtrl.imgSaveCount;
        }
        if( progress )
            subtimer.Start();
        Restart( HList, activeVList, activeConverged, activeEsts );
        if( progress )
            std::cout << "IRA restart: " << subtimer.Stop()
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
        Snapshot( estimates, preimage, numIts, deflate, snapCtrl );
    } 

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );

    return itCounts;
}

template<typename Real>
inline DistMatrix<Int,VR,STAR>
TriangularIRA
( const DistMatrix<Complex<Real>        >& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
        Int basisSize=10, 
  Int maxIts=1000, Real tol=1e-6, bool progress=false, bool deflate=true,
  SnapshotCtrl snapCtrl=SnapshotCtrl() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::TriangularIRA"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();
    const Grid& g = U.Grid();
    if( deflate && U.Grid().Rank() == 0 ) 
        std::cerr << "NOTE: Deflation swaps not yet optimized!" << std::endl;

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

    // Simultaneously run IRA for different shifts
    std::vector<DistMatrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
    {
        VList[j].SetGrid( U.Grid() );
        Zeros( VList[j], n, numShifts );
    }
    Gaussian( VList[0], n, numShifts );
    const Int numMRShifts = VList[0].LocalWidth();
    std::vector<Matrix<Complex<Real>>> HList(numMRShifts);
    std::vector<Real> realComponents;
    std::vector<Complex<Real>> components;

    DistMatrix<Int,MR,STAR> activeConverged(g);
    Zeros( activeConverged, numShifts, 1 );

    snapCtrl.numSaveCount = 0;
    snapCtrl.imgSaveCount = 0;

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
        HList.resize( activeEsts.LocalHeight() );
        for( Int jLoc=0; jLoc<HList.size(); ++jLoc )
            Zeros( HList[jLoc], basisSize+1, basisSize );

        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                timer.Start();
        }
        ColumnNorms( activeVList[0], realComponents );
        InvBetaScale( realComponents, activeVList[0] );
        for( Int j=0; j<basisSize; ++j )
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
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops = 
                        (8.*n*n*numActiveShifts)/(msTime*1.e9);
                    std::cout << "  MultiShiftTrsm's: " << msTime 
                              << " seconds, " << gflops << " GFlops" 
                              << std::endl;
                }
            }

            // Orthogonalize with respect to the old iterate
            if( j > 0 )
            {
                ExtractList( HList, components, j, j-1 );
                // TODO: Conjugate components?
                PlaceList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVList[j-1], activeVList[j+1] );
            }

            // Orthogonalize with respect to the last iterate
            InnerProducts( activeVList[j], activeVList[j+1], components );
            PlaceList( HList, components, j, j );
            ColumnSubtractions
            ( components, activeVList[j], activeVList[j+1] );

            // Explicitly (re)orthogonalize against all previous vectors
            for( Int i=0; i<j-1; ++i )
            {
                InnerProducts
                ( activeVList[i], activeVList[j+1], components );
                PlaceList( HList, components, i, j );
                ColumnSubtractions
                ( components, activeVList[i], activeVList[j+1] );
            }
            if( j > 0 )
            {
                InnerProducts
                ( activeVList[j-1], activeVList[j+1], components );
                UpdateList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVList[j-1], activeVList[j+1] );
            }

            // Compute the norm of what is left
            ColumnNorms( activeVList[j+1], realComponents );
            PlaceList( HList, realComponents, j+1, j );

            // TODO: Handle lucky breakdowns
            InvBetaScale( realComponents, activeVList[j+1] );

            ComputeNewEstimates( HList, activeConverged, activeEsts, j+1 );
            // We will have the same estimate two iterations in a row when
            // restarting
            if( j != 0 ) 
                activeConverged =
                    FindConverged
                    ( lastActiveEsts, activeEsts, activeItCounts, tol );

            if( snapCtrl.numFreq > 0 )
                ++snapCtrl.numSaveCount;
            if( snapCtrl.imgFreq > 0 )
                ++snapCtrl.imgSaveCount;
        }
        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                subtimer.Start();
        }
        Restart( HList, activeConverged, activeVList );
        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                std::cout << "IRA computations: " << subtimer.Stop()
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

        // Save snapshots of the estimates at the requested rate
        Snapshot( estimates, preimage, numIts, deflate, snapCtrl );
    } 

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );

    return itCounts;
}

template<typename Real>
inline DistMatrix<Int,VR,STAR>
HessenbergIRA
( const DistMatrix<Complex<Real>        >& H, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
        Int basisSize=10, 
  Int maxIts=1000, Real tol=1e-6, bool progress=false, bool deflate=true,
  SnapshotCtrl snapCtrl=SnapshotCtrl() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::HessenbergIRA"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = H.Height();
    const Int numShifts = shifts.Height();
    const Grid& g = H.Grid();
    if( deflate && H.Grid().Rank() == 0 ) 
        std::cerr << "NOTE: Deflation swaps not yet optimized!" << std::endl;

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

    DistMatrix<C,VC,STAR> H_VC_STAR( H );
    DistMatrix<C,VC,STAR> HAdj_VC_STAR( H.Grid() );
    Adjoint( H, HAdj_VC_STAR );
    DistMatrix<C,STAR,VR> activeV_STAR_VR( H.Grid() );
    DistMatrix<C,VR,STAR> activeShiftsConj( H.Grid() );

    // Simultaneously run IRA for different shifts
    std::vector<DistMatrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
    {
        VList[j].SetGrid( H.Grid() );
        Zeros( VList[j], n, numShifts );
    }
    Gaussian( VList[0], n, numShifts );
    const Int numMRShifts = VList[0].LocalWidth();
    std::vector<Matrix<Complex<Real>>> HList(numMRShifts);
    std::vector<Real> realComponents;
    std::vector<Complex<Real>> components;

    DistMatrix<Int,MR,STAR> activeConverged(g);
    Zeros( activeConverged, numShifts, 1 );

    snapCtrl.numSaveCount = 0;
    snapCtrl.imgSaveCount = 0;

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
        HList.resize( activeEsts.LocalHeight() );
        for( Int jLoc=0; jLoc<HList.size(); ++jLoc )
            Zeros( HList[jLoc], basisSize+1, basisSize );

        if( progress )
        {
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
                timer.Start();
        }
        ColumnNorms( activeVList[0], realComponents );
        InvBetaScale( realComponents, activeVList[0] );
        for( Int j=0; j<basisSize; ++j )
        {
            lastActiveEsts = activeEsts;
            activeVList[j+1] = activeVList[j];
            if( progress )
            { 
                mpi::Barrier( H.Grid().Comm() );
                if( H.Grid().Rank() == 0 )
                    subtimer.Start();
            }
            // NOTE: This redistribution sequence might not be necessary
            activeV_STAR_VR = activeVList[j+1];
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), H_VC_STAR, activeShifts, 
              activeV_STAR_VR );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), HAdj_VC_STAR, activeShiftsConj,
              activeV_STAR_VR );
            activeVList[j+1] = activeV_STAR_VR;
            if( progress )
            {
                mpi::Barrier( H.Grid().Comm() );
                if( H.Grid().Rank() == 0 )
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
                ExtractList( HList, components, j, j-1 );
                // TODO: Conjugate components?
                PlaceList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVList[j-1], activeVList[j+1] );
            }

            // Orthogonalize with respect to the last iterate
            InnerProducts( activeVList[j], activeVList[j+1], components );
            PlaceList( HList, components, j, j );
            ColumnSubtractions
            ( components, activeVList[j], activeVList[j+1] );

            // Explicitly (re)orthogonalize against all previous vectors
            for( Int i=0; i<j-1; ++i )
            {
                InnerProducts
                ( activeVList[i], activeVList[j+1], components );
                PlaceList( HList, components, i, j );
                ColumnSubtractions
                ( components, activeVList[i], activeVList[j+1] );
            }
            if( j > 0 )
            {
                InnerProducts
                ( activeVList[j-1], activeVList[j+1], components );
                UpdateList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVList[j-1], activeVList[j+1] );
            }

            // Compute the norm of what is left
            ColumnNorms( activeVList[j+1], realComponents );
            PlaceList( HList, realComponents, j+1, j );

            // TODO: Handle lucky breakdowns
            InvBetaScale( realComponents, activeVList[j+1] );

            ComputeNewEstimates( HList, activeConverged, activeEsts, j+1 );
            // We will have the same estimate two iterations in a row when
            // restarting
            if( j != 0 ) 
                activeConverged =
                    FindConverged
                    ( lastActiveEsts, activeEsts, activeItCounts, tol );

            if( snapCtrl.numFreq > 0 )
                ++snapCtrl.numSaveCount;
            if( snapCtrl.imgFreq > 0 )
                ++snapCtrl.imgSaveCount;
        }
        if( progress )
        {
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
                subtimer.Start();
        }
        Restart( HList, activeConverged, activeVList );
        if( progress )
        {
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
                std::cout << "IRA computations: " << subtimer.Stop()
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
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
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
        Snapshot( estimates, preimage, numIts, deflate, snapCtrl );
    } 

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );

    return itCounts;
}

} // namespace pspec
} // namespace elem

#endif // ifndef ELEM_PSEUDOSPECTRUM_IRA_HPP
