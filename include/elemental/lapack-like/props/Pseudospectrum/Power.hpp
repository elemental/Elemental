/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_POWER_HPP
#define ELEM_PSEUDOSPECTRUM_POWER_HPP

#include "./Util.hpp"

namespace elem {
namespace pspec {

template<typename Real>
inline void
Deflate
( Matrix<Complex<Real>>& activeShifts, 
  Matrix<Int          >& activePreimage,
  Matrix<Complex<Real>>& activeX,
  Matrix<Real         >& activeEsts, 
  Matrix<Int          >& activeConverged,
  Matrix<Int          >& activeItCounts,
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
                RowSwap( activeShifts,   swapFrom, swapTo );
                RowSwap( activePreimage, swapFrom, swapTo );
                RowSwap( activeEsts,     swapFrom, swapTo );
                RowSwap( activeItCounts, swapFrom, swapTo );
                ColSwap( activeX,        swapFrom, swapTo );
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
( Matrix<Complex<Real>>& activeShifts,
  Matrix<Int          >& activePreimage,
  Matrix<Real         >& activeXReal,
  Matrix<Real         >& activeXImag,
  Matrix<Real         >& activeEsts,
  Matrix<Int          >& activeConverged,
  Matrix<Int          >& activeItCounts,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Deflate"))
    Timer timer;
    if( progress )
        timer.Start();
    const Int numActive = activeXReal.Width();
    Int swapTo = numActive-1;
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( activeConverged.Get(swapFrom,0) )
        {
            if( swapTo != swapFrom )
            {
                RowSwap( activeShifts,   swapFrom, swapTo );
                RowSwap( activePreimage, swapFrom, swapTo );
                RowSwap( activeEsts,     swapFrom, swapTo );
                RowSwap( activeItCounts, swapFrom, swapTo );
                ColSwap( activeXReal,    swapFrom, swapTo );
                ColSwap( activeXImag,    swapFrom, swapTo );
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
( DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
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
    DistMatrix<Complex<Real>,VC,STAR> XCopy( activeX );

    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( convergedCopy.Get(swapFrom,0) )
        {
            if( swapTo != swapFrom )
            {
                RowSwap( shiftsCopy,    swapFrom, swapTo );
                RowSwap( preimageCopy,  swapFrom, swapTo );
                RowSwap( estimatesCopy, swapFrom, swapTo );
                RowSwap( itCountsCopy,  swapFrom, swapTo );
                ColSwap( XCopy,         swapFrom, swapTo );
            }
            --swapTo;
        }
    }

    activeShifts   = shiftsCopy;
    activePreimage = preimageCopy;
    activeEsts     = estimatesCopy;
    activeItCounts = itCountsCopy;
    activeX        = XCopy;

    if( progress && activeShifts.Grid().Rank() == 0 )
        std::cout << "Deflation took " << timer.Stop() << " seconds"
                  << std::endl;
}

template<typename Real>
inline void
Deflate
( DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
  DistMatrix<Real                 >& activeXReal,
  DistMatrix<Real                 >& activeXImag,
  DistMatrix<Real,         MR,STAR>& activeEsts,
  DistMatrix<Int,          MR,STAR>& activeConverged,
  DistMatrix<Int,          VR,STAR>& activeItCounts,
  bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Deflate"))
    Timer timer;
    if( progress && activeShifts.Grid().Rank() == 0 )
        timer.Start();
    const Int numActive = activeXReal.Width();
    Int swapTo = numActive-1;

    DistMatrix<Complex<Real>,STAR,STAR> shiftsCopy( activeShifts );
    DistMatrix<Int,STAR,STAR> preimageCopy( activePreimage );
    DistMatrix<Real,STAR,STAR> estimatesCopy( activeEsts );
    DistMatrix<Int, STAR,STAR> itCountsCopy( activeItCounts );
    DistMatrix<Int, STAR,STAR> convergedCopy( activeConverged );
    DistMatrix<Real,VC,STAR> XRealCopy( activeXReal ),
                             XImagCopy( activeXImag );

    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( convergedCopy.Get(swapFrom,0) )
        {
            if( swapTo != swapFrom )
            {
                RowSwap( shiftsCopy,    swapFrom, swapTo );
                RowSwap( preimageCopy,  swapFrom, swapTo );
                RowSwap( estimatesCopy, swapFrom, swapTo );
                RowSwap( itCountsCopy,  swapFrom, swapTo );
                ColSwap( XRealCopy,     swapFrom, swapTo );
                ColSwap( XImagCopy,     swapFrom, swapTo );
            }
            --swapTo;
        }
    }

    activeShifts   = shiftsCopy;
    activePreimage = preimageCopy;
    activeEsts     = estimatesCopy;
    activeItCounts = itCountsCopy;
    activeXReal    = XRealCopy;
    activeXImag    = XImagCopy;

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
Power
( const Matrix<Complex<Real>>& U, const Matrix<Complex<Real>>& shifts, 
  Matrix<Real>& invNorms, PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Power"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();

    const Int maxIts = psCtrl.maxIts;
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

    psCtrl.snapCtrl.ResetCounts();

    // The Hessenberg case currently requires explicit access to the adjoint
    Matrix<C> UAdj, activeShiftsConj;
    if( !psCtrl.schur )
        Adjoint( U, UAdj );

    // Simultaneously run inverse iteration for various shifts
    Timer timer;
    Matrix<C> X;
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Int numIts=0, numDone=0;
    Matrix<Real> estimates(numShifts,1);
    Zeros( estimates, numShifts, 1 );
    auto lastActiveEsts = estimates;
    Matrix<Int> activePreimage;
    while( true )
    {
        const Int numActive = ( deflate ? numShifts-numDone : numShifts );
        auto activeShifts = View( pivShifts, 0, 0, numActive, 1 );
        auto activeEsts = View( estimates, 0, 0, numActive, 1 );
        auto activeItCounts = View( itCounts, 0, 0, numActive, 1 );
        auto activeX = View( X, 0, 0, n, numActive );
        if( deflate )
            activePreimage = View( preimage, 0, 0, numActive, 1 );

        if( progress )
            timer.Start(); 

        if( psCtrl.schur )
        {
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeX );
            FixColumns( activeX );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeX );
        }
        else
        {
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), U, activeShifts, activeX );
            FixColumns( activeX );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), UAdj, activeShiftsConj, activeX );
        }

        ColumnNorms( activeX, activeEsts );
        CapEstimates( activeEsts );

        auto activeConverged = 
            FindConverged
            ( lastActiveEsts, activeEsts, activeItCounts, psCtrl.tol );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        if( progress )
        {
            const double iterTime = timer.Stop();
            std::cout << "iteration " << numIts << ": " << iterTime 
                      << " seconds, " << numDone << " of " << numShifts 
                      << " converged" << std::endl;
        }

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate && numActiveDone != 0 )
            Deflate
            ( activeShifts, activePreimage, activeX, activeEsts,
              activeConverged, activeItCounts, progress );

        lastActiveEsts = activeEsts;

        // Save snapshots of the estimates at the requested rate
        psCtrl.snapCtrl.Iterate();
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
Power
( const DistMatrix<Complex<Real>        >& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Power"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();
    const Grid& g = U.Grid();

    const Int maxIts = psCtrl.maxIts;
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

    psCtrl.snapCtrl.ResetCounts();

    // The Hessenberg case currently requires explicit access to the adjoint
    DistMatrix<C,VC,STAR> U_VC_STAR(g), UAdj_VC_STAR(g);
    DistMatrix<C,VR,STAR> activeShiftsConj(g);
    DistMatrix<C,STAR,VR> activeX_STAR_VR(g);
    if( !psCtrl.schur )
    {
        U_VC_STAR = U;
        Adjoint( U, UAdj_VC_STAR );
    }

    // Simultaneously run inverse iteration for various shifts
    Timer timer;
    DistMatrix<C> X(g);
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Int numIts=0, numDone=0;
    DistMatrix<Real,MR,STAR> estimates(g);
    estimates.AlignWith( shifts );
    Zeros( estimates, numShifts, 1 );
    auto lastActiveEsts = estimates;
    DistMatrix<Int,VR,STAR> activePreimage(g);
    while( true )
    {
        const Int numActive = ( deflate ? numShifts-numDone : numShifts );
        auto activeShifts = View( pivShifts, 0, 0, numActive, 1 );
        auto activeEsts = View( estimates, 0, 0, numActive, 1 );
        auto activeItCounts = View( itCounts, 0, 0, numActive, 1 );
        auto activeX = View( X, 0, 0, n, numActive );
        if( deflate )
            activePreimage = View( preimage, 0, 0, numActive, 1 );

        if( progress && g.Rank() == 0 )
            timer.Start();
        if( psCtrl.schur )
        {
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeX );
            FixColumns( activeX );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeX );
        }
        else
        {
            Conjugate( activeShifts, activeShiftsConj );
            activeX_STAR_VR = activeX;
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), U_VC_STAR, activeShifts,
              activeX_STAR_VR );
            FixColumns( activeX_STAR_VR );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), UAdj_VC_STAR, activeShiftsConj,
              activeX_STAR_VR );
            activeX = activeX_STAR_VR;
        }
        ColumnNorms( activeX, activeEsts );
        CapEstimates( activeEsts );

        auto activeConverged =
            FindConverged
            ( lastActiveEsts, activeEsts, activeItCounts, psCtrl.tol );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        if( progress && g.Rank() == 0 )
        {
            const double iterTime = timer.Stop();
            std::cout << "iteration " << numIts << ": " << iterTime 
                      << " seconds, " << numDone << " of " << numShifts 
                      << " converged" << std::endl;
        }

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate && numActiveDone != 0 )
            Deflate
            ( activeShifts, activePreimage, activeX, activeEsts,
              activeConverged, activeItCounts, progress );

        lastActiveEsts = activeEsts;

        // Save snapshots of the estimates at the requested rate
        psCtrl.snapCtrl.Iterate();
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

#endif // ifndef ELEM_PSEUDOSPECTRUM_POWER_HPP
