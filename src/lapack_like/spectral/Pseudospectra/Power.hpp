/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PSEUDOSPECTRA_POWER_HPP
#define EL_PSEUDOSPECTRA_POWER_HPP

#include "./Util.hpp"

namespace El {
namespace pspec {

template<typename Real>
void Deflate
( Matrix<Complex<Real>>& activeShifts,
  Matrix<Int          >& activePreimage,
  Matrix<Complex<Real>>& activeX,
  Matrix<Real         >& activeEsts,
  Matrix<Int          >& activeConverged,
  Matrix<Int          >& activeItCounts,
  bool progress=false )
{
    EL_DEBUG_CSE
    Timer timer;
    if( progress )
        timer.Start();
    const Int numActive = activeX.Width();
    Int swapTo = numActive-1;
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( activeConverged(swapFrom) )
        {
            if( swapTo != swapFrom )
            {
                RowSwap( activeShifts,   swapFrom, swapTo );
                RowSwap( activePreimage, swapFrom, swapTo );
                RowSwap( activeEsts,     swapFrom, swapTo );
                RowSwap( activeItCounts, swapFrom, swapTo );
                // NOTE: We only need to move information one way
                ColSwap( activeX,        swapFrom, swapTo );
            }
            --swapTo;
        }
    }
    if( progress )
        Output("Deflation took ",timer.Stop()," seconds");
}

template<typename Real>
void Deflate
( Matrix<Complex<Real>>& activeShifts,
  Matrix<Int          >& activePreimage,
  Matrix<Real         >& activeXReal,
  Matrix<Real         >& activeXImag,
  Matrix<Real         >& activeEsts,
  Matrix<Int          >& activeConverged,
  Matrix<Int          >& activeItCounts,
  bool progress=false )
{
    EL_DEBUG_CSE
    Timer timer;
    if( progress )
        timer.Start();
    const Int numActive = activeXReal.Width();
    Int swapTo = numActive-1;
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( activeConverged(swapFrom) )
        {
            if( swapTo != swapFrom )
            {
                RowSwap( activeShifts,   swapFrom, swapTo );
                RowSwap( activePreimage, swapFrom, swapTo );
                RowSwap( activeEsts,     swapFrom, swapTo );
                RowSwap( activeItCounts, swapFrom, swapTo );
                // NOTE: We only need to move information one way
                ColSwap( activeXReal,    swapFrom, swapTo );
                ColSwap( activeXImag,    swapFrom, swapTo );
            }
            --swapTo;
        }
    }
    if( progress )
        Output("Deflation took ",timer.Stop()," seconds");
}

template<typename Real>
void Deflate
( DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
  DistMatrix<Complex<Real>        >& activeX,
  DistMatrix<Real,         MR,STAR>& activeEsts,
  DistMatrix<Int,          MR,STAR>& activeConverged,
  DistMatrix<Int,          VR,STAR>& activeItCounts,
  bool progress=false )
{
    EL_DEBUG_CSE
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
    auto& convergedLoc = convergedCopy.Matrix();

    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( convergedLoc(swapFrom) )
        {
            if( swapTo != swapFrom )
            {
                RowSwap( shiftsCopy,    swapFrom, swapTo );
                RowSwap( preimageCopy,  swapFrom, swapTo );
                RowSwap( estimatesCopy, swapFrom, swapTo );
                RowSwap( itCountsCopy,  swapFrom, swapTo );
                // NOTE: We only need to move information one way
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
        Output("Deflation took ",timer.Stop()," seconds");
}

template<typename Real>
void Deflate
( DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
  DistMatrix<Real                 >& activeXReal,
  DistMatrix<Real                 >& activeXImag,
  DistMatrix<Real,         MR,STAR>& activeEsts,
  DistMatrix<Int,          MR,STAR>& activeConverged,
  DistMatrix<Int,          VR,STAR>& activeItCounts,
  bool progress=false )
{
    EL_DEBUG_CSE
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
    auto& convergedLoc = convergedCopy.Matrix();

    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( convergedLoc(swapFrom) )
        {
            if( swapTo != swapFrom )
            {
                RowSwap( shiftsCopy,    swapFrom, swapTo );
                RowSwap( preimageCopy,  swapFrom, swapTo );
                RowSwap( estimatesCopy, swapFrom, swapTo );
                RowSwap( itCountsCopy,  swapFrom, swapTo );
                // NOTE: We only need to move information one way
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
            Output("Deflation took ",timer.Stop()," seconds");
    }
}

template<typename Real>
Matrix<Int>
Power
( const Matrix<Complex<Real>>& U,
  const Matrix<Complex<Real>>& shifts,
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    EL_DEBUG_CSE
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
            preimage(j) = j;
    }

    psCtrl.snapCtrl.ResetCounts();

    // MultiShiftTrsm currently requires write access
    Matrix<C> UCopy( U );

    // The Hessenberg case currently requires explicit access to the adjoint
    Matrix<C> UAdj;
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
        auto activeShifts = pivShifts( IR(0,numActive), ALL );
        auto activeEsts = estimates( IR(0,numActive), ALL );
        auto activeItCounts = itCounts( IR(0,numActive), ALL );
        auto activeX = X( ALL, IR(0,numActive) );
        if( deflate )
            View( activePreimage, preimage, IR(0,numActive), ALL );

        if( progress )
            timer.Start();

        if( psCtrl.schur )
        {
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), UCopy, activeShifts, activeX );
            FixColumns( activeX );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), UCopy, activeShifts, activeX );
        }
        else
        {
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), U, activeShifts, activeX );
            FixColumns( activeX );
            Matrix<C> activeShiftsConj;
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), UAdj, activeShiftsConj, activeX );
        }

        ColumnTwoNorms( activeX, activeEsts );
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
            Output
            ("iteration ",numIts,": ",iterTime," seconds, ",
             numDone," of ",numShifts," converged");
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
DistMatrix<Int,VR,STAR>
Power
( const AbstractDistMatrix<Complex<Real>>& UPre,
  const AbstractDistMatrix<Complex<Real>>& shiftsPre,
        AbstractDistMatrix<Real>& invNormsPre,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    EL_DEBUG_CSE
    using namespace pspec;
    typedef Complex<Real> C;

    DistMatrixReadProxy<C,C,MC,MR> UProx( UPre );
    DistMatrixReadProxy<C,C,VR,STAR> shiftsProx( shiftsPre );
    DistMatrixWriteProxy<Real,Real,VR,STAR> invNormsProx( invNormsPre );
    auto& U = UProx.GetLocked();
    auto& shifts = shiftsProx.GetLocked();
    auto& invNorms = invNormsProx.Get();

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
        auto& preimageLoc = preimage.Matrix();
        const Int numLocShifts = preimage.LocalHeight();
        for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
        {
            const Int i = preimage.GlobalRow(iLoc);
            preimageLoc(iLoc) =  i;
        }
    }

    psCtrl.snapCtrl.ResetCounts();

    // The Hessenberg case currently requires explicit access to the adjoint
    DistMatrix<C,VC,STAR> U_VC_STAR(g), UAdj_VC_STAR(g);
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
        auto activeShifts = pivShifts( IR(0,numActive), ALL );
        auto activeEsts = estimates( IR(0,numActive), ALL );
        auto activeItCounts = itCounts( IR(0,numActive), ALL );
        auto activeX = X( ALL, IR(0,numActive) );
        if( deflate )
            View( activePreimage, preimage, IR(0,numActive), ALL );

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
            DistMatrix<C,STAR,VR> activeX_STAR_VR( activeX );
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), U_VC_STAR, activeShifts,
              activeX_STAR_VR );
            FixColumns( activeX_STAR_VR );
            DistMatrix<C,VR,STAR> activeShiftsConj(g);
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), UAdj_VC_STAR, activeShiftsConj,
              activeX_STAR_VR );
            activeX = activeX_STAR_VR;
        }
        ColumnTwoNorms( activeX, activeEsts );
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
            Output
            ("iteration ",numIts,": ",iterTime," seconds, ",
             numDone," of ",numShifts," converged");
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
} // namespace El

#endif // ifndef EL_PSEUDOSPECTRA_POWER_HPP
