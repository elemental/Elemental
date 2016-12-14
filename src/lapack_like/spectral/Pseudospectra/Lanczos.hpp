/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PSEUDOSPECTRA_LANCZOS_HPP
#define EL_PSEUDOSPECTRA_LANCZOS_HPP

#include "./Power.hpp"

namespace El {
namespace pspec {

const Int HCapacityInit = 10;

template<typename Real>
void ComputeNewEstimates
( const vector<Matrix<Real>>& HDiagList,
  const vector<Matrix<Real>>& HSubdiagList,
  Matrix<Real>& activeEsts )
{
    EL_DEBUG_CSE
    const Real normCap = NormCap<Real>();
    const Int numShifts = activeEsts.Height();
    if( numShifts == 0 )
        return;
    const BlasInt krylovSize = BlasInt(HDiagList[0].Height());

    // We are only requesting the largest eigenvalue
    HermitianTridiagEigCtrl<Real> ctrl;
    ctrl.subset.indexSubset = true;
    ctrl.subset.lowerIndex = krylovSize-1;
    ctrl.subset.upperIndex = krylovSize-1;

    Matrix<Real> HDiag, HSubdiag, w;
    for( Int j=0; j<numShifts; ++j )
    {
        HDiag = HDiagList[j];
        HSubdiag = HSubdiagList[j];
        if( !HasNan(HDiag) && !HasNan(HSubdiag) )
        {
            HermitianTridiagEig( HDiag, HSubdiag, w, ctrl );
            const Real est = Sqrt(w(0));
            activeEsts(j) = Min(est,normCap);
        }
        else
            activeEsts(j) = normCap;
    }
}

template<typename Real>
void ComputeNewEstimates
( const vector<Matrix<Real>>& HDiagList,
  const vector<Matrix<Real>>& HSubdiagList,
  DistMatrix<Real,MR,STAR>& activeEsts )
{
    EL_DEBUG_CSE
    ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts.Matrix() );
}

template<typename Real>
void Deflate
( vector<Matrix<Real>>& HDiagList,
  vector<Matrix<Real>>& HSubdiagList,
  Matrix<Complex<Real>>& activeShifts,
  Matrix<Int          >& activePreimage,
  Matrix<Complex<Real>>& activeXOld,
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
                std::swap( HDiagList[swapFrom], HDiagList[swapTo] );
                std::swap( HSubdiagList[swapFrom], HSubdiagList[swapTo] );
                RowSwap( activeShifts,   swapFrom, swapTo );
                RowSwap( activePreimage, swapFrom, swapTo );
                RowSwap( activeEsts,     swapFrom, swapTo );
                RowSwap( activeItCounts, swapFrom, swapTo );
                ColSwap( activeXOld,     swapFrom, swapTo );
                ColSwap( activeX,        swapFrom, swapTo );
            }
            --swapTo;
        }
    }
    if( progress )
        cout << "Deflation took " << timer.Stop() << " seconds" << endl;
}

template<typename Real>
void Deflate
( vector<Matrix<Real>>& HDiagList,
  vector<Matrix<Real>>& HSubdiagList,
  DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
  DistMatrix<Complex<Real>        >& activeXOld,
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
    DistMatrix<Complex<Real>,VC,STAR> XOldCopy( activeXOld ), XCopy( activeX );

    auto& convergedLoc = convergedCopy.Matrix();

    const Int n = ( activeX.LocalWidth()>0 ? HDiagList[0].Height() : 0 );
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( convergedLoc(swapFrom) )
        {
            if( swapTo != swapFrom )
            {
                // TODO: Avoid this large latency penalty
                if( activeX.IsLocalCol(swapFrom) &&
                    activeX.IsLocalCol(swapTo) )
                {
                    const Int localFrom = activeX.LocalCol(swapFrom);
                    const Int localTo = activeX.LocalCol(swapTo);
                    EL_DEBUG_ONLY(
                      if( HDiagList[localFrom].Height() != n )
                          LogicError("Invalid HDiagList size");
                      if( HDiagList[localTo].Height() != n )
                          LogicError("Invalid HDiagList size");
                      if( HSubdiagList[localFrom].Height() != n )
                          LogicError("Invalid HSubdiagList size");
                      if( HSubdiagList[localTo].Height() != n )
                          LogicError("Invalid HSubdiagList size");
                    )
                    std::swap( HDiagList[localFrom], HDiagList[localTo] );
                    std::swap( HSubdiagList[localFrom], HSubdiagList[localTo] );
                }
                else if( activeX.IsLocalCol(swapFrom) )
                {
                    const Int localFrom = activeX.LocalCol(swapFrom);
                    EL_DEBUG_ONLY(
                      if( HDiagList[localFrom].Height() != n )
                          LogicError("Invalid HDiagList size");
                      if( HSubdiagList[localFrom].Height() != n )
                          LogicError("Invalid HSubdiagList size");
                    )
                    const Int partner = activeX.ColOwner(swapTo);
                    mpi::TaggedSendRecv
                    ( HDiagList[localFrom].Buffer(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                    mpi::TaggedSendRecv
                    ( HSubdiagList[localFrom].Buffer(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                }
                else if( activeX.IsLocalCol(swapTo) )
                {
                    const Int localTo = activeX.LocalCol(swapTo);
                    EL_DEBUG_ONLY(
                      if( HDiagList[localTo].Height() != n )
                          LogicError("Invalid HDiagList size");
                      if( HSubdiagList[localTo].Height() != n )
                          LogicError("Invalid HSubdiagList size");
                    )
                    const Int partner = activeX.ColOwner(swapFrom);
                    mpi::TaggedSendRecv
                    ( HDiagList[localTo].Buffer(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                    mpi::TaggedSendRecv
                    ( HSubdiagList[localTo].Buffer(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                }

                RowSwap( shiftsCopy,    swapFrom, swapTo );
                RowSwap( preimageCopy,  swapFrom, swapTo );
                RowSwap( estimatesCopy, swapFrom, swapTo );
                RowSwap( itCountsCopy,  swapFrom, swapTo );
                ColSwap( XOldCopy,      swapFrom, swapTo );
                ColSwap( XCopy,         swapFrom, swapTo );
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
            cout << "Deflation took " << timer.Stop() << " seconds" << endl;
    }
}

template<typename Real>
Matrix<Int>
Lanczos
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

    // MultiShiftTrsm requires write access for now
    Matrix<C> UCopy( U );

    // The Hessenberg case currently requires explicit access to the adjoint
    Matrix<C> UAdj;
    if( !psCtrl.schur )
        Adjoint( U, UAdj );

    // Simultaneously run Lanczos for various shifts
    Matrix<C> XOld, X, XNew;
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    vector<Matrix<Real>> HDiagList( numShifts ),
                         HSubdiagList( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        HDiagList[j].Resize( 0, 1, Max(HCapacityInit,1) );
        HSubdiagList[j].Resize( 0, 1, Max(HCapacityInit-1,1) );
    }

    psCtrl.snapCtrl.ResetCounts();

    Timer timer, subtimer;
    Int numIts=0, numDone=0;
    Matrix<Real> estimates(numShifts,1);
    Zeros( estimates, numShifts, 1 );
    auto lastActiveEsts = estimates;
    Matrix<Int> activePreimage;
    Matrix<Real> components;
    Matrix<Real> colNorms;
    while( true )
    {
        const Int numActive = ( deflate ? numShifts-numDone : numShifts );
        auto activeShifts = pivShifts( IR(0,numActive), ALL );
        auto activeEsts = estimates( IR(0,numActive), ALL );
        auto activeItCounts = itCounts( IR(0,numActive), ALL );
        auto activeXOld = XOld( ALL, IR(0,numActive) );
        auto activeX    = X(    ALL, IR(0,numActive) );
        auto activeXNew = XNew( ALL, IR(0,numActive) );
        if( deflate )
            View( activePreimage, preimage, IR(0,numActive), ALL );
        HDiagList.resize( numActive );
        HSubdiagList.resize( numActive );

        if( progress )
            timer.Start();
        activeXNew = activeX;
        if( psCtrl.schur )
        {
            if( progress )
                subtimer.Start();
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), UCopy, activeShifts, activeXNew );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), UCopy, activeShifts, activeXNew );
            if( progress )
            {
                const double msTime = subtimer.Stop();
                const Int numActiveShifts = activeShifts.Height();
                const double gflops = (8.*n*n*numActiveShifts)/(msTime*1.e9);
                cout << "  MultiShiftTrsm's: " << msTime << " seconds, "
                     << gflops << " GFlops" << endl;
            }
        }
        else
        {
            if( progress )
                subtimer.Start();
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), U, activeShifts, activeXNew );
            Matrix<C> activeShiftsConj;
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), UAdj, activeShiftsConj, activeXNew );
            if( progress )
            {
                const double msTime = subtimer.Stop();
                const Int numActiveShifts = activeShifts.Height();
                const double gflops = (32.*n*n*numActiveShifts)/(msTime*1.e9);
                cout << "  MultiShiftHessSolve's: " << msTime
                     << " seconds, " << gflops << " GFlops" << endl;
            }
        }

        // Orthogonalize with respect to the old iterate
        if( numIts > 0 )
        {
            ExtractList( HSubdiagList, components, numIts-1 );
            ColumnSubtractions( components, activeXOld, activeXNew );
        }

        // Orthogonalize with respect to the last iterate
        InnerProducts( activeX, activeXNew, components );
        PushBackList( HDiagList, components );
        ColumnSubtractions( components, activeX, activeXNew );

        // Compute the norm of what is left
        ColumnTwoNorms( activeXNew, colNorms );
        PushBackList( HSubdiagList, colNorms );

        activeXOld = activeX;
        activeX    = activeXNew;
        DiagonalSolve( RIGHT, NORMAL, colNorms, activeX );
        if( progress )
            subtimer.Start();
        ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts );
        if( progress )
            cout << "  Ritz computations: " << subtimer.Stop()
                 << " seconds" << endl;

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
            cout << "iteration " << numIts << ": " << iterTime
                 << " seconds, " << numDone << " of " << numShifts
                 << " converged" << endl;
        }

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate && numActiveDone != 0 )
            Deflate
            ( HDiagList, HSubdiagList, activeShifts, activePreimage, activeXOld,
              activeX, activeEsts, activeConverged, activeItCounts, progress );

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
Lanczos
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

    if( deflate && g.Rank() == 0 )
        cerr << "NOTE: Deflation swaps not yet optimized!" << endl;

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
            preimageLoc(iLoc) = i;
        }
    }

    // The Hessenberg case currently requires explicit access to the adjoint
    DistMatrix<C,VC,STAR> U_VC_STAR(g), UAdj_VC_STAR(g);
    if( !psCtrl.schur )
    {
        U_VC_STAR = U;
        Adjoint( U, UAdj_VC_STAR );
    }

    // Simultaneously run Lanczos for various shifts
    DistMatrix<C> XOld(g), X(g), XNew(g);
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    vector<Matrix<Real>> HDiagList( X.LocalWidth() ),
                         HSubdiagList( X.LocalWidth() );
    for( Int j=0; j<X.LocalWidth(); ++j )
    {
        HDiagList[j].Resize( 0, 1, Max(HCapacityInit,1) );
        HSubdiagList[j].Resize( 0, 1, Max(HCapacityInit-1,1) );
    }

    psCtrl.snapCtrl.ResetCounts();

    Timer timer, subtimer;
    Int numIts=0, numDone=0;
    DistMatrix<Real,MR,STAR> estimates(g);
    estimates.AlignWith( shifts );
    Zeros( estimates, numShifts, 1 );
    auto lastActiveEsts = estimates;
    DistMatrix<Int,VR,STAR> activePreimage(g);
    Matrix<Real> components;
    DistMatrix<Real,MR,STAR> colNorms(g);
    while( true )
    {
        const Int numActive = ( deflate ? numShifts-numDone : numShifts );
        auto activeShifts = pivShifts( IR(0,numActive), ALL );
        auto activeEsts = estimates( IR(0,numActive), ALL );
        auto activeItCounts = itCounts( IR(0,numActive), ALL );
        auto activeXOld = XOld( ALL, IR(0,numActive) );
        auto activeX    = X(    ALL, IR(0,numActive) );
        auto activeXNew = XNew( ALL, IR(0,numActive) );
        if( deflate )
            View( activePreimage, preimage, IR(0,numActive), ALL );
        HDiagList.resize( activeX.LocalWidth() );
        HSubdiagList.resize( activeX.LocalWidth() );

        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                timer.Start();
        }
        activeXNew = activeX;
        if( psCtrl.schur )
        {
            if( progress )
            {
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                    subtimer.Start();
            }
            MultiShiftTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeXNew );
            MultiShiftTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeXNew );
            if( progress )
            {
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops = (8.*n*n*numActiveShifts)/(msTime*1e9);
                    cout << "  MultiShiftTrsm's: " << msTime
                         << " seconds, " << gflops << " GFlops" << endl;
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
            DistMatrix<C,STAR,VR> activeXNew_STAR_VR( activeXNew );
            MultiShiftHessSolve
            ( UPPER, NORMAL, C(1), U_VC_STAR, activeShifts,
              activeXNew_STAR_VR );
            DistMatrix<C,VR,STAR> activeShiftsConj(g);
            Conjugate( activeShifts, activeShiftsConj );
            MultiShiftHessSolve
            ( LOWER, NORMAL, C(1), UAdj_VC_STAR, activeShiftsConj,
              activeXNew_STAR_VR );
            activeXNew = activeXNew_STAR_VR;
            if( progress )
            {
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops =
                        (32.*n*n*numActiveShifts)/(msTime*1.e9);
                    cout << "  MultiShiftHessSolve's: " << msTime
                         << " seconds, " << gflops << " GFlops" << endl;
                }
            }
        }

        // Orthogonalize with respect to the old iterate
        if( numIts > 0 )
        {
            ExtractList( HSubdiagList, components, numIts-1 );
            ColumnSubtractions( components, activeXOld, activeXNew );
        }

        // Orthogonalize with respect to the last iterate
        InnerProducts( activeX, activeXNew, components );
        PushBackList( HDiagList, components );
        ColumnSubtractions( components, activeX, activeXNew );

        // Compute the norm of what is left
        ColumnTwoNorms( activeXNew, colNorms );
        PushBackList( HSubdiagList, colNorms.Matrix() );

        activeXOld = activeX;
        activeX    = activeXNew;
        DiagonalSolve( RIGHT, NORMAL, colNorms, activeX );
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                subtimer.Start();
        }
        ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts );
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                cout << "  Ritz computations: " << subtimer.Stop()
                     << " seconds" << endl;
        }

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
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
            {
                const double iterTime = timer.Stop();
                cout << "iteration " << numIts << ": " << iterTime
                     << " seconds, " << numDone << " of " << numShifts
                     << " converged" << endl;
            }
        }

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate && numActiveDone != 0 )
            Deflate
            ( HDiagList, HSubdiagList, activeShifts, activePreimage, activeXOld,
              activeX, activeEsts, activeConverged, activeItCounts, progress );

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

#endif // ifndef EL_PSEUDOSPECTRA_LANCZOS_HPP
