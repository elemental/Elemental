/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PSEUDOSPECTRA_IRL_HPP
#define EL_PSEUDOSPECTRA_IRL_HPP

#include "./Lanczos.hpp"

namespace El {
namespace pspec {

template<typename Real>
void ComputeNewEstimates
( const vector<Matrix<Real>>& HDiagList,
  const vector<Matrix<Real>>& HSubdiagList,
  const Matrix<Int>& activeConverged,
  Matrix<Real>& activeEsts )
{
    DEBUG_CSE
    const Real normCap = NormCap<Real>();
    const Int numShifts = activeEsts.Height();
    if( numShifts == 0 )
        return;
    const Int basisSize = HDiagList[0].size();
    Matrix<Real> HDiag, HSubdiag, w(basisSize);
    for( Int j=0; j<numShifts; ++j )
    {
        HDiag = HDiagList[j];
        HSubdiag = HSubdiagList[j];
        if( !activeConverged(j) )
        {
            if( !HasNan(HDiag) && !HasNan(HSubdiag) )
            {
                lapack::SymmetricTridiagEig
                ( basisSize, HDiag.Buffer(), HSubdiag.Buffer(), w.data(),
                  basisSize-1, basisSize-1 );
                const Real est = Sqrt(w[0]);
                activeEsts(j) = Min(est,normCap);
            }
            else
               activeEsts(j) = normCap;
        }
    }
}

template<typename Real>
void ComputeNewEstimates
( const vector<Matrix<Real>>& HDiagList,
  const vector<Matrix<Real>>& HSubdiagList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
        DistMatrix<Real,MR,STAR>& activeEsts )
{
    DEBUG_CSE
    ComputeNewEstimates
    ( HDiagList, HSubdiagList, activeConverged.LockedMatrix(), 
      activeEsts.Matrix() );
}

template<typename Real>
void Restart
( const vector<Matrix<Real>>& HDiagList,
  const vector<Matrix<Real>>& HSubdiagList,
  const Matrix<Int>& activeConverged,
  vector<Matrix<Complex<Real>>>& VList )
{
    DEBUG_CSE
    const Int n = VList[0].Height();
    const Int numShifts = VList[0].Width();
    if( numShifts == 0 )
        return;
    const Int basisSize = HDiagList[0].size();
    Matrix<Real> HDiag, HSubdiag, w(basisSize);
    Matrix<Real> q(basisSize,1);
    Matrix<Complex<Real>> u(n,1);
    for( Int j=0; j<numShifts; ++j )
    {
        HDiag = HDiagList[j];
        HSubdiag = HSubdiagList[j];

        if( !activeConverged(j) )
        {
            if( !HasNan(HDiag) && !HasNan(HSubdiag) )
            {
                lapack::SymmetricTridiagEig
                ( basisSize, HDiag.Buffer(), HSubdiag.Buffer(), w.data(), 
                  q.Buffer(), basisSize, basisSize-1, basisSize-1 );

                Zeros( u, n, 1 );
                for( Int k=0; k<basisSize; ++k )
                {
                    const Matrix<Complex<Real>>& V = VList[k];
                    auto v = V( ALL, IR(j) ); 
                    Axpy( q(k), v, u );
                }
                Matrix<Complex<Real>>& V = VList[0];
                auto v = V( ALL, IR(j) );
                v = u;
            }
        }
    }
}

template<typename Real>
void Restart
( const vector<Matrix<Real>>& HDiagList,
  const vector<Matrix<Real>>& HSubdiagList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
  vector<DistMatrix<Complex<Real>>>& VList )
{
    DEBUG_CSE
    const Int basisSize = HDiagList[0].size();
    vector<Matrix<Complex<Real>>> VLocList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        View( VLocList[j], VList[j].Matrix() );
    Restart
    ( HDiagList, HSubdiagList, activeConverged.LockedMatrix(), VLocList );
}

template<typename Real>
Matrix<Int>
IRL
( const Matrix<Complex<Real>>& U,
  const Matrix<Complex<Real>>& shifts, 
        Matrix<Real>& invNorms,
        PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_CSE
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
            preimage(j) = j;
    }

    // MultiShiftTrsm currently requires write access
    Matrix<C> UCopy( U );

    // The Hessenberg algorithm currently needs explicit access to the adjoint
    Matrix<C> UAdj;
    if( !psCtrl.schur )
        Adjoint( U, UAdj );

    // Simultaneously run IRL for different shifts
    vector<Matrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        Zeros( VList[j], n, numShifts );
    Gaussian( VList[0], n, numShifts );
    vector<Matrix<Real>> HDiagList(numShifts), HSubdiagList(numShifts);
    Matrix<Real> realComponents;
    Matrix<Complex<Real>> components;
    Matrix<Real> colNorms; 

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
        auto activeShifts = pivShifts( IR(0,numActive), ALL );
        auto activeEsts = estimates( IR(0,numActive), ALL );
        auto activeItCounts = itCounts( IR(0,numActive), ALL );
        for( Int j=0; j<basisSize+1; ++j )
            View( activeVList[j], VList[j], ALL, IR(0,numActive) ); 
        if( deflate )
        {
            View( activePreimage, preimage, IR(0,numActive), ALL );
            Zeros( activeConverged, numActive, 1 );
        }

        // Reset the Rayleigh quotients
        for( Int j=0; j<numActive; ++j )
        {
            HDiagList[j].Resize( 0, 1, basisSize );
            HSubdiagList[j].Resize( 0, 1, basisSize );
        }

        if( progress )
            timer.Start();
        Matrix<Real> colNorms;
        ColumnTwoNorms( activeVList[0], colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, activeVList[0] );
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
                  C(1), UCopy, activeShifts, activeVList[j+1] );
                MultiShiftTrsm
                ( LEFT, UPPER, ADJOINT, 
                  C(1), UCopy, activeShifts, activeVList[j+1] );
                if( progress )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops = 
                        (8.*n*n*numActiveShifts)/(msTime*1.e9);
                    cout << "  MultiShiftTrsm's: " << msTime 
                         << " seconds, " << gflops << " GFlops" << endl;
                }
            }
            else
            {
                if( progress )
                    subtimer.Start();
                MultiShiftHessSolve
                ( UPPER, NORMAL, 
                  C(1), U, activeShifts, activeVList[j+1] );
                Matrix<C> activeShiftsConj;
                Conjugate( activeShifts, activeShiftsConj );
                MultiShiftHessSolve
                ( LOWER, NORMAL, 
                  C(1), UAdj, activeShiftsConj, activeVList[j+1] );
                if( progress )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops = 
                        (32.*n*n*numActiveShifts)/(msTime*1.e9);
                    cout << "  MultiShiftHessSolve's: " << msTime
                         << " seconds, " << gflops << " GFlops" << endl;
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
            ColumnTwoNorms( activeVList[j+1], colNorms );
            PushBackList( HSubdiagList, colNorms );
            // TODO: Handle lucky breakdowns
            DiagonalSolve( RIGHT, NORMAL, colNorms, activeVList[j+1] );

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
            cout << "IRL restart: " << subtimer.Stop()
                 << " seconds" << endl;

        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        numIts += basisSize;
        if( progress )
        {
            const double iterTime = timer.Stop();
            cout << "iteration " << numIts << ": " << iterTime
                 << " seconds, " << numDone << " of " << numShifts
                 << " converged" << endl;
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
DistMatrix<Int,VR,STAR>
IRL
( const ElementalMatrix<Complex<Real>>& UPre, 
  const ElementalMatrix<Complex<Real>>& shiftsPre, 
        ElementalMatrix<Real>& invNormsPre, 
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    DEBUG_CSE
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
    if( !psCtrl.schur )
    {
        U_VC_STAR = U;
        Adjoint( U, UAdj_VC_STAR );
    }

    // Simultaneously run IRL for different shifts
    vector<DistMatrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
    {
        VList[j].SetGrid( g );
        Zeros( VList[j], n, numShifts );
    }
    Gaussian( VList[0], n, numShifts );
    const Int numMRShifts = VList[0].LocalWidth();
    vector<Matrix<Real>> HDiagList(numMRShifts), 
                         HSubdiagList(numMRShifts);
    Matrix<Real> realComponents;
    Matrix<Complex<Real>> components;
    DistMatrix<Real,MR,STAR> colNorms(g);

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
        auto activeShifts = pivShifts( IR(0,numActive), ALL );
        auto activeEsts = estimates( IR(0,numActive), ALL );
        auto activeItCounts = itCounts( IR(0,numActive), ALL );
        for( Int j=0; j<basisSize+1; ++j )
            View( activeVList[j], VList[j], ALL, IR(0,numActive) ); 
        if( deflate )
        {
            View( activePreimage, preimage, IR(0,numActive), ALL );
            Zeros( activeConverged, numActive, 1 );
        }

        // Reset the Rayleigh quotients
        const Int numActiveMR = estimates.LocalHeight();
        for( Int jLoc=0; jLoc<numActiveMR; ++jLoc )
        {
            HDiagList[jLoc].Resize( 0, 1, basisSize );
            HSubdiagList[jLoc].Resize( 0, 1, basisSize );
        }

        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                timer.Start();
        }
        DistMatrix<Real,MR,STAR> colNorms(g);
        ColumnTwoNorms( activeVList[0], colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, activeVList[0] );
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
                DistMatrix<C,STAR,VR> activeV_STAR_VR( activeVList[j+1] );
                MultiShiftHessSolve
                ( UPPER, NORMAL, C(1), U_VC_STAR, activeShifts,
                  activeV_STAR_VR );
                DistMatrix<C,VR,STAR> activeShiftsConj(g);
                Conjugate( activeShifts, activeShiftsConj );
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
                        cout << "  MultiShiftHessSolve's: " << msTime
                             << " seconds, " << gflops << " GFlops" << endl;
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
            ColumnTwoNorms( activeVList[j+1], colNorms );
            PushBackList( HSubdiagList, colNorms.Matrix() );
            // TODO: Handle lucky breakdowns
            DiagonalSolve( RIGHT, NORMAL, colNorms, activeVList[j+1] );

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
                cout << "IRL computations: " << subtimer.Stop()
                     << " seconds" << endl;
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
                cout << "iteration " << numIts << ": " << iterTime
                     << " seconds, " << numDone << " of " << numShifts
                     << " converged" << endl;
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
} // namespace El

#endif // ifndef EL_PSEUDOSPECTRA_IRL_HPP
