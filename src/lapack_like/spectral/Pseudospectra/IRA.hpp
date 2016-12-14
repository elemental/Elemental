/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_PSEUDOSPECTRA_IRA_HPP
#define EL_PSEUDOSPECTRA_IRA_HPP

#include "./Lanczos.hpp"

namespace El {
namespace pspec {

template<typename Real>
void ComputeNewEstimates
( const vector<Matrix<Complex<Real>>>& HList,
  const Matrix<Int>& activeConverged,
        Matrix<Real>& activeEsts,
        Int n )
{
    EL_DEBUG_CSE
    const Real normCap = NormCap<Real>();
    const Int numShifts = activeEsts.Height();
    if( numShifts == 0 )
        return;

    HessenbergSchurCtrl ctrl;
    ctrl.fullTriangle = false;

    Matrix<Complex<Real>> H, HTL;
    Matrix<Complex<Real>> w(n,1);
    for( Int j=0; j<numShifts; ++j )
    {
        H = HList[j];
        HTL = H( IR(0,n), IR(0,n) );
        if( !activeConverged(j) )
        {
            if( !HasNan(HTL) )
            {
                HessenbergSchur( HTL, w, ctrl );
                Real estSquared=0;
                for( Int k=0; k<n; ++k )
                    if( RealPart(w(k)) > estSquared )
                        estSquared = RealPart(w(k));
                activeEsts(j) = Min(Sqrt(estSquared),normCap);
            }
            else
               activeEsts(j) = normCap;
        }
    }
}

template<typename Real>
void ComputeNewEstimates
( const vector<Matrix<Complex<Real>>>& HList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
        DistMatrix<Real,MR,STAR>& activeEsts,
        Int n )
{
    EL_DEBUG_CSE
    ComputeNewEstimates
    ( HList, activeConverged.LockedMatrix(), activeEsts.Matrix(), n );
}

template<typename Real>
void Restart
( const vector<Matrix<Complex<Real>>>& HList,
  const Matrix<Int>& activeConverged,
        vector<Matrix<Complex<Real>>>& VList )
{
    EL_DEBUG_CSE
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
        HTL = H( IR(0,basisSize), IR(0,basisSize) );

        if( !activeConverged(j) )
        {
            if( !HasNan(HTL) )
            {
                HessenbergSchur( HTL, w, Q );

                Real maxReal=0;
                Int maxIdx=0;
                for( Int k=0; k<basisSize; ++k )
                {
                    if( RealPart(w(k)) > maxReal )
                    {
                        maxReal = RealPart(w(k));
                        maxIdx = k;
                    }
                }

                Zeros( u, n, 1 );
                for( Int k=0; k<basisSize; ++k )
                {
                    const Matrix<Complex<Real>>& V = VList[k];
                    auto v = V( ALL, IR(j) );
                    Axpy( Q(k,maxIdx), v, u );
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
( const vector<Matrix<Complex<Real>>>& HList,
  const Matrix<Int>& activeConverged,
        vector<Matrix<Real>>& VRealList,
        vector<Matrix<Real>>& VImagList )
{
    EL_DEBUG_CSE
    const Int n = VRealList[0].Height();
    const Int numShifts = VRealList[0].Width();
    if( numShifts == 0 )
        return;
    const Int basisSize = HList[0].Width();
    Matrix<Complex<Real>> H, HTL, Q(basisSize,basisSize);
    Matrix<Complex<Real>> w(basisSize,1), u(n,1), v(n,1);
    for( Int j=0; j<numShifts; ++j )
    {
        H = HList[j];
        HTL = H( IR(0,basisSize), IR(0,basisSize) );

        if( !activeConverged(j) )
        {
            if( !HasNan(HTL) )
            {
                HessenbergSchur( HTL, w, Q );

                Real maxReal=0;
                Int maxIdx=0;
                for( Int k=0; k<basisSize; ++k )
                {
                    if( RealPart(w(k)) > maxReal )
                    {
                        maxReal = RealPart(w(k));
                        maxIdx = k;
                    }
                }

                Zeros( u, n, 1 );
                for( Int k=0; k<basisSize; ++k )
                {
                    const Matrix<Real>& VReal = VRealList[k];
                    const Matrix<Real>& VImag = VImagList[k];
                    auto vReal = VReal( ALL, IR(j) );
                    auto vImag = VImag( ALL, IR(j) );
                    for( Int i=0; i<n; ++i )
                        v(i) = Complex<Real>(vReal(i),vImag(i));
                    Axpy( Q(k,maxIdx), v, u );
                }
                Matrix<Real>& VReal = VRealList[0];
                Matrix<Real>& VImag = VImagList[0];
                auto vReal = VReal( ALL, IR(j) );
                auto vImag = VImag( ALL, IR(j) );
                for( Int i=0; i<n; ++i )
                {
                    vReal(i) = RealPart(u(i));
                    vImag(i) = ImagPart(u(i));
                }
            }
        }
    }
}

template<typename Real>
void Restart
( const vector<Matrix<Complex<Real>>>& HList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
        vector<DistMatrix<Complex<Real>>>& VList )
{
    EL_DEBUG_CSE
    const Int basisSize = HList[0].Width();
    vector<Matrix<Complex<Real>>> VLocList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        View( VLocList[j], VList[j].Matrix() );
    Restart( HList, activeConverged.LockedMatrix(), VLocList );
}

template<typename Real>
void Restart
( const vector<Matrix<Complex<Real>>>& HList,
  const DistMatrix<Int,MR,STAR>& activeConverged,
        vector<DistMatrix<Real>>& VRealList,
        vector<DistMatrix<Real>>& VImagList )
{
    EL_DEBUG_CSE
    const Int basisSize = HList[0].Width();
    vector<Matrix<Real>> VRealLocList(basisSize+1),
                         VImagLocList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
    {
        View( VRealLocList[j], VRealList[j].Matrix() );
        View( VImagLocList[j], VImagList[j].Matrix() );
    }
    Restart
    ( HList, activeConverged.LockedMatrix(), VRealLocList, VImagLocList );
}

template<typename Real>
Matrix<Int> IRA
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
    const Int basisSize = psCtrl.basisSize;
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

    // MultiShiftTrm requires write access for now...
    Matrix<C> UCopy( U );

    // The Hessenberg variant currently requires explicit access to the adjoint
    Matrix<C> UAdj;
    if( !psCtrl.schur )
        Adjoint( U, UAdj );

    // Simultaneously run IRA for different shifts
    vector<Matrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
        Zeros( VList[j], n, numShifts );
    Gaussian( VList[0], n, numShifts );
    vector<Matrix<Complex<Real>>> HList(numShifts);
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
        HList.resize( numActive );
        for( Int j=0; j<numActive; ++j )
            Zeros( HList[j], basisSize+1, basisSize );

        if( progress )
            timer.Start();
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
                    const double gflops = (8.*n*n*numActiveShifts)/(msTime*1e9);
                    Output
                    ("  MultiShiftTrsm's: ",msTime," seconds, ",
                     gflops," GFlops");
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
                    Output
                    ("  MultiShiftHessSolve's: ",msTime," seconds, ",
                     gflops," GFlops");
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
            ColumnTwoNorms( activeVList[j+1], colNorms );
            PlaceList( HList, colNorms, j+1, j );
            // TODO: Handle lucky breakdowns
            DiagonalSolve( RIGHT, NORMAL, colNorms, activeVList[j+1] );

            ComputeNewEstimates( HList, activeConverged, activeEsts, j+1 );
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
        Restart( HList, activeConverged, activeVList );
        if( progress )
            cout << "IRA restart: " << subtimer.Stop()
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
Matrix<Int>
IRA
( const Matrix<Real>& U,
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
    const Int basisSize = psCtrl.basisSize;
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

    // Simultaneously run IRA for different shifts
    vector<Matrix<Real>> VRealList(basisSize+1),
                         VImagList(basisSize+1),
                         activeVRealList(basisSize+1),
                         activeVImagList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
    {
        Zeros( VRealList[j], n, numShifts );
        Zeros( VImagList[j], n, numShifts );
    }
    // The variance will be off from that of the usual complex case
    Gaussian( VRealList[0], n, numShifts );
    Gaussian( VImagList[0], n, numShifts );
    vector<Matrix<Complex<Real>>> HList(numShifts);
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
        {
            activeVRealList[j] = VRealList[j]( ALL, IR(0,numActive) );
            activeVImagList[j] = VImagList[j]( ALL, IR(0,numActive) );
        }
        if( deflate )
        {
            activePreimage = preimage( IR(0,numActive), ALL );
            Zeros( activeConverged, numActive, 1 );
        }
        HList.resize( numActive );
        for( Int j=0; j<numActive; ++j )
            Zeros( HList[j], basisSize+1, basisSize );

        if( progress )
            timer.Start();
        ColumnTwoNorms( activeVRealList[0], activeVImagList[0], colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, activeVRealList[0] );
        DiagonalSolve( RIGHT, NORMAL, colNorms, activeVImagList[0] );
        for( Int j=0; j<basisSize; ++j )
        {
            lastActiveEsts = activeEsts;
            activeVRealList[j+1] = activeVRealList[j];
            activeVImagList[j+1] = activeVImagList[j];
            if( progress )
                subtimer.Start();
            MultiShiftQuasiTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts,
              activeVRealList[j+1], activeVImagList[j+1] );
            MultiShiftQuasiTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts,
              activeVRealList[j+1], activeVImagList[j+1] );
            if( progress )
            {
                const double msTime = subtimer.Stop();
                const Int numActiveShifts = activeShifts.Height();
                const double gflops = (4.*n*n*numActiveShifts)/(msTime*1.e9);
                cout << "  MultiShiftQuasiTrsm's: " << msTime
                     << " seconds, " << gflops << " GFlops" << endl;
            }

            // Orthogonalize with respect to the old iterate
            if( j > 0 )
            {
                ExtractList( HList, components, j, j-1 );
                // TODO: Conjugate components?
                PlaceList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVRealList[j-1], activeVImagList[j-1],
                              activeVRealList[j+1], activeVImagList[j+1] );
            }

            // Orthogonalize with respect to the last iterate
            InnerProducts
            ( activeVRealList[j  ], activeVImagList[j  ],
              activeVRealList[j+1], activeVImagList[j+1], components );
            PlaceList( HList, components, j, j );
            ColumnSubtractions
            ( components, activeVRealList[j  ], activeVImagList[j  ],
                          activeVRealList[j+1], activeVImagList[j+1] );

            // Explicitly (re)orthogonalize against all previous vectors
            for( Int i=0; i<j-1; ++i )
            {
                InnerProducts
                ( activeVRealList[i  ], activeVImagList[i  ],
                  activeVRealList[j+1], activeVImagList[j+1], components );
                PlaceList( HList, components, i, j );
                ColumnSubtractions
                ( components, activeVRealList[i  ], activeVImagList[i  ],
                              activeVRealList[j+1], activeVImagList[j+1] );
            }
            if( j > 0 )
            {
                InnerProducts
                ( activeVRealList[j-1], activeVImagList[j-1],
                  activeVRealList[j+1], activeVImagList[j+1], components );
                UpdateList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVRealList[j-1], activeVImagList[j-1],
                              activeVRealList[j+1], activeVImagList[j+1] );
            }

            // Compute the norm of what is left
            ColumnTwoNorms
            ( activeVRealList[j+1], activeVImagList[j+1], colNorms );
            PlaceList( HList, colNorms, j+1, j );
            // TODO: Handle lucky breakdowns
            DiagonalSolve( RIGHT, NORMAL, colNorms, activeVRealList[j+1] );
            DiagonalSolve( RIGHT, NORMAL, colNorms, activeVImagList[j+1] );

            ComputeNewEstimates( HList, activeConverged, activeEsts, j+1 );
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
        Restart( HList, activeConverged, activeVRealList, activeVImagList );
        if( progress )
            cout << "IRA restart: " << subtimer.Stop()
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
            ( activeShifts, activePreimage,
              activeVRealList[0], activeVImagList[0], activeEsts,
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
IRA
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
    const Int basisSize = psCtrl.basisSize;
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

    // The Hessenberg case currently requires explicit access to the adjoint
    DistMatrix<C,VC,STAR> U_VC_STAR(g), UAdj_VC_STAR(g);
    if( !psCtrl.schur )
    {
        U_VC_STAR = U;
        Adjoint( U, UAdj_VC_STAR );
    }

    // Simultaneously run IRA for different shifts
    vector<DistMatrix<C>> VList(basisSize+1), activeVList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
    {
        VList[j].SetGrid( g );
        Zeros( VList[j], n, numShifts );
    }
    Gaussian( VList[0], n, numShifts );
    const Int numMRShifts = VList[0].LocalWidth();
    vector<Matrix<Complex<Real>>> HList(numMRShifts);
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
        HList.resize( activeEsts.LocalHeight() );
        for( size_t jLoc=0; jLoc<HList.size(); ++jLoc )
            Zeros( HList[jLoc], basisSize+1, basisSize );

        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                timer.Start();
        }
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
            ColumnTwoNorms( activeVList[j+1], colNorms );
            PlaceList( HList, colNorms.Matrix(), j+1, j );
            // TODO: Handle lucky breakdowns
            DiagonalSolve( RIGHT, NORMAL, colNorms, activeVList[j+1] );

            ComputeNewEstimates( HList, activeConverged, activeEsts, j+1 );
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
        Restart( HList, activeConverged, activeVList );
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                cout << "IRA computations: " << subtimer.Stop()
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

template<typename Real>
DistMatrix<Int,VR,STAR>
IRA
( const AbstractDistMatrix<Real>& UPre,
  const AbstractDistMatrix<Complex<Real>>& shiftsPre,
        AbstractDistMatrix<Real>& invNormsPre,
  PseudospecCtrl<Real> psCtrl=PseudospecCtrl<Real>() )
{
    EL_DEBUG_CSE
    using namespace pspec;
    typedef Complex<Real> C;

    DistMatrixReadProxy<Real,Real,MC,MR> UProx( UPre );
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

    // Simultaneously run IRA for different shifts
    vector<DistMatrix<Real>> VRealList(basisSize+1),
                             VImagList(basisSize+1),
                             activeVRealList(basisSize+1),
                             activeVImagList(basisSize+1);
    for( Int j=0; j<basisSize+1; ++j )
    {
        VRealList[j].SetGrid( g );
        VImagList[j].SetGrid( g );
        Zeros( VRealList[j], n, numShifts );
        Zeros( VImagList[j], n, numShifts );
    }
    // The variance will be off from that of the usual complex case
    Gaussian( VRealList[0], n, numShifts );
    Gaussian( VImagList[0], n, numShifts );
    const Int numMRShifts = VRealList[0].LocalWidth();
    vector<Matrix<Complex<Real>>> HList(numMRShifts);
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
        {
            View( activeVRealList[j], VRealList[j], ALL, IR(0,numActive) );
            View( activeVImagList[j], VImagList[j], ALL, IR(0,numActive) );
        }
        if( deflate )
        {
            View( activePreimage, preimage, IR(0,numActive), ALL );
            Zeros( activeConverged, numActive, 1 );
        }
        HList.resize( activeEsts.LocalHeight() );
        for( size_t jLoc=0; jLoc<HList.size(); ++jLoc )
            Zeros( HList[jLoc], basisSize+1, basisSize );

        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                timer.Start();
        }
        ColumnTwoNorms( activeVRealList[0], activeVImagList[0], colNorms );
        DiagonalSolve( RIGHT, NORMAL, colNorms, activeVRealList[0] );
        DiagonalSolve( RIGHT, NORMAL, colNorms, activeVImagList[0] );
        for( Int j=0; j<basisSize; ++j )
        {
            lastActiveEsts = activeEsts;
            activeVRealList[j+1] = activeVRealList[j];
            activeVImagList[j+1] = activeVImagList[j];
            if( progress )
            {
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                    subtimer.Start();
            }
            MultiShiftQuasiTrsm
            ( LEFT, UPPER, NORMAL, C(1), U, activeShifts,
              activeVRealList[j+1], activeVImagList[j+1] );
            MultiShiftQuasiTrsm
            ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts,
              activeVRealList[j+1], activeVImagList[j+1] );
            if( progress )
            {
                mpi::Barrier( g.Comm() );
                if( g.Rank() == 0 )
                {
                    const double msTime = subtimer.Stop();
                    const Int numActiveShifts = activeShifts.Height();
                    const double gflops =
                        (4.*n*n*numActiveShifts)/(msTime*1.e9);
                    cout << "  MultiShiftQuasiTrsm's: " << msTime
                         << " seconds, " << gflops << " GFlops" << endl;
                }
            }

            // Orthogonalize with respect to the old iterate
            if( j > 0 )
            {
                ExtractList( HList, components, j, j-1 );
                // TODO: Conjugate components?
                PlaceList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVRealList[j-1], activeVImagList[j-1],
                              activeVRealList[j+1], activeVImagList[j+1] );
            }

            // Orthogonalize with respect to the last iterate
            InnerProducts
            ( activeVRealList[j+0], activeVImagList[j+0],
              activeVRealList[j+1], activeVImagList[j+1], components );
            PlaceList( HList, components, j, j );
            ColumnSubtractions
            ( components, activeVRealList[j+0], activeVImagList[j+0],
                          activeVRealList[j+1], activeVImagList[j+1] );

            // Explicitly (re)orthogonalize against all previous vectors
            for( Int i=0; i<j-1; ++i )
            {
                InnerProducts
                ( activeVRealList[i  ], activeVImagList[i  ],
                  activeVRealList[j+1], activeVImagList[j+1], components );
                PlaceList( HList, components, i, j );
                ColumnSubtractions
                ( components, activeVRealList[i  ], activeVImagList[i  ],
                              activeVRealList[j+1], activeVImagList[j+1] );
            }
            if( j > 0 )
            {
                InnerProducts
                ( activeVRealList[j-1], activeVImagList[j-1],
                  activeVRealList[j+1], activeVImagList[j+1], components );
                UpdateList( HList, components, j-1, j );
                ColumnSubtractions
                ( components, activeVRealList[j-1], activeVImagList[j-1],
                              activeVRealList[j+1], activeVImagList[j+1] );
            }

            // Compute the norm of what is left
            ColumnTwoNorms
            ( activeVRealList[j+1], activeVImagList[j+1], colNorms );
            PlaceList( HList, colNorms.Matrix(), j+1, j );
            // TODO: Handle lucky breakdowns
            DiagonalSolve( RIGHT, NORMAL, colNorms, activeVRealList[j+1] );
            DiagonalSolve( RIGHT, NORMAL, colNorms, activeVImagList[j+1] );

            ComputeNewEstimates( HList, activeConverged, activeEsts, j+1 );
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
        Restart( HList, activeConverged, activeVRealList, activeVImagList );
        if( progress )
        {
            mpi::Barrier( g.Comm() );
            if( g.Rank() == 0 )
                cout << "IRA computations: " << subtimer.Stop()
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
            ( activeShifts, activePreimage,
              activeVRealList[0], activeVImagList[0], activeEsts,
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

#endif // ifndef EL_PSEUDOSPECTRA_IRA_HPP
