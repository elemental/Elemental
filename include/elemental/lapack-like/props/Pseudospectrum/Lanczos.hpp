/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PSEUDOSPECTRUM_LANCZOS_HPP
#define ELEM_PSEUDOSPECTRUM_LANCZOS_HPP

#include "./Power.hpp"

namespace elem {
namespace pspec {

const Int HCapacityInit = 10;

template<typename Real>
inline bool HasNan
( const std::vector<Real>& HDiag, 
  const std::vector<Real>& HSubdiag )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::HasNan"))
    bool hasNan = false;
    const Int n = HDiag.size();
    for( Int j=0; j<n; ++j )
        if( std::isnan( HDiag[j] ) )
            hasNan = true;
    for( Int j=0; j<n-1; ++j )
        if( std::isnan( HSubdiag[j] ) )
            hasNan = true;
    return hasNan;
}

template<typename F,typename FComp>
inline void
ColumnSubtractions
( const std::vector<FComp>& components,
  const Matrix<F>& X,
        Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnSubtractions"))
    const Int numShifts = Y.Width();
    if( numShifts == 0 )
        return;
    const Int m = Y.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const F gamma = components[j];
        blas::Axpy( m, -gamma, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
    }
}

template<typename F,typename FComp>
inline void
ColumnSubtractions
( const std::vector<FComp>& components,
  const DistMatrix<F>& X,
        DistMatrix<F>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnSubtractions");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    ColumnSubtractions( components, X.LockedMatrix(), Y.Matrix() );
}

template<typename F,typename FComp>
inline void
ColumnSubtractions
( const std::vector<std::vector<FComp> >& diags,
  const Matrix<F>& X, 
        Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnSubtractions"))
    const Int numShifts = Y.Width();
    if( numShifts == 0 )
        return;
    const Int m = Y.Height();
    const Int n = diags[0].size();
    for( Int j=0; j<numShifts; ++j )
    {
        const F gamma = diags[j][n-1];
        blas::Axpy( m, -gamma, X.LockedBuffer(0,j), 1, Y.Buffer(0,j), 1 );
    }
}

template<typename F,typename FComp>
inline void
ColumnSubtractions
( const std::vector<std::vector<FComp> >& diags,
  const DistMatrix<F>& X, 
        DistMatrix<F>& Y )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnSubtractions");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    ColumnSubtractions( diags, X.LockedMatrix(), Y.Matrix() );
}

template<typename F>
inline void
InnerProducts
( const Matrix<F>& X,
  const Matrix<F>& Y,
        std::vector<BASE(F)>& innerProds )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InnerProducts"))
    typedef Base<F> Real;
    const Int numShifts = X.Width();
    const Int m = X.Height();
    innerProds.resize( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        const Real alpha =
            RealPart(blas::Dot( m, X.LockedBuffer(0,j), 1,
                                   Y.LockedBuffer(0,j), 1 ));
        innerProds[j] = alpha;
    }
}

template<typename F>
inline void
InnerProducts
( const Matrix<F>& X,
  const Matrix<F>& Y,
        std::vector<F>& innerProds )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InnerProducts"))
    const Int numShifts = X.Width();
    const Int m = X.Height();
    innerProds.resize( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        const F alpha =
            blas::Dot( m, X.LockedBuffer(0,j), 1,
                          Y.LockedBuffer(0,j), 1 );
        innerProds[j] = alpha;
    }
}

template<typename F>
inline void
InnerProducts
( const DistMatrix<F>& X,
  const DistMatrix<F>& Y,
        std::vector<BASE(F)>& innerProds )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnNorms");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.data(), numLocShifts, mpi::SUM, X.ColComm() );
}

template<typename F>
inline void
InnerProducts
( const DistMatrix<F>& X,
  const DistMatrix<F>& Y,
        std::vector<F>& innerProds )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnNorms");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() )
            LogicError("X and Y should have been aligned");
    )
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.data(), numLocShifts, mpi::SUM, X.ColComm() );
}

template<typename F>
inline void
InnerProducts
( const Matrix<F>& X,
  const Matrix<F>& Y,
        std::vector<std::vector<BASE(F)> >& diags )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InnerProducts"))
    typedef Base<F> Real;
    const Int numShifts = X.Width();
    const Int m = X.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Real alpha = 
            RealPart(blas::Dot( m, X.LockedBuffer(0,j), 1, 
                                   Y.LockedBuffer(0,j), 1 ));
        diags[j].push_back( alpha );
    }
}

template<typename F>
inline void
InnerProducts
( const Matrix<F>& X,
  const Matrix<F>& Y,
        std::vector<std::vector<F> >& diags )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InnerProducts"))
    const Int numShifts = X.Width();
    const Int m = X.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const F alpha = 
            blas::Dot( m, X.LockedBuffer(0,j), 1, 
                          Y.LockedBuffer(0,j), 1 );
        diags[j].push_back( alpha );
    }
}

template<typename F>
inline void
InnerProducts
( const DistMatrix<F>& X,
  const DistMatrix<F>& Y,
        std::vector<std::vector<BASE(F)> >& diags )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnNorms");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() ) 
            LogicError("X and Y should have been aligned");
    )
    std::vector<Base<F>> innerProds;
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.data(), numLocShifts, mpi::SUM, X.ColComm() );
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        diags[jLoc].push_back( innerProds[jLoc] );
}

template<typename F>
inline void
InnerProducts
( const DistMatrix<F>& X,
  const DistMatrix<F>& Y,
        std::vector<std::vector<F> >& diags )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnNorms");
        if( X.ColAlign() != Y.ColAlign() || X.RowAlign() != Y.RowAlign() ) 
            LogicError("X and Y should have been aligned");
    )
    std::vector<F> innerProds;
    InnerProducts( X.LockedMatrix(), Y.LockedMatrix(), innerProds );
    const Int numLocShifts = X.LocalWidth();
    mpi::AllReduce( innerProds.data(), numLocShifts, mpi::SUM, X.ColComm() );
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        diags[jLoc].push_back( innerProds[jLoc] );
}

template<typename F>
inline void
ColumnNorms
( const Matrix<F>& X,
        std::vector<std::vector<BASE(F)> >& diags )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    typedef Base<F> Real;
    Matrix<Real> norms;
    ColumnNorms( X, norms );

    const Int numShifts = X.Width();
    for( Int j=0; j<numShifts; ++j )
        diags[j].push_back( norms.Get(j,0) );
}

template<typename F>
inline void
ColumnNorms
( const DistMatrix<F>& X,
        std::vector<std::vector<BASE(F)> >& diags )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    typedef Base<F> Real;
    DistMatrix<Real,MR,STAR> norms( X.Grid() );
    ColumnNorms( X, norms );

    const Int numLocShifts = X.LocalWidth();
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        diags[jLoc].push_back( norms.GetLocal(jLoc,0) );
}

template<typename F>
inline void
InvBetaScale( const Matrix<BASE(F)>& scales, Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InvBetaScale"))
    typedef Base<F> Real;
    const Int numShifts = Y.Width();
    if( numShifts == 0 )
        return;
    const Int m = Y.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Real beta = scales.Get(j,0);
        blas::Scal( m, F(1)/beta, Y.Buffer(0,j), 1 );
    }
}

template<typename F>
inline void
InvBetaScale( const DistMatrix<BASE(F),MR,STAR>& scales, DistMatrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InvBetaScale"))
    InvBetaScale( scales.LockedMatrix(), Y.Matrix() );
}

template<typename F>
inline void
InvBetaScale
( const std::vector<std::vector<BASE(F)> >& HSubdiagList,
        Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InvBetaScale"))
    typedef Base<F> Real;
    const Int numShifts = Y.Width();
    if( numShifts == 0 )
        return;
    const Int m = Y.Height();
    const Int n = HSubdiagList[0].size();
    for( Int j=0; j<numShifts; ++j )
    {
        const Real beta = HSubdiagList[j][n-1];
        blas::Scal( m, F(1)/beta, Y.Buffer(0,j), 1 );
    }
}

template<typename F>
inline void
InvBetaScale
( const std::vector<std::vector<BASE(F)> >& HSubdiagList,
        DistMatrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InvBetaScale"))
    InvBetaScale( HSubdiagList, Y.Matrix() );
}

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real> >& HDiagList, 
  const std::vector<std::vector<Real> >& HSubdiagList,
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

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real> >& HDiagList, 
  const std::vector<std::vector<Real> >& HSubdiagList,
  DistMatrix<Real,MR,STAR>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts.Matrix() );
}

template<typename Real>
inline void
Deflate
( std::vector<std::vector<Real> >& HDiagList,
  std::vector<std::vector<Real> >& HSubdiagList,
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
                std::swap( HDiagList[swapFrom], HDiagList[swapTo] );
                std::swap( HSubdiagList[swapFrom], HSubdiagList[swapTo] );
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
( std::vector<std::vector<Real> >& HDiagList,
  std::vector<std::vector<Real> >& HSubdiagList,
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

    const Int n = ( activeX.LocalWidth()>0 ? HDiagList[0].size() : 0 );
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
                    DEBUG_ONLY(
                        if( HDiagList[localFrom].size() != n )
                            LogicError("Invalid HDiagList size");
                        if( HDiagList[localTo].size() != n )
                            LogicError("Invalid HDiagList size");
                        if( HSubdiagList[localFrom].size() != n )
                            LogicError("Invalid HSubdiagList size");
                        if( HSubdiagList[localTo].size() != n )
                            LogicError("Invalid HSubdiagList size");
                    )
                    std::swap( HDiagList[localFrom], HDiagList[localTo] );
                    std::swap( HSubdiagList[localFrom], HSubdiagList[localTo] );
                }
                else if( activeX.IsLocalCol(swapFrom) )
                {
                    const Int localFrom = activeX.LocalCol(swapFrom);
                    DEBUG_ONLY(
                        if( HDiagList[localFrom].size() != n )
                            LogicError("Invalid HDiagList size");
                        if( HSubdiagList[localFrom].size() != n )
                            LogicError("Invalid HSubdiagList size");
                    )
                    const Int partner = activeX.ColOwner(swapTo);
                    mpi::TaggedSendRecv
                    ( HDiagList[localFrom].data(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                    mpi::TaggedSendRecv
                    ( HSubdiagList[localFrom].data(), n, 
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                }
                else if( activeX.IsLocalCol(swapTo) )
                {
                    const Int localTo = activeX.LocalCol(swapTo);
                    DEBUG_ONLY(
                        if( HDiagList[localTo].size() != n )
                            LogicError("Invalid HDiagList size");
                        if( HSubdiagList[localTo].size() != n )
                            LogicError("Invalid HSubdiagList size");
                    )
                    const Int partner = activeX.ColOwner(swapFrom);
                    mpi::TaggedSendRecv
                    ( HDiagList[localTo].data(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                    mpi::TaggedSendRecv
                    ( HSubdiagList[localTo].data(), n, 
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
TriangularLanczos
( const Matrix<Complex<Real> >& U, const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, 
  bool deflate=true, Int maxIts=1000, Real tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::TriangularLanczos"))
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

    // Simultaneously run Lanczos for various shifts
    Matrix<C> XOld, X, XNew;
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    std::vector<std::vector<Real>> HDiagList( numShifts ),
                                   HSubdiagList( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        HDiagList[j].reserve( HCapacityInit );
        HSubdiagList[j].reserve( HCapacityInit-1 );
    }

    Timer timer, subtimer;
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
        auto activeXOld = View( XOld, 0, 0, n, numActive );
        auto activeX    = View( X,    0, 0, n, numActive );
        auto activeXNew = View( XNew, 0, 0, n, numActive );
        if( deflate )
            activePreimage = View( preimage, 0, 0, numActive, 1 );

        if( progress )
            timer.Start();
        activeXNew = activeX;
        if( progress )
            subtimer.Start();
        MultiShiftTrsm
        ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeXNew );
        MultiShiftTrsm
        ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeXNew );
        if( progress )
        {
            const double msTime = subtimer.Stop();
            const double gflops = (8.*n*n*numShifts)/(msTime*1.e9);
            std::cout << "  MultiShiftTrsm's: " << msTime << " seconds, "
                      << gflops << " GFlops" << std::endl;
        }
        ColumnSubtractions( HSubdiagList, activeXOld, activeXNew );
        InnerProducts( activeX, activeXNew, HDiagList );
        ColumnSubtractions( HDiagList, activeX, activeXNew );
        ColumnNorms( activeXNew, HSubdiagList );
        activeXOld = activeX;
        activeX    = activeXNew; 
        InvBetaScale( HSubdiagList, activeX );
        if( progress )
            subtimer.Start();
        ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts );
        if( progress )
            std::cout << "  Ritz computations: " << subtimer.Stop() 
                      << " seconds" << std::endl;

        auto activeConverged = 
            FindConverged( lastActiveEsts, activeEsts, activeItCounts, tol );
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
            ( HDiagList, HSubdiagList, activeShifts, activePreimage, activeXOld,
              activeX, activeEsts, activeConverged, activeItCounts, progress );

        lastActiveEsts = activeEsts;
    } 
    if( numDone != numShifts )
        RuntimeError("Two-norm estimates did not converge in time");

    invNorms = estimates;
    if( deflate )
        RestoreOrdering( preimage, invNorms, itCounts );

    return itCounts;
}

template<typename Real>
inline Matrix<Int>
HessenbergLanczos
( const Matrix<Complex<Real> >& H, const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, 
  bool deflate=true, Int maxIts=1000, Real tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::HessenbergLanczos"))
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

    // Simultaneously run Lanczos for various shifts
    Matrix<C> XOld, X, XNew;
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    std::vector<std::vector<Real>> HDiagList( numShifts ),
                                   HSubdiagList( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        HDiagList[j].reserve( HCapacityInit );
        HSubdiagList[j].reserve( HCapacityInit-1 );
    }

    Timer timer, subtimer;
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
        auto activeXOld = View( XOld, 0, 0, n, numActive );
        auto activeX    = View( X,    0, 0, n, numActive );
        auto activeXNew = View( XNew, 0, 0, n, numActive );
        if( deflate )
            activePreimage = View( preimage, 0, 0, numActive, 1 );

        if( progress )
            timer.Start();
        activeXNew = activeX;
        if( progress )
            subtimer.Start();
        Conjugate( activeShifts, activeShiftsConj );
        MultiShiftHessSolve
        ( UPPER, NORMAL, C(1), H, activeShifts, activeXNew );
        MultiShiftHessSolve
        ( LOWER, NORMAL, C(1), HAdj, activeShiftsConj, activeXNew );
        if( progress )
        {
            const double msTime = subtimer.Stop();
            const double gflops = (32.*n*n*numShifts)/(msTime*1.e9);
            std::cout << "  MultiShiftHessSolve's: " << msTime << " seconds, "
                      << gflops << " GFlops" << std::endl;
        }
        ColumnSubtractions( HSubdiagList, activeXOld, activeXNew );
        InnerProducts( activeX, activeXNew, HDiagList );
        ColumnSubtractions( HDiagList, activeX, activeXNew );
        ColumnNorms( activeXNew, HSubdiagList );
        activeXOld = activeX;
        activeX    = activeXNew; 
        InvBetaScale( HSubdiagList, activeX );
        if( progress )
            subtimer.Start();
        ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts );
        if( progress )
            std::cout << "  Ritz computations: " << subtimer.Stop() 
                      << " seconds" << std::endl;

        auto activeConverged = 
            FindConverged( lastActiveEsts, activeEsts, activeItCounts, tol );
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
            ( HDiagList, HSubdiagList, activeShifts, activePreimage, activeXOld,
              activeX, activeEsts, activeConverged, activeItCounts, progress );

        lastActiveEsts = activeEsts;
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
TriangularLanczos
( const DistMatrix<Complex<Real>        >& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
  bool deflate=true, Int maxIts=1000, Real tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::TriangularLanczos"))
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
        for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
        {
            const Int i = preimage.GlobalRow(iLoc);
            preimage.SetLocal( iLoc, 0, i );
        }
    }

    // Simultaneously run Lanczos for various shifts
    DistMatrix<C> XOld(g), X(g), XNew(g);
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    std::vector<std::vector<Real>> HDiagList( X.LocalWidth() ),
                                   HSubdiagList( X.LocalWidth() );
    for( Int j=0; j<X.LocalWidth(); ++j )
    {
        HDiagList[j].reserve( HCapacityInit );
        HSubdiagList[j].reserve( HCapacityInit-1 );
    }

    Timer timer, subtimer;
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
        auto activeXOld = View( XOld, 0, 0, n, numActive );
        auto activeX    = View( X,    0, 0, n, numActive );
        auto activeXNew = View( XNew, 0, 0, n, numActive );
        if( deflate )
            activePreimage = View( preimage, 0, 0, numActive, 1 );

        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                timer.Start();
        }
        activeXNew = activeX;
        if( progress )
        { 
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                subtimer.Start();
        }
        MultiShiftTrsm
        ( LEFT, UPPER, NORMAL, C(1), U, activeShifts, activeXNew );
        MultiShiftTrsm
        ( LEFT, UPPER, ADJOINT, C(1), U, activeShifts, activeXNew );
        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
            {
                const double msTime = subtimer.Stop();
                const double gflops = (8.*n*n*numShifts)/(msTime*1.e9);
                std::cout << "  MultiShiftTrsm's: " << msTime << " seconds, "
                          << gflops << " GFlops" << std::endl;
            }
        }
        ColumnSubtractions( HSubdiagList, activeXOld, activeXNew );
        InnerProducts( activeX, activeXNew, HDiagList );
        ColumnSubtractions( HDiagList, activeX, activeXNew );
        ColumnNorms( activeXNew, HSubdiagList );
        activeXOld = activeX;
        activeX    = activeXNew;
        InvBetaScale( HSubdiagList, activeX );
        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                subtimer.Start();
        }
        ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts );
        if( progress )
        {
            mpi::Barrier( U.Grid().Comm() );
            if( U.Grid().Rank() == 0 )
                std::cout << "  Ritz computations: " << subtimer.Stop() 
                          << " seconds" << std::endl;
        }

        auto activeConverged =
            FindConverged( lastActiveEsts, activeEsts, activeItCounts, tol );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
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
HessenbergLanczos
( const DistMatrix<Complex<Real>        >& H, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
  bool deflate=true, Int maxIts=1000, Real tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::HessenbergLanczos"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = H.Height();
    const Int numShifts = shifts.Height();
    const Grid& g = H.Grid();
    if( deflate && H.Grid().Rank() == 0 ) 
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
        for( Int iLoc=0; iLoc<numLocShifts; ++iLoc )
        {
            const Int i = preimage.GlobalRow(iLoc);
            preimage.SetLocal( iLoc, 0, i );
        }
    }

    DistMatrix<C,VC,STAR> H_VC_STAR( H );
    DistMatrix<C,VC,STAR> HAdj_VC_STAR( H.Grid() );
    Adjoint( H, HAdj_VC_STAR );
    DistMatrix<C,STAR,VR> activeXNew_STAR_VR( H.Grid() );
    DistMatrix<C,VR,STAR> activeShiftsConj( H.Grid() );

    // Simultaneously run Lanczos for various shifts
    DistMatrix<C> XOld(g), X(g), XNew(g);
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    std::vector<std::vector<Real>> HDiagList( X.LocalWidth() ),
                                   HSubdiagList( X.LocalWidth() );
    for( Int j=0; j<X.LocalWidth(); ++j )
    {
        HDiagList[j].reserve( HCapacityInit );
        HSubdiagList[j].reserve( HCapacityInit-1 );
    }

    Timer timer, subtimer;
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
        auto activeXOld = View( XOld, 0, 0, n, numActive );
        auto activeX    = View( X,    0, 0, n, numActive );
        auto activeXNew = View( XNew, 0, 0, n, numActive );
        if( deflate )
            activePreimage = View( preimage, 0, 0, numActive, 1 );

        if( progress )
        {
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
                timer.Start();
        }
        activeXNew = activeX;
        if( progress )
        { 
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
                subtimer.Start();
        }
        // NOTE: This redistribution sequence might not be necessary
        activeXNew_STAR_VR = activeXNew;
        Conjugate( activeShifts, activeShiftsConj );
        MultiShiftHessSolve
        ( UPPER, NORMAL, C(1), H_VC_STAR, activeShifts, 
          activeXNew_STAR_VR );
        MultiShiftHessSolve
        ( LOWER, NORMAL, C(1), HAdj_VC_STAR, activeShiftsConj, 
          activeXNew_STAR_VR );
        activeXNew = activeXNew_STAR_VR;
        if( progress )
        {
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
            {
                const double msTime = subtimer.Stop();
                const double gflops = (32.*n*n*numShifts)/(msTime*1.e9);
                std::cout << "  MultiShiftHessSolve's: " << msTime 
                          << " seconds, " << gflops << " GFlops" << std::endl;
            }
        }
        ColumnSubtractions( HSubdiagList, activeXOld, activeXNew );
        InnerProducts( activeX, activeXNew, HDiagList );
        ColumnSubtractions( HDiagList, activeX, activeXNew );
        ColumnNorms( activeXNew, HSubdiagList );
        activeXOld = activeX;
        activeX    = activeXNew;
        InvBetaScale( HSubdiagList, activeX );
        if( progress )
        {
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
                subtimer.Start();
        }
        ComputeNewEstimates( HDiagList, HSubdiagList, activeEsts );
        if( progress )
        {
            mpi::Barrier( H.Grid().Comm() );
            if( H.Grid().Rank() == 0 )
                std::cout << "  Ritz computations: " << subtimer.Stop() 
                          << " seconds" << std::endl;
        }

        auto activeConverged =
            FindConverged( lastActiveEsts, activeEsts, activeItCounts, tol );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
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

#endif // ifndef ELEM_PSEUDOSPECTRUM_LANCZOS_HPP
