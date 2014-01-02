/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_PSEUDOSPECTRUM_POWER_HPP
#define ELEM_LAPACK_PSEUDOSPECTRUM_POWER_HPP

#include "elemental/lapack-like/Norm/Zero.hpp"
#include "elemental/matrices/Ones.hpp"

#include "./ShiftedTrsm.hpp"

namespace elem {
namespace pspec {

template<typename F>
inline Base<F> NormCap()
{ return Base<F>(1)/lapack::MachineEpsilon<Base<F>>(); }

template<typename F>
inline void
ColumnNorms( const Matrix<F>& X, Matrix<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    typedef Base<F> Real;
    const Int m = X.Height();
    const Int n = X.Width();
    norms.ResizeTo( n, 1 );
    for( Int j=0; j<n; ++j )
    {
        const Real alpha = blas::Nrm2( m, X.LockedBuffer(0,j), 1 );
        norms.Set( j, 0, alpha );
    }
}

template<typename F>
inline void
ColumnNorms( const DistMatrix<F>& X, DistMatrix<BASE(F),MR,STAR>& norms )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::ColumnNorms");
        if( X.RowAlign() != norms.ColAlign() )
            LogicError("Invalid norms alignment");
    )
    typedef Base<F> Real;
    const Int n = X.Width();
    const Int mLocal = X.LocalHeight();
    const Int nLocal = X.LocalWidth();
    const Grid& g = X.Grid();

    // TODO: Switch to more stable parallel norm computation using scaling
    norms.ResizeTo( n, 1 ); 
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Base<F> localNorm = blas::Nrm2(mLocal,X.LockedBuffer(0,jLoc),1);
        norms.SetLocal( jLoc, 0, localNorm*localNorm );
    }

    mpi::AllReduce( norms.Buffer(), nLocal, mpi::SUM, X.ColComm() );
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Real alpha = norms.GetLocal(jLoc,0);
        norms.SetLocal( jLoc, 0, Sqrt(alpha) );
    }
}

template<typename F>
inline void
FixColumns( Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FixColumns"))
    typedef Base<F> Real;
    Matrix<Real> norms;
    ColumnNorms( X, norms );
    const Int m = X.Height();
    const Int n = X.Width();
    for( Int j=0; j<n; ++j )
    {
        auto x = View( X, 0, j, m, 1 );
        Real norm = norms.Get(j,0);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        Scale( Real(1)/norm, x );
    }
}

template<typename F>
inline void
FixColumns( DistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FixColumns"))
    typedef Base<F> Real;
    DistMatrix<Real,MR,STAR> norms( X.Grid() );
    ColumnNorms( X, norms );
    const Int m = X.Height();
    const Int nLocal = X.LocalWidth();
    for( Int jLoc=0; jLoc<nLocal; ++jLoc )
    {
        const Int j = X.RowShift() + jLoc*X.RowStride();
        auto x = View( X, 0, j, m, 1 );
        Real norm = norms.GetLocal(jLoc,0);
        if( norm == Real(0) )
        {
            MakeGaussian( x );
            norm = FrobeniusNorm( x );
        }
        Scale( Real(1)/norm, x );
    }
}

template<typename Real>
inline void CapEstimates( Matrix<Real>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::CapEstimates"))
    const Real normCap = NormCap<Real>();
    const Int n = activeEsts.Height();
    for( Int j=0; j<n; ++j )
    {
        Real alpha = activeEsts.Get(j,0);
        if( std::isnan(alpha) || alpha >= normCap )
            alpha = normCap;
        activeEsts.Set( j, 0, alpha );
    }
}

template<typename Real>
inline void CapEstimates( DistMatrix<Real,MR,STAR>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::CapEstimates"))
    CapEstimates( activeEsts.Matrix() );
}

template<typename Real>
inline Matrix<Int>
FindConverged
( const Matrix<Real>& lastActiveEsts, 
  const Matrix<Real>& activeEsts,
        Matrix<Int >& activeItCounts,
        Real maxDiff )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::FindConverged"))
    const Real normCap = NormCap<Real>();

    const Int numActiveShifts=activeEsts.Height();
    Matrix<Int> activeConverged;
    Zeros( activeConverged, numActiveShifts, 1 );

    for( Int j=0; j<numActiveShifts; ++j )
    {
        const Real lastEst = lastActiveEsts.Get(j,0);
        const Real currEst = activeEsts.Get(j,0);
        bool converged = false;
        if( currEst >= normCap )
            converged = true;
        else if( Abs(currEst) > 0 )
            converged = (Abs(lastEst-currEst)/Abs(currEst) <= maxDiff);

        if( converged )
            activeConverged.Set( j, 0, 1 );
        else 
            activeItCounts.Update( j, 0, 1 );
    }
    return activeConverged;
}

template<typename Real>
inline DistMatrix<Int,MR,STAR>
FindConverged
( const DistMatrix<Real,MR,STAR>& lastActiveEsts,
  const DistMatrix<Real,MR,STAR>& activeEsts,
        DistMatrix<Int, VR,STAR>& activeItCounts,
        Real maxDiff )
{
    DEBUG_ONLY(
        CallStackEntry cse("pspec::NumConverged");
        if( activeItCounts.ColAlign()%activeEsts.ColStride() !=
            activeEsts.ColAlign() )
            LogicError("Invalid column alignment");
    )
    const Real normCap = NormCap<Real>();

    DistMatrix<Int,MR,STAR> activeConverged( activeEsts.Grid() );
    activeConverged.AlignWith( activeEsts );
    Zeros( activeConverged, activeEsts.Height(), 1 );

    const Int numLocShifts=activeEsts.LocalHeight();
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
    {
        const Real lastEst = lastActiveEsts.GetLocal(jLoc,0);
        const Real currEst = activeEsts.GetLocal(jLoc,0);
        bool converged = false;
        if( currEst >= normCap )
            converged = true;
        else if( Abs(currEst) > 0 )
            converged = (Abs(lastEst-currEst)/Abs(currEst) <= maxDiff);

        if( converged )
            activeConverged.SetLocal( jLoc, 0, 1 );
        else 
        {
            const Int j = activeEsts.ColShift()+jLoc*activeEsts.ColStride();
            activeItCounts.Update( j, 0, 1 );
        }
    }

    return activeConverged;
}

template<typename Real>
inline void
Deflate
( Matrix<Complex<Real> >& activeShifts, 
  Matrix<Int           >& activePreimage,
  Matrix<Complex<Real> >& activeX,
  Matrix<Real          >& activeEsts, 
  Matrix<Int           >& activeConverged,
  Matrix<Int           >& activeItCounts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Deflate"))
    const Int numActive = activeX.Width(); 
    Int swapTo = numActive-1;
    for( Int swapFrom=numActive-1; swapFrom>=0; --swapFrom )
    {
        if( activeConverged.Get(swapFrom,0) )
        {
            if( swapTo != swapFrom )
            {
                RowSwap( activeShifts, swapFrom, swapTo );
                RowSwap( activePreimage, swapFrom, swapTo );
                RowSwap( activeEsts, swapFrom, swapTo );
                RowSwap( activeItCounts, swapFrom, swapTo );
                ColumnSwap( activeX, swapFrom, swapTo );
            }
            --swapTo;
        }
    }
}

template<typename Real>
inline void
Deflate
( DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
  DistMatrix<Complex<Real>        >& activeX,
  DistMatrix<Real,         MR,STAR>& activeEsts,
  DistMatrix<Int,          MR,STAR>& activeConverged,
  DistMatrix<Int,          VR,STAR>& activeItCounts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::Deflate"))
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
                RowSwap( shiftsCopy, swapFrom, swapTo );
                RowSwap( preimageCopy, swapFrom, swapTo );
                RowSwap( estimatesCopy, swapFrom, swapTo );
                RowSwap( itCountsCopy, swapFrom, swapTo );
                ColumnSwap( XCopy, swapFrom, swapTo );
            }
            --swapTo;
        }
    }

    activeShifts   = shiftsCopy;
    activePreimage = preimageCopy;
    activeEsts     = estimatesCopy;
    activeItCounts = itCountsCopy;
    activeX        = XCopy;
}

template<typename Real>
inline void
RestoreOrdering
( const Matrix<Int>& preimage, Matrix<Real>& invNorms, Matrix<Int>& itCounts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    auto invNormsCopy = invNorms;
    auto itCountsCopy = itCounts;
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimage.Get(j,0);
        invNorms.Set( dest, 0, invNormsCopy.Get(j,0) );
        itCounts.Set( dest, 0, itCountsCopy.Get(j,0) );
    }
}

template<typename Real>
inline void
RestoreOrdering
( const DistMatrix<Int, VR,STAR>& preimage, 
        DistMatrix<Real,VR,STAR>& invNorms,
        DistMatrix<Int, VR,STAR>& itCounts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::RestoreOrdering"))
    DistMatrix<Int, STAR,STAR> preimageCopy( preimage );
    DistMatrix<Real,STAR,STAR> invNormsCopy( invNorms );
    DistMatrix<Int, STAR,STAR> itCountsCopy( itCounts );
    const Int numShifts = preimage.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Int dest = preimageCopy.Get(j,0);
        invNorms.Set( dest, 0, invNormsCopy.Get(j,0) );
        itCounts.Set( dest, 0, itCountsCopy.Get(j,0) );
    }
}

template<typename Real>
inline Matrix<Int>
TriangularPower
( const Matrix<Complex<Real> >& U, const Matrix<Complex<Real> >& shifts, 
  Matrix<Real>& invNorms, 
  bool deflate=true, Int maxIts=1000, Real tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::TriangularPower"))
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
        preimage.ResizeTo( numShifts, 1 );
        for( Int j=0; j<numShifts; ++j )
            preimage.Set( j, 0, j );
    }

    // Simultaneously run inverse iteration for various shifts
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

        ShiftedTrsmLUN( U, activeShifts, activeX );
        FixColumns( activeX );
        ShiftedTrsmLUT( U, activeShifts, activeX );
        ColumnNorms( activeX, activeEsts );
        CapEstimates( activeEsts );

        auto activeConverged = 
            FindConverged( lastActiveEsts, activeEsts, activeItCounts, tol );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        if( progress )
            std::cout << numDone << " of " << numShifts << " converged"
                      << std::endl;

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate )
            Deflate
            ( activeShifts, activePreimage, activeX, activeEsts,
              activeConverged, activeItCounts );

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
TriangularPower
( const DistMatrix<Complex<Real>        >& U, 
  const DistMatrix<Complex<Real>,VR,STAR>& shifts, 
        DistMatrix<Real,         VR,STAR>& invNorms, 
  bool deflate=true, Int maxIts=1000, Real tol=1e-6, bool progress=false )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::TriangularPower"))
    using namespace pspec;
    typedef Complex<Real> C;
    const Int n = U.Height();
    const Int numShifts = shifts.Height();
    const Grid& g = U.Grid();

    // Keep track of the number of iterations per shift
    DistMatrix<Int,VR,STAR> itCounts(g);
    Ones( itCounts, numShifts, 1 );

    // Keep track of the pivoting history if deflation is requested
    DistMatrix<Int,VR,STAR> preimage(g);
    DistMatrix<C,  VR,STAR> pivShifts( shifts );
    if( deflate )
    {
        preimage.AlignWith( shifts );
        preimage.ResizeTo( numShifts, 1 );
        const Int numLocShifts = preimage.LocalHeight();
        for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        {
            const Int j = preimage.ColShift() + jLoc*preimage.ColStride();
            preimage.SetLocal( jLoc, 0, j );
        }
    }

    // Simultaneously run inverse iteration for various shifts
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

        ShiftedTrsmLUN( U, activeShifts, activeX );
        FixColumns( activeX );
        ShiftedTrsmLUT( U, activeShifts, activeX );
        ColumnNorms( activeX, activeEsts );
        CapEstimates( activeEsts );

        auto activeConverged =
            FindConverged( lastActiveEsts, activeEsts, activeItCounts, tol );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        if( progress && U.Grid().Rank() == 0 )
            std::cout << numDone << " of " << numShifts << " converged"
                      << std::endl;

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate )
            Deflate
            ( activeShifts, activePreimage, activeX, activeEsts,
              activeConverged, activeItCounts );

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

#endif // ifndef ELEM_LAPACK_PSEUDOSPECTRUM_POWER_HPP
