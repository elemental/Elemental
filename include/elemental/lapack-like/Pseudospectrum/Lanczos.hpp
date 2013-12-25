/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_PSEUDOSPECTRUM_LANCZOS_HPP
#define ELEM_LAPACK_PSEUDOSPECTRUM_LANCZOS_HPP

#include "elemental/lapack-like/Norm/Zero.hpp"

#include "./ShiftedTrsm.hpp"

namespace elem {
namespace pspec {

const Int HCapacityInit = 10;

template<typename F>
inline void
ColumnSubtractions
( const std::vector<std::vector<BASE(F)> >& diags,
  const Matrix<F>& X, 
        Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnSubtractions"))
    typedef Base<F> Real;
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

template<typename F>
inline void
ColumnSubtractions
( const std::vector<std::vector<BASE(F)> >& diags,
  const DistMatrix<F>& X, 
        DistMatrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnSubtractions"))
    ColumnSubtractions( diags, X.LockedMatrix(), Y.Matrix() );
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
( const DistMatrix<F>& X,
  const DistMatrix<F>& Y,
        std::vector<std::vector<BASE(F)> >& diags )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    const Int mLocal = X.LocalHeight();
    const Int numLocShifts = X.LocalWidth();
    std::vector<F> innerProds( numLocShifts );
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        innerProds[jLoc] = 
            blas::Dot( mLocal, X.LockedBuffer(0,jLoc), 1, 
                               Y.LockedBuffer(0,jLoc), 1 );
    mpi::AllReduce( innerProds.data(), numLocShifts, mpi::SUM, X.ColComm() );
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        diags[jLoc].push_back( RealPart(innerProds[jLoc]) );
}

template<typename F>
inline void
ColumnNorms
( const Matrix<F>& X,
        std::vector<std::vector<BASE(F)> >& diags )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    typedef Base<F> Real;
    const Int numShifts = X.Width();
    const Int m = X.Height();
    for( Int j=0; j<numShifts; ++j )
    {
        const Real beta = blas::Nrm2( m, X.LockedBuffer(0,j), 1 );
        diags[j].push_back( beta );
    }
}

template<typename F>
inline void
ColumnNorms
( const DistMatrix<F>& X,
        std::vector<std::vector<BASE(F)> >& diags )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ColumnNorms"))
    typedef Base<F> Real;
    const Int mLocal = X.LocalHeight();
    const Int numLocShifts = X.LocalWidth();
    std::vector<Real> normsSquared( numLocShifts );
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
    {
        const Real beta = blas::Nrm2( mLocal, X.LockedBuffer(0,jLoc), 1 );
        normsSquared[jLoc] = beta*beta;
    }
    mpi::AllReduce( normsSquared.data(), numLocShifts, mpi::SUM, X.ColComm() );
    for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        diags[jLoc].push_back( Sqrt(normsSquared[jLoc]) );
}

template<typename F>
inline void
InvBetaScale
( const std::vector<std::vector<BASE(F)> >& HSubdiags,
        Matrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InvBetaScale"))
    typedef Base<F> Real;
    const Int numShifts = Y.Width();
    if( numShifts == 0 )
        return;
    const Int m = Y.Height();
    const Int n = HSubdiags[0].size();
    for( Int j=0; j<numShifts; ++j )
    {
        const Real beta = HSubdiags[j][n-1];
        blas::Scal( m, F(1)/beta, Y.Buffer(0,j), 1 );
    }
}

template<typename F>
inline void
InvBetaScale
( const std::vector<std::vector<BASE(F)> >& HSubdiags,
        DistMatrix<F>& Y )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::InvBetaScale"))
    InvBetaScale( HSubdiags, Y.Matrix() );
}

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real> >& HDiags, 
  const std::vector<std::vector<Real> >& HSubdiags,
  Matrix<Real>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    const Int numProbs = activeEsts.Height();
    if( numProbs == 0 )
        return;
    const Int n = HDiags[0].size();
    std::vector<Real> HDiag, HSubdiag, w(n);
    for( Int j=0; j<numProbs; ++j )
    {
        HDiag = HDiags[j]; 
        HSubdiag = HSubdiags[j];
        lapack::SymmetricTridiagonalEig     
        ( 'N', 'I', n, HDiag.data(), HSubdiag.data(), 0, 0, n, n, 0, w.data(),
          0, 1 );
        activeEsts.Set( j, 0, Sqrt(w[0]) );
    }
}

template<typename Real>
inline void
ComputeNewEstimates
( const std::vector<std::vector<Real> >& HDiags, 
  const std::vector<std::vector<Real> >& HSubdiags,
  DistMatrix<Real,MR,STAR>& activeEsts )
{
    DEBUG_ONLY(CallStackEntry cse("pspec::ComputeNewEstimates"))
    ComputeNewEstimates( HDiags, HSubdiags, activeEsts.Matrix() );
}

template<typename Real>
inline void
Deflate
( std::vector<std::vector<Real> >& HDiags,
  std::vector<std::vector<Real> >& HSubdiags,
  Matrix<Complex<Real> >& activeShifts, 
  Matrix<Int           >& activePreimage,
  Matrix<Complex<Real> >& activeXOld,
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
                std::swap( HDiags[swapFrom], HDiags[swapTo] );
                std::swap( HSubdiags[swapFrom], HSubdiags[swapTo] );
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
}

template<typename Real>
inline void
Deflate
( std::vector<std::vector<Real> >& HDiags,
  std::vector<std::vector<Real> >& HSubdiags,
  DistMatrix<Complex<Real>,VR,STAR>& activeShifts,
  DistMatrix<Int,          VR,STAR>& activePreimage,
  DistMatrix<Complex<Real>        >& activeXOld,
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
    DistMatrix<Complex<Real>,VC,STAR> XOldCopy( activeXOld ), XCopy( activeX );

    const Int n = ( activeX.LocalWidth()>0 ? HDiags[0].size() : 0 );
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
                        if( HDiags[localFrom].size() != n )
                            LogicError("Invalid HDiags size");
                        if( HDiags[localTo].size() != n )
                            LogicError("Invalid HDiags size");
                        if( HSubdiags[localFrom].size() != n )
                            LogicError("Invalid HSubdiags size");
                        if( HSubdiags[localTo].size() != n )
                            LogicError("Invalid HSubdiags size");
                    )
                    std::swap( HDiags[localFrom], HDiags[localTo] );
                    std::swap( HSubdiags[localFrom], HSubdiags[localTo] );
                }
                else if( activeX.IsLocalCol(swapFrom) )
                {
                    const Int localFrom = activeX.LocalCol(swapFrom);
                    DEBUG_ONLY(
                        if( HDiags[localFrom].size() != n )
                            LogicError("Invalid HDiags size");
                        if( HSubdiags[localFrom].size() != n )
                            LogicError("Invalid HSubdiags size");
                    )
                    const Int partner = activeX.ColOwner(swapTo);
                    mpi::TaggedSendRecv
                    ( HDiags[localFrom].data(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                    mpi::TaggedSendRecv
                    ( HSubdiags[localFrom].data(), n, 
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                }
                else if( activeX.IsLocalCol(swapTo) )
                {
                    const Int localTo = activeX.LocalCol(swapTo);
                    DEBUG_ONLY(
                        if( HDiags[localTo].size() != n )
                            LogicError("Invalid HDiags size");
                        if( HSubdiags[localTo].size() != n )
                            LogicError("Invalid HSubdiags size");
                    )
                    const Int partner = activeX.ColOwner(swapFrom);
                    mpi::TaggedSendRecv
                    ( HDiags[localTo].data(), n,
                      partner, swapFrom, partner, swapFrom, activeX.RowComm() );
                    mpi::TaggedSendRecv
                    ( HSubdiags[localTo].data(), n, 
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
    const int numShifts = shifts.Height();

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

    // Simultaneously run Lanczos for various shifts
    Matrix<C> XOld, X, XNew;
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    std::vector<std::vector<Real>> HDiags( numShifts ),
                                   HSubdiags( numShifts );
    for( Int j=0; j<numShifts; ++j )
    {
        HDiags[j].reserve( HCapacityInit );
        HSubdiags[j].reserve( HCapacityInit-1 );
    }

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

        activeXNew = activeX;
        ShiftedTrsmLUN( U, activeShifts, activeXNew );
        ShiftedTrsmLUT( U, activeShifts, activeXNew );
        ColumnSubtractions( HSubdiags, activeXOld, activeXNew );
        InnerProducts( activeXNew, activeX, HDiags );
        ColumnSubtractions( HDiags, activeX, activeXNew );
        ColumnNorms( activeXNew, HSubdiags );
        activeXOld = activeX;
        activeX    = activeXNew; 
        InvBetaScale( HSubdiags, activeX );
        ComputeNewEstimates( HDiags, HSubdiags, activeEsts );

        auto activeConverged = 
            FindConverged( lastActiveEsts, activeEsts, activeItCounts, tol*n );
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
            ( HDiags, HSubdiags, activeShifts, activePreimage, activeXOld, 
              activeX, activeEsts, activeConverged, activeItCounts );

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
    const int numShifts = shifts.Height();
    const Grid& g = U.Grid();
    if( deflate && mpi::WorldRank() == 0 ) 
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
        preimage.ResizeTo( numShifts, 1 );
        const Int numLocShifts = preimage.LocalHeight();
        for( Int jLoc=0; jLoc<numLocShifts; ++jLoc )
        {
            const Int j = preimage.ColShift() + jLoc*preimage.ColStride();
            preimage.SetLocal( jLoc, 0, j );
        }
    }

    // Simultaneously run Lanczos for various shifts
    DistMatrix<C> XOld(g), X(g), XNew(g);
    Zeros( XOld, n, numShifts );
    Gaussian( X, n, numShifts );
    FixColumns( X );
    Zeros( XNew, n, numShifts );
    std::vector<std::vector<Real>> HDiags( X.LocalWidth() ),
                                   HSubdiags( X.LocalWidth() );
    for( Int j=0; j<X.LocalWidth(); ++j )
    {
        HDiags[j].reserve( HCapacityInit );
        HSubdiags[j].reserve( HCapacityInit-1 );
    }

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

        activeXNew = activeX;
        ShiftedTrsmLUN( U, activeShifts, activeXNew );
        ShiftedTrsmLUT( U, activeShifts, activeXNew );
        ColumnSubtractions( HSubdiags, activeXOld, activeXNew );
        InnerProducts( activeXNew, activeX, HDiags );
        ColumnSubtractions( HDiags, activeX, activeXNew );
        ColumnNorms( activeXNew, HSubdiags );
        activeXOld = activeX;
        activeX    = activeXNew;
        InvBetaScale( HSubdiags, activeX );
        ComputeNewEstimates( HDiags, HSubdiags, activeEsts );

        auto activeConverged =
            FindConverged( lastActiveEsts, activeEsts, activeItCounts, tol*n );
        const Int numActiveDone = ZeroNorm( activeConverged );
        if( deflate )
            numDone += numActiveDone;
        else
            numDone = numActiveDone;
        if( progress && mpi::WorldRank() == 0 )
            std::cout << numDone << " of " << numShifts << " converged"
                      << std::endl;

        ++numIts;
        if( numIts >= maxIts )
            break;

        if( numDone == numShifts )
            break;
        else if( deflate )
            Deflate
            ( HDiags, HSubdiags, activeShifts, activePreimage, activeXOld, 
              activeX, activeEsts, activeConverged, activeItCounts );

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

#endif // ifndef ELEM_LAPACK_PSEUDOSPECTRUM_LANCZOS_HPP
