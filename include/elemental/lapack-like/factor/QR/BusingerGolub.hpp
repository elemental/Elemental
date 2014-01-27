/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_QR_BUSINGERGOLUB_HPP
#define ELEM_QR_BUSINGERGOLUB_HPP

#include ELEM_GEMV_INC
#include ELEM_GER_INC

#include ELEM_REFLECTOR_INC

#include ELEM_ZEROS_INC

#include <algorithm>

namespace elem {
namespace qr {

template<typename F>
inline BASE(F)
ColumnNorms( const Matrix<F>& A, std::vector<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ColumnNorms"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    Real maxNorm = 0;
    norms.resize( n );
    for( Int j=0; j<n; ++j )
    {
        norms[j] = blas::Nrm2( m, A.LockedBuffer(0,j), 1 );
        maxNorm = std::max( maxNorm, norms[j] );
    }
    return maxNorm;
}

template<typename Real>
inline ValueInt<Real>
FindPivot( const std::vector<Real>& norms, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("qr::FindPivot"))
    const auto maxNorm = std::max_element( norms.begin()+col, norms.end() );
    ValueInt<Real> pivot;
    pivot.value = *maxNorm;
    pivot.index = maxNorm - norms.begin();
    return pivot;
}

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Int>& p,
  Int maxSteps, BASE(F) tol, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
        if( maxSteps > Min(m,n) )
            LogicError("Too many steps requested");
    )
    p.Resize( maxSteps, 1 );
    t.Resize( maxSteps, 1 );

    Matrix<F> z;

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the 
    // updated quantities are no longer accurate.
    std::vector<Real> origNorms;
    const Real maxOrigNorm = ColumnNorms( A, origNorms );
    std::vector<Real> norms = origNorms;
    const Real updateTol = Sqrt(lapack::MachineEpsilon<Real>());

    for( Int k=0; k<maxSteps; ++k )
    {
        auto alpha11   = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12       = ViewRange( A, k,   k+1, k+1, n   );
        auto a21       = ViewRange( A, k+1, k,   m,   k+1 );
        auto aLeftCol  = ViewRange( A, k,   k,   m,   k+1 );
        auto ARightPan = ViewRange( A, k,   k+1, m,   n   );

        // Find the next column pivot
        const ValueInt<Real> pivot = FindPivot( norms, k );
        if( pivot.value <= tol*maxOrigNorm )
        {
            p.Resize( k, 1 );
            t.Resize( k, 1 );
            break;
        }
        p.Set( k, 0, pivot.index );
 
        // Perform the swap
        const Int jPiv = pivot.index;
        if( jPiv != k )
        {
            blas::Swap( m, A.Buffer(0,k), 1, A.Buffer(0,jPiv), 1 );
            norms[jPiv] = norms[k];
            origNorms[jPiv] = origNorms[k];
        }

        // Compute and apply the Householder reflector for this column
        const F tau = Reflector( alpha11, a21 );
        t.Set( k, 0, tau );
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Zeros( z, ARightPan.Width(), 1 );
        Gemv( ADJOINT, F(1), ARightPan, aLeftCol, F(0), z );
        Ger( -Conj(tau), aLeftCol, z, ARightPan );
        alpha11.Set(0,0,alpha);

        // Update the column norm estimates in the same manner as LAWN 176
        for( Int j=k+1; j<n; ++j )
        {
            if( norms[j] != Real(0) )
            {
                Real gamma = Abs(A.Get(k,j)) / norms[j];
                gamma = std::max( Real(0), (Real(1)-gamma)*(Real(1)+gamma) );

                const Real ratio = norms[j] / origNorms[j];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol || alwaysRecompute )
                {
                    norms[j] = blas::Nrm2( m-(k+1), A.Buffer(k+1,j), 1 );
                    origNorms[j] = norms[j];
                }
                else
                    norms[j] *= Sqrt(gamma);
            }
        }
    }
}

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<Int>& p,
  Int maxSteps, BASE(F) tol, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    Matrix<F> t;
    BusingerGolub( A, t, p, maxSteps, tol, alwaysRecompute );
}

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Int>& p,
  Int numSteps, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    // Use a tolerance of -1 so that we always perform numSteps iterations
    BusingerGolub( A, t, p, numSteps, BASE(F)(-1), alwaysRecompute );
}

// If we don't need 't' from the above routine
template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<Int>& p, Int numSteps, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    Matrix<F> t;
    BusingerGolub( A, t, p, numSteps, alwaysRecompute );
}

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Int>& p, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    const Int numSteps = Min(A.Height(),A.Width());
    BusingerGolub( A, t, p, numSteps, alwaysRecompute );
}

// If we don't need 't' from the above routine
template<typename F> 
inline void
BusingerGolub( Matrix<F>& A, Matrix<Int>& p, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    Matrix<F> t;
    BusingerGolub( A, t, p, alwaysRecompute );
}

template<typename F>
inline ValueInt<BASE(F)>
FindColumnPivot
( const DistMatrix<F>& A, const std::vector<BASE(F)>& norms, Int col )
{
    DEBUG_ONLY(CallStackEntry cse("qr::FindColumnPivot"))
    typedef Base<F> Real;
    const Int rowShift = A.RowShift();
    const Int rowStride = A.RowStride();
    const Int localColsBefore = Length( col, rowShift, rowStride );
    const ValueInt<Real> localPivot = FindPivot( norms, localColsBefore );
    ValueInt<Real> pivot;
    pivot.value = localPivot.value;
    pivot.index = rowShift + localPivot.index*rowStride;
    return mpi::AllReduce( pivot, mpi::MaxLocOp<Real>(), A.Grid().RowComm() );
}

template<typename F>
inline BASE(F)
ColumnNorms( const DistMatrix<F>& A, std::vector<BASE(F)>& norms )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ColumnNorms"))
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    mpi::Comm colComm = A.Grid().ColComm();
    mpi::Comm rowComm = A.Grid().RowComm();

    // Carefully perform the local portion of the computation
    std::vector<Real> localScales(localWidth,0), 
                      localScaledSquares(localWidth,1);
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Real alphaAbs = Abs(A.GetLocal(iLoc,jLoc));    
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= localScales[jLoc] )
                {
                    const Real relScale = alphaAbs/localScales[jLoc];
                    localScaledSquares[jLoc] += relScale*relScale;
                }
                else
                {
                    const Real relScale = localScales[jLoc]/alphaAbs;
                    localScaledSquares[jLoc] = 
                        localScaledSquares[jLoc]*relScale*relScale + 1;
                    localScales[jLoc] = alphaAbs;
                }
            }
        }
    }

    // Find the maximum relative scales 
    std::vector<Real> scales(localWidth);
    mpi::AllReduce
    ( localScales.data(), scales.data(), localWidth, mpi::MAX, colComm );

    // Equilibrate the local scaled sums to the maximum scale
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        if( scales[jLoc] != 0 )
        {
            const Real relScale = localScales[jLoc]/scales[jLoc];
            localScaledSquares[jLoc] *= relScale*relScale;
        }
    }

    // Now sum the local contributions (can ignore results where scale is 0)
    std::vector<Real> scaledSquares(localWidth); 
    mpi::AllReduce
    ( localScaledSquares.data(), scaledSquares.data(), localWidth, colComm );

    // Finish the computation
    Real maxLocalNorm = 0;
    norms.resize( localWidth );
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        if( scales[jLoc] != 0 )
            norms[jLoc] = scales[jLoc]*Sqrt(scaledSquares[jLoc]);
        else
            norms[jLoc] = 0;
        maxLocalNorm = std::max( maxLocalNorm, norms[jLoc] );
    }
    return mpi::AllReduce( maxLocalNorm, mpi::MAX, rowComm );
}

template<typename F>
inline void
ReplaceColumnNorms
( const DistMatrix<F>& A, std::vector<Int>& inaccurateNorms, 
  std::vector<BASE(F)>& norms, std::vector<BASE(F)>& origNorms )
{
    DEBUG_ONLY(CallStackEntry cse("qr::ReplaceColumnNorms"))
    typedef Base<F> Real;
    const Int localHeight = A.LocalHeight();
    const Int numInaccurate = inaccurateNorms.size();
    mpi::Comm colComm = A.Grid().ColComm();

    // Carefully perform the local portion of the computation
    std::vector<Real> localScales(numInaccurate,0), 
                      localScaledSquares(numInaccurate,1);
    for( Int s=0; s<numInaccurate; ++s )
    {
        const Int jLoc = inaccurateNorms[s];
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Real alphaAbs = Abs(A.GetLocal(iLoc,jLoc));    
            if( alphaAbs != 0 )
            {
                if( alphaAbs <= localScales[s] )
                {
                    const Real relScale = alphaAbs/localScales[s];
                    localScaledSquares[s] += relScale*relScale;
                }
                else
                {
                    const Real relScale = localScales[s]/alphaAbs;
                    localScaledSquares[s] = 
                        localScaledSquares[s]*relScale*relScale + 1;
                    localScales[s] = alphaAbs;
                }
            }
        }
    }

    // Find the maximum relative scales 
    std::vector<Real> scales(numInaccurate);
    mpi::AllReduce
    ( localScales.data(), scales.data(), numInaccurate, mpi::MAX, colComm );

    // Equilibrate the local scaled sums to the maximum scale
    for( Int s=0; s<numInaccurate; ++s )
    {
        if( scales[s] != 0 )
        {
            const Real relScale = localScales[s]/scales[s];
            localScaledSquares[s] *= relScale*relScale;
        }
    }

    // Now sum the local contributions (can ignore results where scale is 0)
    std::vector<Real> scaledSquares(numInaccurate); 
    mpi::AllReduce
    ( localScaledSquares.data(), scaledSquares.data(), numInaccurate, colComm );

    // Finish the computation
    for( Int s=0; s<numInaccurate; ++s )
    {
        const Int jLoc = inaccurateNorms[s];
        if( scales[s] != 0 )
            norms[jLoc] = scales[s]*Sqrt(scaledSquares[s]);
        else
            norms[jLoc] = 0;
        origNorms[jLoc] = norms[jLoc];
    }
}

template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Int,VR,STAR>& p,
  Int maxSteps, BASE(F) tol, bool alwaysRecompute=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("qr::BusingerGolub");
        if( A.Grid() != p.Grid() || A.Grid() != t.Grid() )
            LogicError("A, t, and p must have the same grid");
    )
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
        if( maxSteps > Min(m,n) )
            LogicError("Too many steps requested");
    )
    t.Resize( maxSteps, 1 );
    p.Resize( maxSteps, 1 );

    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const Int rowAlign = A.RowAlign();
    const Int rowShift = A.RowShift();
    const Int rowStride = A.RowStride();

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the 
    // updated quantities are no longer accurate.
    std::vector<Real> origNorms( nLocal );
    const Real maxOrigNorm = ColumnNorms( A, origNorms );
    std::vector<Real> norms = origNorms;
    const Real updateTol = Sqrt(lapack::MachineEpsilon<Real>());
    std::vector<Int> inaccurateNorms;

    const Grid& g = A.Grid();
    DistMatrix<F> z(g);
    DistMatrix<F,MC,STAR> aLeftCol_MC_STAR(g);
    DistMatrix<F,MR,STAR> z_MR_STAR(g);
    DistMatrix<F,STAR,MR> a12_STAR_MR(g);

    for( Int k=0; k<maxSteps; ++k )
    {
        auto alpha11   = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12       = ViewRange( A, k,   k+1, k+1, n   );
        auto a21       = ViewRange( A, k+1, k,   m,   k+1 );
        auto aLeftCol  = ViewRange( A, k,   k,   m,   k+1 );
        auto ARightPan = ViewRange( A, k,   k+1, m,   n   );

        // Find the next column pivot
        const ValueInt<Real> pivot = FindColumnPivot( A, norms, k );
        if( pivot.value <= tol*maxOrigNorm )
        {
            p.Resize( k, 1 );
            t.Resize( k, 1 );
            break;
        }
        p.Set( k, 0, pivot.index );

        // Perform the swap
        const Int jPiv = pivot.index;
        const Int curOwner = (k+rowAlign) % rowStride;
        const Int pivOwner = (jPiv+rowAlign) % rowStride;
        const bool myCur = ( g.Col() == curOwner );
        const bool myPiv = ( g.Col() == pivOwner );
        if( jPiv != k )
        {
            if( myCur && myPiv )
            {
                const Int kLoc    = (k   -rowShift) / rowStride;
                const Int jPivLoc = (jPiv-rowShift) / rowStride;
                blas::Swap
                ( mLocal, A.Buffer(0,kLoc), 1, A.Buffer(0,jPivLoc), 1 );
                norms[jPivLoc] = norms[kLoc];
                origNorms[jPivLoc] = origNorms[kLoc];
            }
            else if( myCur )
            {
                const Int kLoc = (k-rowShift) / rowStride;
                mpi::SendRecv
                ( A.Buffer(0,kLoc), mLocal, pivOwner, pivOwner, g.RowComm() );
                mpi::Send( norms[kLoc], pivOwner, g.RowComm() );
            }
            else if( myPiv )
            {
                const Int jPivLoc = (jPiv-rowShift) / rowStride;
                mpi::SendRecv
                ( A.Buffer(0,jPivLoc), mLocal, 
                  curOwner, curOwner, g.RowComm() );
                norms[jPivLoc] = mpi::Recv<Real>( curOwner, g.RowComm() );
            }
        }

        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a21 );
        t.Set( k, 0, tau );

        // Apply the Householder reflector
        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlign() &&
                                       g.Col() == alpha11.RowAlign() );
        F alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }
        aLeftCol_MC_STAR.AlignWith( ARightPan );
        aLeftCol_MC_STAR = aLeftCol;
        z_MR_STAR.AlignWith( ARightPan );
        Zeros( z_MR_STAR, ARightPan.Width(), 1 );
        LocalGemv
        ( ADJOINT, F(1), ARightPan, aLeftCol_MC_STAR, F(0), z_MR_STAR );
        z_MR_STAR.SumOverCol();
        Ger
        ( -Conj(tau),
          aLeftCol_MC_STAR.LockedMatrix(),
          z_MR_STAR.LockedMatrix(),
          ARightPan.Matrix() );
        if( myDiagonalEntry )
            alpha11.SetLocal(0,0,alpha);

        // Update the column norm estimates in the same manner as LAWN 176.
        // However, we do so in two steps in order to lower the communication
        // latency:
        //   1) Each process first computes which of its column norms are
        //      too inaccurate and need to be recomputed.
        //   2) Each process communicates within its process column in order
        //      to replace the inaccurate column norms.
        // Step 1: Perform all of the easy updates and mark inaccurate norms
        a12_STAR_MR = a12;
        const Int a12LocalWidth = a12_STAR_MR.LocalWidth();
        const Int a12RowShift = a12_STAR_MR.RowShift();
        inaccurateNorms.resize(0);
        for( Int jLoc12=0; jLoc12<a12LocalWidth; ++jLoc12 )
        {
            const Int j = (k+1) + a12RowShift + jLoc12*rowStride;
            const Int jLoc = (j-rowShift) / rowStride;
            if( norms[jLoc] != Real(0) )
            {
                const Real beta = Abs(a12_STAR_MR.GetLocal(0,jLoc12));
                Real gamma = beta / norms[jLoc];
                gamma = std::max( Real(0), (Real(1)-gamma)*(Real(1)+gamma) );

                const Real ratio = norms[jLoc] / origNorms[jLoc];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol || alwaysRecompute )
                    inaccurateNorms.push_back( jLoc );
                else
                    norms[jLoc] *= Sqrt(gamma);
            }
        }
        // Step 2: Compute the replacement norms and also reset origNorms
        ReplaceColumnNorms( A, inaccurateNorms, norms, origNorms );
    }
}

// If we don't need 't' from the above routine
template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<Int,VR,STAR>& p,
  Int maxSteps, BASE(F) tol, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    DistMatrix<F,MD,STAR> t( A.Grid() );
    BusingerGolub( A, t, p, maxSteps, tol, alwaysRecompute );
}

template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Int,VR,STAR>& p,
  Int numSteps, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    // Use a tolerance of -1 so that we always perform numSteps iterations
    BusingerGolub( A, t, p, numSteps, BASE(F)(-1), alwaysRecompute );
}

// If we don't need 't' from the above routine
template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<Int,VR,STAR>& p,
  Int numSteps, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    DistMatrix<F,MD,STAR> t( A.Grid() );
    BusingerGolub( A, t, p, numSteps, alwaysRecompute );
}

template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Int,VR,STAR>& p,
  bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    const Int numSteps = Min(A.Height(),A.Width());
    BusingerGolub( A, t, p, numSteps, alwaysRecompute );
}

// If we don't need 't' from the above routine
template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<Int,VR,STAR>& p, bool alwaysRecompute=false )
{
    DEBUG_ONLY(CallStackEntry cse("qr::BusingerGolub"))
    DistMatrix<F,MD,STAR> t( A.Grid() );
    BusingerGolub( A, t, p, alwaysRecompute );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_QR_BUSINGERGOLUB_HPP
