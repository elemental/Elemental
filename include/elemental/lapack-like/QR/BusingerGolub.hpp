/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_QR_BUSINGERGOLUB_HPP
#define ELEM_LAPACK_QR_BUSINGERGOLUB_HPP

#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

#include <algorithm>

namespace elem {
namespace qr {

template<typename F>
inline BASE(F)
ColumnNorms( const Matrix<F>& A, std::vector<BASE(F)>& norms )
{
#ifndef RELEASE
    CallStackEntry entry("qr::ColumnNorms");
#endif
    typedef BASE(F) Real;
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
#ifndef RELEASE
    CallStackEntry entry("qr::FindPivot");
#endif
    const Int n = norms.size();
    const Real* maxNorm = std::max_element( &norms[col], &norms[0]+n );
    ValueInt<Real> pivot;
    pivot.value = *maxNorm;
    pivot.index = maxNorm - &norms[0];
    return pivot;
}

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Int>& p,
  Int maxSteps, BASE(F) tol, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
    if( maxSteps > std::min(A.Height(),A.Width()) )
        LogicError("Too many steps requested");
#endif
    typedef BASE(F) Real;
    p.ResizeTo( maxSteps, 1 );
    t.ResizeTo( maxSteps, 1 );

    Matrix<F>
        ATL, ATR,  A00, a01,     A02,  aLeftCol, ARightPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;
    Matrix<F> z;

    const Int m = A.Height();
    std::vector<F> swapBuf( m );

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the 
    // updated quantities are no longer accurate.
    std::vector<Real> origNorms;
    const Real maxOrigNorm = ColumnNorms( A, origNorms );
    std::vector<Real> norms = origNorms;
    const Real updateTol = Sqrt(lapack::MachineEpsilon<Real>());

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    for( Int col=0; col<maxSteps; ++col )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        View2x1( aLeftCol, alpha11,
                           a21 );

        View2x1( ARightPan, a12,
                            A22 );

        //--------------------------------------------------------------------//
        // Find the next column pivot
        const ValueInt<Real> pivot = FindPivot( norms, col );
        if( pivot.value <= tol*maxOrigNorm )
        {
            p.ResizeTo( col, 1 );
            t.ResizeTo( col, 1 );
            break;
        }
        p.Set( col, 0, pivot.index );
 
        // Perform the swap
        if( col != pivot.index )
        {
            MemSwap( A.Buffer(0,col), A.Buffer(0,pivot.index), &swapBuf[0], m );
            norms[pivot.index] = norms[col];
            origNorms[pivot.index] = origNorms[col];
        }

        // Compute and apply the Householder reflector for this column
        const F tau = Reflector( alpha11, a21 );
        t.Set( col, 0, tau );
        const F alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Zeros( z, ARightPan.Width(), 1 );
        Gemv( ADJOINT, F(1), ARightPan, aLeftCol, F(0), z );
        Ger( -Conj(tau), aLeftCol, z, ARightPan );
        alpha11.Set(0,0,alpha);

        // Update the column norm estimates in the same manner as LAWN 176
        for( Int k=0; k<a12.Width(); ++k )
        {
            const Int j = k + col+1;    
            if( norms[j] != Real(0) )
            {
                Real gamma = Abs(a12.Get(0,k)) / norms[j];
                gamma = std::max( Real(0), (Real(1)-gamma)*(Real(1)+gamma) );
                const Real ratio = norms[j] / origNorms[j];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol || alwaysRecompute )
                {
                    norms[j] = blas::Nrm2( m-(col+1), A.Buffer(col+1,j), 1 );
                    origNorms[j] = norms[j];
                }
                else
                    norms[j] *= Sqrt(gamma);
            }
        }
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
}

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<Int>& p,
  Int maxSteps, BASE(F) tol, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    Matrix<F> t;
    BusingerGolub( A, t, p, maxSteps, tol, alwaysRecompute );
}

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Int>& p,
  Int numSteps, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    // Use a tolerance of -1 so that we always perform numSteps iterations
    BusingerGolub( A, t, p, numSteps, BASE(F)(-1), alwaysRecompute );
}

// If we don't need 't' from the above routine
template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<Int>& p, Int numSteps, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    Matrix<F> t;
    BusingerGolub( A, t, p, numSteps, alwaysRecompute );
}

template<typename F> 
inline void
BusingerGolub
( Matrix<F>& A, Matrix<F>& t, Matrix<Int>& p, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    const Int numSteps = std::min(A.Height(),A.Width());
    BusingerGolub( A, t, p, numSteps, alwaysRecompute );
}

// If we don't need 't' from the above routine
template<typename F> 
inline void
BusingerGolub( Matrix<F>& A, Matrix<Int>& p, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    Matrix<F> t;
    BusingerGolub( A, t, p, alwaysRecompute );
}

template<typename F>
inline ValueInt<BASE(F)>
FindColumnPivot
( const DistMatrix<F>& A, const std::vector<BASE(F)>& norms, Int col )
{
#ifndef RELEASE
    CallStackEntry entry("qr::FindColumnPivot");
#endif
    typedef BASE(F) Real;
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
#ifndef RELEASE
    CallStackEntry entry("qr::ColumnNorms");
#endif
    typedef BASE(F) Real;
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
    ( &localScales[0], &scales[0], localWidth, mpi::MAX, colComm );

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
    ( &localScaledSquares[0], &scaledSquares[0], localWidth, colComm );

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
#ifndef RELEASE
    CallStackEntry entry("qr::ReplaceColumnNorms");
#endif
    typedef BASE(F) Real;
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
    ( &localScales[0], &scales[0], numInaccurate, mpi::MAX, colComm );

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
    ( &localScaledSquares[0], &scaledSquares[0], numInaccurate, colComm );

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
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
    if( maxSteps > std::min(A.Height(),A.Width()) )
        LogicError("Too many steps requested");
    if( A.Grid() != p.Grid() || A.Grid() != t.Grid() )
        LogicError("A, t, and p must have the same grid");
#endif
    typedef BASE(F) Real;
    const Grid& g = A.Grid();
    t.ResizeTo( maxSteps, 1 );
    p.ResizeTo( maxSteps, 1 );

    DistMatrix<F>
        ATL(g), ATR(g),  A00(g), a01(g),     A02(g),  aLeftCol(g), ARightPan(g),
        ABL(g), ABR(g),  a10(g), alpha11(g), a12(g),
                         A20(g), a21(g),     A22(g);
    DistMatrix<F> z(g);

    const Int mLocal = A.LocalHeight();
    const Int nLocal = A.LocalWidth();
    const Int rowAlign = A.RowAlignment();
    const Int rowShift = A.RowShift();
    const Int rowStride = A.RowStride();
    std::vector<F> swapBuf( mLocal );

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the 
    // updated quantities are no longer accurate.
    std::vector<Real> origNorms( nLocal );
    const Real maxOrigNorm = ColumnNorms( A, origNorms );
    std::vector<Real> norms = origNorms;
    const Real updateTol = Sqrt(lapack::MachineEpsilon<Real>());
    std::vector<Int> inaccurateNorms;

    // Temporary distributions
    DistMatrix<F,MC,STAR> aLeftCol_MC_STAR(g);
    DistMatrix<F,MR,STAR> z_MR_STAR(g);
    DistMatrix<F,STAR,MR> a12_STAR_MR(g);

    PartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    for( Int col=0; col<maxSteps; ++col )
    {
        RepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ a01,     A02,
         /*************/ /**********************/
               /**/       a10, /**/ alpha11, a12,
          ABL, /**/ ABR,  A20, /**/ a21,     A22, 1 );

        View2x1( aLeftCol, alpha11,
                           a21 );

        View2x1( ARightPan, a12,
                            A22 );

        aLeftCol_MC_STAR.AlignWith( ARightPan );
        z_MR_STAR.AlignWith( ARightPan );
        //--------------------------------------------------------------------//
        // Find the next column pivot
        const ValueInt<Real> pivot = FindColumnPivot( A, norms, col );
        if( pivot.value <= tol*maxOrigNorm )
        {
            p.ResizeTo( col, 1 );
            t.ResizeTo( col, 1 );
            break;
        }
        p.Set( col, 0, pivot.index );

        // Perform the swap
        const Int colOwner = (col+rowAlign) % rowStride;
        const Int pivotColOwner = (pivot.index+rowAlign) % rowStride;
        const bool myCol = ( g.Col() == colOwner );
        const bool myPivotCol = ( g.Col() == pivotColOwner );
        if( col != pivot.index )
        {
            if( myCol && myPivotCol )
            {
                const Int colLocal = (col-rowShift) / rowStride;
                const Int pivotColLocal = (pivot.index-rowShift) / rowStride;
                MemSwap
                ( A.Buffer(0,colLocal), A.Buffer(0,pivotColLocal),
                  &swapBuf[0], mLocal );
                norms[pivotColLocal] = norms[colLocal];
                origNorms[pivotColLocal] = origNorms[colLocal];
            }
            else if( myCol )
            {
                const Int colLocal = (col-rowShift) / rowStride;
                mpi::SendRecv
                ( A.Buffer(0,colLocal), mLocal, 
                  pivotColOwner, pivotColOwner, g.RowComm() );
                mpi::Send( norms[colLocal], pivotColOwner, g.RowComm() );
            }
            else if( myPivotCol )
            {
                const Int pivotColLocal = (pivot.index-rowShift) / rowStride;
                mpi::SendRecv
                ( A.Buffer(0,pivotColLocal), mLocal,
                  colOwner, colOwner, g.RowComm() );
                norms[pivotColLocal] = mpi::Recv<Real>( colOwner, g.RowComm() );
            }
        }

        // Compute the Householder reflector
        const F tau = Reflector( alpha11, a21 );
        t.Set( A00.Width(), 0, tau );

        // Apply the Householder reflector
        const bool myDiagonalEntry = ( g.Row() == alpha11.ColAlignment() &&
                                       g.Col() == alpha11.RowAlignment() );
        F alpha = 0;
        if( myDiagonalEntry )
        {
            alpha = alpha11.GetLocal(0,0);
            alpha11.SetLocal(0,0,1);
        }
        aLeftCol_MC_STAR = aLeftCol;
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
        for( Int kLocal=0; kLocal<a12LocalWidth; ++kLocal )
        {
            const Int k = a12RowShift + kLocal*rowStride;
            const Int j = k + col+1;
            const Int jLoc = (j-rowShift) / rowStride;
            if( norms[jLoc] != Real(0) )
            {
                const Real beta = Abs(a12_STAR_MR.GetLocal(0,kLocal));
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
        //--------------------------------------------------------------------//

        SlidePartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, a01,     /**/ A02,
               /**/       a10, alpha11, /**/ a12,
         /*************/ /**********************/
          ABL, /**/ ABR,  A20, a21,     /**/ A22 );
    }
}

// If we don't need 't' from the above routine
template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<Int,VR,STAR>& p,
  Int maxSteps, BASE(F) tol, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    DistMatrix<F,MD,STAR> t( A.Grid() );
    BusingerGolub( A, t, p, maxSteps, tol, alwaysRecompute );
}

template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Int,VR,STAR>& p,
  Int numSteps, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
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
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    DistMatrix<F,MD,STAR> t( A.Grid() );
    BusingerGolub( A, t, p, numSteps, alwaysRecompute );
}

template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<F,MD,STAR>& t, DistMatrix<Int,VR,STAR>& p,
  bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    const Int numSteps = std::min(A.Height(),A.Width());
    BusingerGolub( A, t, p, numSteps, alwaysRecompute );
}

// If we don't need 't' from the above routine
template<typename F>
inline void
BusingerGolub
( DistMatrix<F>& A, DistMatrix<Int,VR,STAR>& p, bool alwaysRecompute=false )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
#endif
    DistMatrix<F,MD,STAR> t( A.Grid() );
    BusingerGolub( A, t, p, alwaysRecompute );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_LAPACK_QR_BUSINGERGOLUB_HPP
