/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_QR_BUSINGERGOLUB_HPP
#define LAPACK_QR_BUSINGERGOLUB_HPP

#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/blas-like/level2/Ger.hpp"
#include "elemental/lapack-like/Reflector.hpp"
#include "elemental/matrices/Zeros.hpp"

#include <algorithm>

namespace elem {
namespace qr {

template<typename F>
inline void
ColumnNorms( const Matrix<F>& A, std::vector<BASE(F)>& norms )
{
#ifndef RELEASE
    CallStackEntry entry("qr::ColumnNorms");
#endif
    const int m = A.Height();
    const int n = A.Width();
    norms.resize( n );
    for( int j=0; j<n; ++j )
        norms[j] = blas::Nrm2( m, A.LockedBuffer(0,j), 1 );
}

template<typename Real>
inline int
FindPivot( const std::vector<Real>& norms, int col )
{
#ifndef RELEASE
    CallStackEntry entry("qr::FindPivot");
#endif
    const int n = norms.size();
    const Real* maxNorm = std::max_element( &norms[col], &norms[n] );
    return maxNorm - &norms[0];
}

template<typename Real>
inline void
BusingerGolub( Matrix<Real>& A, Matrix<int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
    if( p.Viewing() && 
        (p.Height() != std::min(A.Height(),A.Width()) || p.Width() != 1) )
        throw std::logic_error
        ("p must be a vector of the same height as the min dimension of A");
#endif
    if( !p.Viewing() )
        p.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    Matrix<Real>
        ATL, ATR,  A00, a01,     A02,  aLeftCol, ARightPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;

    Matrix<Real> z;

    const int m = A.Height();
    const int n = A.Width();
    std::vector<Real> swapBuf( m );

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the 
    // updated quantities are no longer accurate.
    std::vector<Real> origNorms;
    ColumnNorms( A, origNorms );
    std::vector<Real> norms = origNorms;
    const Real updateTol = Sqrt(lapack::MachineEpsilon<Real>());

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < m && ATL.Width() < n )
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
        const int col = A00.Width();
        const int pivotCol = FindPivot( norms, col );
        p.Set( col, 0, pivotCol );

        // Perform the swap
        if( col != pivotCol )
        {
            MemSwap( A.Buffer(0,col), A.Buffer(0,pivotCol), &swapBuf[0], m );
            norms[pivotCol] = norms[col];
            origNorms[pivotCol] = origNorms[col];
        }

        // Compute and apply the Householder reflector for this column
        const Real tau = Reflector( alpha11, a21 );
        const Real alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Zeros( z, ARightPan.Width(), 1 );
        Gemv( TRANSPOSE, Real(1), ARightPan, aLeftCol, Real(0), z );
        Ger( -tau, aLeftCol, z, ARightPan );
        alpha11.Set(0,0,alpha);

        // Update the column norm estimates in the same manner as LAWN 176
        for( int k=0; k<a12.Width(); ++k )
        {
            const int j = k + col+1;    
            if( norms[j] != Real(0) )
            {
                Real gamma = Abs(a12.Get(0,k)) / norms[j];
                gamma = std::max( Real(0), (Real(1)-gamma)*(Real(1)+gamma) );
                const Real ratio = norms[j] / origNorms[j];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol )
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

template<typename Real> 
inline void
BusingerGolub
( Matrix<Complex<Real> >& A,
  Matrix<Complex<Real> >& t,
  Matrix<int>& p )
{
#ifndef RELEASE
    CallStackEntry entry("qr::BusingerGolub");
    if( t.Viewing() && 
        (t.Height() != std::min(A.Height(),A.Width()) || t.Width() != 1) )
        throw std::logic_error
        ("t must be a vector of the same height as the min dimension of A");
    if( p.Viewing() && 
        (p.Height() != std::min(A.Height(),A.Width()) || p.Width() != 1) )
        throw std::logic_error
        ("p must be a vector of the same height as the min dimension of A");
#endif
    if( !p.Viewing() )
        p.ResizeTo( std::min(A.Height(),A.Width()), 1 );
    if( !t.Viewing() )
        t.ResizeTo( std::min(A.Height(),A.Width()), 1 );

    typedef Complex<Real> C;

    Matrix<C>
        ATL, ATR,  A00, a01,     A02,  aLeftCol, ARightPan,
        ABL, ABR,  a10, alpha11, a12,
                   A20, a21,     A22;
    Matrix<C> z;

    const int m = A.Height();
    const int n = A.Width();
    std::vector<C> swapBuf( m );

    // Initialize two copies of the column norms, one will be consistently
    // updated, but the original copy will be kept to determine when the 
    // updated quantities are no longer accurate.
    std::vector<Real> origNorms;
    ColumnNorms( A, origNorms );
    std::vector<Real> norms = origNorms;
    const Real updateTol = Sqrt(lapack::MachineEpsilon<Real>());

    PartitionDownLeftDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    while( ATL.Height() < A.Height() && ATL.Width() < A.Width() )
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
        const int col = A00.Width();
        const int pivotCol = FindPivot( norms, col );
        p.Set( col, 0, pivotCol );
 
        // Perform the swap
        if( col != pivotCol )
        {
            MemSwap( A.Buffer(0,col), A.Buffer(0,pivotCol), &swapBuf[0], m );
            norms[pivotCol] = norms[col];
            origNorms[pivotCol] = origNorms[col];
        }

        // Compute and apply the Householder reflector for this column
        const C tau = Reflector( alpha11, a21 );
        t.Set( col, 0, tau );
        const C alpha = alpha11.Get(0,0);
        alpha11.Set(0,0,1);
        Zeros( z, ARightPan.Width(), 1 );
        Gemv( ADJOINT, C(1), ARightPan, aLeftCol, C(0), z );
        Ger( -Conj(tau), aLeftCol, z, ARightPan );
        alpha11.Set(0,0,alpha);

        // Update the column norm estimates in the same manner as LAWN 176
        for( int k=0; k<a12.Width(); ++k )
        {
            const int j = k + col+1;    
            if( norms[j] != Real(0) )
            {
                Real gamma = Abs(a12.Get(0,k)) / norms[j];
                gamma = std::max( Real(0), (Real(1)-gamma)*(Real(1)+gamma) );
                const Real ratio = norms[j] / origNorms[j];
                const Real phi = gamma*(ratio*ratio);
                if( phi <= updateTol )
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
ColumnNorms( const DistMatrix<F>& A, std::vector<BASE(F)>& norms )
{
#ifndef RELEASE
    CallStackEntry entry("qr::ColumnNorms");
#endif
    throw std::logic_error("Not yet written");
}

template<typename F>
inline int
FindColumnPivot
( const DistMatrix<F>& A, const std::vector<BASE(F)>& norms, int col )
{
#ifndef RELEASE
    CallStackEntry entry("qr::FindColumnPivot");
#endif
    typedef BASE(F) Real;
    const int rowShift = A.RowShift();
    const int rowStride = A.RowStride();
    const int localColsBefore = Length( col, rowShift, rowStride );
    const int localPivot = FindPivot( norms, localColsBefore );
    mpi::ValueInt<Real> pivotInfo;
    pivotInfo.value = norms[localPivot];
    pivotInfo.index = rowShift+localPivot*rowStride;
    mpi::AllReduce( &pivotInfo, 1, mpi::MAXLOC, A.Grid().RowComm() );
}

// TODO: Parallel versions

} // namespace qr
} // namespace elem

#endif // ifndef LAPACK_QR_BUSINGERGOLUB_HPP
