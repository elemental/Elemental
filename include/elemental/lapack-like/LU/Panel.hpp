/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_LU_PANEL_HPP
#define ELEM_LAPACK_LU_PANEL_HPP

#include "elemental/blas-like/level1/Scale.hpp"
#include "elemental/blas-like/level1/Swap.hpp"
#include "elemental/blas-like/level2/Geru.hpp"

namespace elem {
namespace lu {

template<typename F>
inline void
Panel( Matrix<F>& A, Matrix<Int>& p, Int pivotOffset=0 )
{
#ifndef RELEASE
    CallStackEntry entry("lu::Panel");
#endif
    const Int m = A.Height();
    const Int n = A.Width();
#ifndef RELEASE
    if( m < n )
        LogicError("Must be a column panel");
#endif
    p.ResizeTo( n, 1 );

    for( Int k=0; k<n; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );

        // Find the index and value of the pivot candidate
        auto pivot = VectorMax( ViewRange(A,k,k,m,k+1) );
        const Int iPiv = pivot.index + k;
        p.Set( k, 0, iPiv+pivotOffset );

        // Swap the pivot row and current row
        if( iPiv != k )
        {
            auto aCurRow = ViewRange( A, k,    0, k+1,    n );
            auto aPivRow = ViewRange( A, iPiv, 0, iPiv+1, n );
            Swap( NORMAL, aCurRow, aPivRow );
        }

        // Now we can perform the update of the current panel
        const F alpha = alpha11.Get(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        Scale( alpha11Inv, a21 );
        Geru( F(-1), a21, a12, A22 );
    }
}

template<typename F>
inline void
Panel
( DistMatrix<F,  STAR,STAR>& A, 
  DistMatrix<F,  MC,  STAR>& B, 
  DistMatrix<Int,STAR,STAR>& p, 
  Int pivotOffset=0 )
{
#ifndef RELEASE
    CallStackEntry entry("lu::Panel");
    if( A.Grid() != p.Grid() || p.Grid() != B.Grid() )
        LogicError("Matrices must be distributed over the same grid");
    if( A.Width() != B.Width() )
        LogicError("A and B must be the same width");
    if( A.Height() != p.Height() || p.Width() != 1 )
        LogicError("p must be a vector that conforms with A");
#endif
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    const Int r = g.Height();
    const Int colShift = B.ColShift();
    const Int colAlign = B.ColAlign();

    // For packing rows of data for pivoting
    const Int n = A.Width();
    const Int mB = B.Height();
    const Int nB = B.Width();
    std::vector<F> pivotBuffer( n );

    for( Int k=0; k<n; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   n,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, n,   n   );
        auto b1      = ViewRange( B, 0,   k,   mB,  k+1 );
        auto B2      = ViewRange( B, 0,   k+1, mB,  nB  );

        // Store the index/value of the local pivot candidate
        ValueInt<Real> localPivot;
        localPivot.value = FastAbs(alpha11.GetLocal(0,0));
        localPivot.index = k;
        for( Int i=0; i<a21.Height(); ++i )
        {
            const Real value = FastAbs(a21.GetLocal(i,0));
            if( value > localPivot.value )
            {
                localPivot.value = value;
                localPivot.index = k + i + 1;
            }
        }
        for( Int i=0; i<B.LocalHeight(); ++i )
        {
            const Real value = FastAbs(b1.GetLocal(i,0));
            if( value > localPivot.value )
            {
                localPivot.value = value;
                localPivot.index = n + colShift + i*r;
            }
        }

        // Compute and store the location of the new pivot
        const ValueInt<Real> pivot = 
            mpi::AllReduce( localPivot, mpi::MaxLocOp<Real>(), g.ColComm() );
        const Int iPiv = pivot.index;
        p.SetLocal( k, 0, iPiv+pivotOffset );

        // Perform the pivot within this panel
        if( iPiv < n )
        {
            // Pack pivot into temporary
            for( Int j=0; j<n; ++j )
                pivotBuffer[j] = A.GetLocal( iPiv, j );
            // Replace pivot with current
            for( Int j=0; j<n; ++j )
                A.SetLocal( iPiv, j, A.GetLocal(k,j) );
        }
        else
        {
            // The owning row of the pivot row packs it into the row buffer
            // and then overwrites with the current row
            const Int relIndex = iPiv - n;
            const Int ownerRow = (colAlign+relIndex) % r;
            if( g.Row() == ownerRow )
            {
                const int iLoc = (relIndex-colShift) / r;
                for( Int j=0; j<n; ++j )
                    pivotBuffer[j] = B.GetLocal( iLoc, j );
                for( Int j=0; j<n; ++j )
                    B.SetLocal( iLoc, j, A.GetLocal(k,j) );
            }
            // The owning row broadcasts within process columns
            mpi::Broadcast( pivotBuffer.data(), n, ownerRow, g.ColComm() );
        }
        // Overwrite the current row with the pivot row
        for( Int j=0; j<n; ++j )
            A.SetLocal( k, j, pivotBuffer[j] );

        // Now we can perform the update of the current panel
        const F alpha = alpha11.GetLocal(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        Scale( alpha11Inv, a21 );
        Scale( alpha11Inv, b1  );
        Geru( F(-1), a21.Matrix(), a12.Matrix(), A22.Matrix() );
        Geru( F(-1), b1.Matrix(), a12.Matrix(), B2.Matrix() );
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LAPACK_LU_PANEL_HPP
