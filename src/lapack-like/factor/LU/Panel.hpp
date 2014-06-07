/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LU_PANEL_HPP
#define EL_LU_PANEL_HPP

namespace El {
namespace lu {

template<typename F>
void Panel( Matrix<F>& A, Matrix<Int>& pivots )
{
    DEBUG_ONLY(CallStackEntry cse("lu::Panel"))
    const Int m = A.Height();
    const Int n = A.Width();
    DEBUG_ONLY(
        if( m < n )
            LogicError("Must be a column panel");
    )
    pivots.Resize( n, 1 );

    for( Int k=0; k<n; ++k )
    {
        auto alpha11 = ViewRange( A, k,   k,   k+1, k+1 );
        auto a12     = ViewRange( A, k,   k+1, k+1, n   );
        auto a21     = ViewRange( A, k+1, k,   m,   k+1 );
        auto A22     = ViewRange( A, k+1, k+1, m,   n   );

        // Find the index and value of the pivot candidate
        auto pivot = VectorMaxAbs( ViewRange(A,k,k,m,k+1) );
        const Int iPiv = pivot.index + k;
        pivots.Set( k, 0, iPiv );

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
void Panel
( DistMatrix<F,  STAR,STAR>& A, 
  DistMatrix<F,  MC,  STAR>& B, 
  DistMatrix<Int,STAR,STAR>& pivots )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::Panel");
        if( A.Grid() != pivots.Grid() || pivots.Grid() != B.Grid() )
            LogicError("Matrices must be distributed over the same grid");
        if( A.Width() != B.Width() )
            LogicError("A and B must be the same width");
    )
    typedef Base<F> Real;

    // For packing rows of data for pivoting
    const Int n = A.Width();
    const Int mB = B.Height();
    const Int nB = B.Width();
    std::vector<F> pivotBuffer( n );

    pivots.Resize( n, 1 );

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
        for( Int iLoc=0; iLoc<B.LocalHeight(); ++iLoc )
        {
            const Real value = FastAbs(b1.GetLocal(iLoc,0));
            if( value > localPivot.value )
            {
                localPivot.value = value;
                localPivot.index = n + B.GlobalRow(iLoc);
            }
        }

        // Compute and store the location of the new pivot
        const ValueInt<Real> pivot = 
            mpi::AllReduce( localPivot, mpi::MaxLocOp<Real>(), B.ColComm() );
        const Int iPiv = pivot.index;
        pivots.SetLocal( k, 0, iPiv );

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
            const Int ownerRow = B.RowOwner(relIndex);
            if( B.IsLocalRow(relIndex) )
            {
                const Int iLoc = B.LocalRow(relIndex);
                for( Int j=0; j<n; ++j )
                    pivotBuffer[j] = B.GetLocal( iLoc, j );
                for( Int j=0; j<n; ++j )
                    B.SetLocal( iLoc, j, A.GetLocal(k,j) );
            }
            // The owning row broadcasts within process columns
            mpi::Broadcast( pivotBuffer.data(), n, ownerRow, B.ColComm() );
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
} // namespace El

#endif // ifndef EL_LU_PANEL_HPP
