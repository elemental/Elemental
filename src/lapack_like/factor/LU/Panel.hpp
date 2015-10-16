/*
   Copyright (c) 2009-2015, Jack Poulson
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
    const Int n = A.Width();
    DEBUG_ONLY(
      CSE cse("lu::Panel");
      if( A.Height() < n )
          LogicError("Must be a column panel");
    )
    pivots.Resize( n, 1 );

    for( Int k=0; k<n; ++k )
    {
        const Range<Int> ind1( k ), ind2( k+1, END );

        auto alpha11 = A( ind1, ind1 );
        auto a12     = A( ind1, ind2 );
        auto a21     = A( ind2, ind1 );
        auto A22     = A( ind2, ind2 );

        // Find the index and value of the pivot candidate
        auto pivot = VectorMaxAbs( A(IR(k,END),IR(k)) );
        const Int iPiv = pivot.index + k;
        pivots.Set( k, 0, iPiv );

        // Swap the pivot row and current row
        if( iPiv != k )
        {
            auto aCurRow = A( ind1,     ALL );
            auto aPivRow = A( IR(iPiv), ALL );
            Swap( NORMAL, aCurRow, aPivRow );
        }

        // Now we can perform the update of the current panel
        const F alpha = alpha11.Get(0,0);
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        a21 *= alpha11Inv;
        Geru( F(-1), a21, a12, A22 );
    }
}

template<typename F>
void Panel
( DistMatrix<F,  STAR,STAR>& A, 
  DistMatrix<F,  MC,  STAR>& B, 
  DistMatrix<Int,STAR,STAR>& pivots,
  vector<F>& pivotBuffer )
{
    typedef Base<F> Real;
    const Int n = A.Width();
    const Int BLocHeight = B.LocalHeight();
    F* ABuf = A.Buffer();
    F* BBuf = B.Buffer();
    const Int ALDim = A.LDim();
    const Int BLDim = B.LDim();
    DEBUG_ONLY(
      CSE cse("lu::Panel");
      AssertSameGrids( A, B, pivots );
      if( n != B.Width() )
          LogicError("A and B must be the same width");
    )

    pivots.Resize( n, 1 );
    pivotBuffer.resize( n );
    for( Int k=0; k<n; ++k )
    {
        const Range<Int> ind1( k ), ind2( k+1, END );
        const Int a21Height = n-k-1;

        auto a12 = A( ind1, ind2 );
        auto a21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto b1 = B( ALL, ind1 );
        auto B2 = B( ALL, ind2 );

        // Store the index/value of the local pivot candidate
        ValueInt<Real> localPivot;
        localPivot.value = FastAbs(ABuf[k+k*ALDim]);
        localPivot.index = k;
        for( Int i=0; i<a21Height; ++i )
        {
            const Real value = FastAbs(ABuf[(i+k+1)+k*ALDim]);
            if( value > localPivot.value )
            {
                localPivot.value = value;
                localPivot.index = k + i + 1;
            }
        }
        for( Int iLoc=0; iLoc<BLocHeight; ++iLoc )
        {
            const Real value = FastAbs(BBuf[iLoc+k*BLDim]);
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
                pivotBuffer[j] = ABuf[iPiv+j*ALDim];
            // Replace pivot with current
            for( Int j=0; j<n; ++j )
                ABuf[iPiv+j*ALDim] = ABuf[k+j*ALDim];
        }
        else
        {
            // The owning row of the pivot row packs it into the row buffer
            // and then overwrites with the current row
            const Int relIndex = iPiv - n;
            const int ownerRow = B.RowOwner(relIndex);
            if( B.IsLocalRow(relIndex) )
            {
                const Int iLoc = B.LocalRow(relIndex);
                for( Int j=0; j<n; ++j )
                    pivotBuffer[j] = BBuf[iLoc+j*BLDim];
                for( Int j=0; j<n; ++j )
                    BBuf[iLoc+j*BLDim] = ABuf[k+j*ALDim];
            }
            // The owning row broadcasts within process columns
            mpi::Broadcast( pivotBuffer.data(), n, ownerRow, B.ColComm() );
        }
        // Overwrite the current row with the pivot row
        for( Int j=0; j<n; ++j )
            ABuf[k+j*ALDim] = pivotBuffer[j];

        // Now we can perform the update of the current panel
        const F alpha = ABuf[k+k*ALDim];
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        a21 *= alpha11Inv;
        b1 *= alpha11Inv;
        Geru( F(-1), a21.Matrix(), a12.Matrix(), A22.Matrix() );
        Geru( F(-1), b1.Matrix(), a12.Matrix(), B2.Matrix() );
    }
}

} // namespace lu
} // namespace El

#endif // ifndef EL_LU_PANEL_HPP
