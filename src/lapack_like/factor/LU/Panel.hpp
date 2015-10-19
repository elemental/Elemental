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
    F* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    DEBUG_ONLY(
      CSE cse("lu::Panel");
      if( A.Height() < n )
          LogicError("Must be a column panel");
    )
    pivots.Resize( n, 1 );

    for( Int k=0; k<n; ++k )
    {
        const Int ind2Size = n - (k+1);
        const F* a12Buf = &ABuf[ k    + (k+1)*ALDim];
        const F* aB1Buf = &ABuf[ k    +  k   *ALDim];
              F* a21Buf = &ABuf[(k+1) +  k   *ALDim];
              F* A22Buf = &ABuf[(k+1) + (k+1)*ALDim];

        // Find the index and value of the pivot candidate
        const Int maxInd = blas::MaxInd( ind2Size+1, aB1Buf, 1 );
        const Int maxAbs = Abs(aB1Buf[maxInd]);
        const Int iPiv = maxInd + k;
        pivots.Set( k, 0, iPiv );

        // Swap the pivot row and current row
        if( iPiv != k )
            blas::Swap( n, &ABuf[k], ALDim, &ABuf[iPiv], ALDim );

        // Now we can perform the update of the current panel
        const F alpha = ABuf[k+k*ALDim];
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        blas::Scal( ind2Size, alpha11Inv, a21Buf, 1 );
        blas::Geru
        ( ind2Size, ind2Size, F(-1), a21Buf, 1, a12Buf, ALDim, A22Buf, ALDim );
    }
}

// NOTE: It is assumed that the local buffers of A[*,*] and B[MC,*] can be
//       verticially stacked, so that the top-left local entry of B is 
//       the n'th local entry of A[*,*]'s local buffer.
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
    mpi::Comm BColComm = B.ColComm();
    mpi::Op maxLocOp = mpi::MaxLocOp<Real>();
    DEBUG_ONLY(
      CSE cse("lu::Panel");
      AssertSameGrids( A, B, pivots );
      if( n != B.Width() )
          LogicError("A and B must be the same width");
      if( A.Buffer()+n != B.Buffer() )
          LogicError("Buffers of A and B did not properly align");
    )

    pivots.Resize( n, 1 );
    pivotBuffer.resize( n );
    for( Int k=0; k<n; ++k )
    {
        const Int ind2Size = n-k-1;
        const F* a12Buf = &ABuf[ k    + (k+1)*ALDim];
        const F* aB1Buf = &ABuf[ k    +  k   *ALDim];
              F* a21Buf = &ABuf[(k+1) +  k   *ALDim];
              F* A22Buf = &ABuf[(k+1) + (k+1)*ALDim];

        // Store the index/value of the local pivot candidate
        const Int aB1LocalInd =
          blas::MaxInd( ind2Size+1+BLocHeight, aB1Buf, 1 ); 
        const Real aB1LocalVal = Abs(aB1Buf[aB1LocalInd]);
        ValueInt<Real> localPivot;
        localPivot.value = aB1LocalVal;
        if( aB1LocalInd+k < n )
        {
            localPivot.index = aB1LocalInd + k;
        }
        else
        {
            const Int b1LocalInd = aB1LocalInd-(n-k);
            localPivot.index = B.GlobalRow(b1LocalInd) + n;
        }

        // Compute and store the location of the new pivot
        const auto pivot = mpi::AllReduce( localPivot, maxLocOp, BColComm );
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
            mpi::Broadcast( pivotBuffer.data(), n, ownerRow, BColComm );
        }
        // Overwrite the current row with the pivot row
        for( Int j=0; j<n; ++j )
            ABuf[k+j*ALDim] = pivotBuffer[j];

        // Now we can perform the update of the current panel
        const F alpha = aB1Buf[0];
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        blas::Scal( ind2Size+BLocHeight, alpha11Inv, a21Buf, 1 );
        blas::Geru
        ( ind2Size+BLocHeight, ind2Size, F(-1),
          a21Buf, 1, a12Buf, ALDim, A22Buf, ALDim );
    }
}

} // namespace lu
} // namespace El

#endif // ifndef EL_LU_PANEL_HPP
