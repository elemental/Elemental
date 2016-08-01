/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LU_PANEL_HPP
#define EL_LU_PANEL_HPP

namespace El {
namespace lu {

template<typename F>
void Panel( Matrix<F>& A, Permutation& P, Permutation& PB, Int offset )
{
    DEBUG_CSE
    const Int m = A.Height(); 
    const Int n = A.Width();
    F* ABuf = A.Buffer();
    const Int ALDim = A.LDim();
    DEBUG_ONLY(
      if( A.Height() < n )
          LogicError("Must be a column panel");
    )

    PB.MakeIdentity( A.Height() );
    PB.ReserveSwaps( n );

    for( Int k=0; k<n; ++k )
    {
        const Int ind2HorzSize = n - (k+1);
        const Int ind2VertSize = m - (k+1);
        const F* a12Buf = &ABuf[ k    + (k+1)*ALDim];
        const F* aB1Buf = &ABuf[ k    +  k   *ALDim];
              F* a21Buf = &ABuf[(k+1) +  k   *ALDim];
              F* A22Buf = &ABuf[(k+1) + (k+1)*ALDim];

        // Find the index and value of the pivot candidate
        const Int maxInd = blas::MaxInd( ind2VertSize+1, aB1Buf, 1 );
        const Int iPiv = maxInd + k;
        P.RowSwap( k+offset, iPiv+offset );
        PB.RowSwap( k, iPiv );

        // Swap the pivot row and current row
        if( iPiv != k )
            blas::Swap( n, &ABuf[k], ALDim, &ABuf[iPiv], ALDim );

        // Now we can perform the update of the current panel
        const F alpha = ABuf[k+k*ALDim];
        if( alpha == F(0) )
            throw SingularMatrixException();
        const F alpha11Inv = F(1) / alpha;
        blas::Scal( ind2VertSize, alpha11Inv, a21Buf, 1 );
        blas::Geru
        ( ind2VertSize, ind2HorzSize,
          F(-1), a21Buf, 1, a12Buf, ALDim, A22Buf, ALDim );
    }
}

// NOTE: It is assumed that the local buffers of A[*,*] and B[MC,*] can be
//       verticially stacked, so that the top-left local entry of B is 
//       the n'th local entry of A[*,*]'s local buffer.
//       Also, on entry, it is only required that process row 0 has the correct
//       data for A.
template<typename F>
void Panel
( DistMatrix<F,  STAR,STAR>& A, 
  DistMatrix<F,  MC,  STAR>& B, 
  DistPermutation& P,
  DistPermutation& PB,
  Int offset,
  vector<F>& pivotBuffer )
{
    DEBUG_CSE
    typedef Base<F> Real;
    const Int n = A.Width();
    const Int BLocHeight = B.LocalHeight();
    F* ABuf = A.Buffer();
    F* BBuf = B.Buffer();
    const Int ALDim = A.LDim();
    const Int BLDim = B.LDim();
    mpi::Comm colComm = B.ColComm();
    mpi::Op maxLocOp = mpi::MaxLocOp<Real>();
    DEBUG_ONLY(
      AssertSameGrids( A, B );
      if( n != B.Width() )
          LogicError("A and B must be the same width");
      if( A.Buffer()+n != B.Buffer() )
          LogicError("Buffers of A and B did not properly align");
    )

    PB.MakeIdentity( A.Height()+B.Height() );
    PB.ReserveSwaps( n );

    pivotBuffer.resize( n );
    for( Int k=0; k<n; ++k )
    {
        const Int ind2Size = n-k-1;
        const F* a12Buf = &ABuf[ k    + (k+1)*ALDim];
        const F* aB1Buf = &ABuf[ k    +  k   *ALDim];
              F* a21Buf = &ABuf[(k+1) +  k   *ALDim];
              F* A22Buf = &ABuf[(k+1) + (k+1)*ALDim];

        // Store the index/value of the local pivot candidate
        Int aB1LocalInd = blas::MaxInd( ind2Size+1+BLocHeight, aB1Buf, 1 ); 
        ValueInt<Real> localPivot;
        localPivot.value = Abs(aB1Buf[aB1LocalInd]);
        if( aB1LocalInd+k < n )
            localPivot.index = aB1LocalInd + k;
        else
            localPivot.index = B.GlobalRow(aB1LocalInd-(n-k)) + n;

        // Compute and store the location of the new pivot
        const auto pivot = mpi::AllReduce( localPivot, maxLocOp, colComm );
        const Int iPiv = pivot.index;
        P.RowSwap( k+offset, iPiv+offset );
        PB.RowSwap( k, iPiv );

        // Perform the pivot within this panel
        if( iPiv < n )
        {
            blas::Swap( n, &ABuf[iPiv], ALDim, &ABuf[k], ALDim );
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
            mpi::Broadcast( pivotBuffer.data(), n, ownerRow, colComm );
            // Overwrite the current row with the pivot row
            for( Int j=0; j<n; ++j )
                ABuf[k+j*ALDim] = pivotBuffer[j];
        }

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
