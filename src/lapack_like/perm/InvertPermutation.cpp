/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

void InvertPermutation( const Matrix<Int>& p, Matrix<Int>& pInv )
{
    DEBUG_ONLY(
      CallStackEntry cse("InvertPermutation");
      if( p.Width() != 1 )
          LogicError("p must be a column vector");
    )
    const Int n = p.Height();
    pInv.Resize( n, 1 );
    if( n == 0 )
        return;

    DEBUG_ONLY(
      // This is obviously necessary but not sufficient for 'p' to contain
      // a reordering of (0,1,...,n-1).
      const Int range = MaxNorm( p ) + 1;
      if( range != n )
          LogicError("Invalid putation range");
    )

    for( Int i=0; i<n; ++i )
        pInv.Set( p.Get(i,0), 0, i );
}

void InvertPermutation
( const AbstractDistMatrix<Int>& pPre, AbstractDistMatrix<Int>& pInvPre )
{
    DEBUG_ONLY(
      CallStackEntry cse("InvertPermutation");
      if( pPre.Width() != 1 )
          LogicError("p must be a column vector");
    )

    const Int n = pPre.Height();
    pInvPre.AlignWith( pPre, false );
    pInvPre.Resize( n, 1 );
    if( n == 0 )
        return;

    auto pPtr = ReadProxy<Int,VC,STAR>( &pPre ); 
    auto& p = *pPtr;

    auto pInvPtr = WriteProxy<Int,VC,STAR>( &pInvPre ); 
    auto& pInv = *pInvPtr;

    DEBUG_ONLY(
      // This is obviously necessary but not sufficient for 'p' to contain
      // a reordering of (0,1,...,n-1).
      const Int range = MaxNorm( p ) + 1;
      if( range != n )
          LogicError("Invalid putation range");
    )

    // Compute the send counts
    const mpi::Comm colComm = p.ColComm();
    const Int commSize = mpi::Size( colComm );
    vector<int> sendSizes(commSize,0), recvSizes(commSize,0);
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
    {
        const Int iDest = p.GetLocal(iLoc,0);
        const int owner = pInv.RowOwner(iDest);
        sendSizes[owner] += 2; // we'll send the global index and the value
    }
    // Perform a small AllToAll to get the receive counts
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, colComm );
    vector<int> sendOffs, recvOffs;
    const int sendTotal = Scan( sendSizes, sendOffs );
    const int recvTotal = Scan( recvSizes, recvOffs );

    // Pack the send data
    vector<Int> sendBuf(sendTotal);
    auto offsets = sendOffs;
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
    {
        const Int i     = p.GlobalRow(iLoc);
        const Int iDest = p.GetLocal(iLoc,0);
        const int owner = pInv.RowOwner(iDest);
        sendBuf[offsets[owner]++] = iDest;
        sendBuf[offsets[owner]++] = i;
    }

    // Perform the actual exchange
    vector<Int> recvBuf(recvTotal);
    mpi::AllToAll
    ( sendBuf.data(), sendSizes.data(), sendOffs.data(),
      recvBuf.data(), recvSizes.data(), recvOffs.data(), colComm );
    SwapClear( sendBuf );
    SwapClear( sendSizes );
    SwapClear( sendOffs );

    // Unpack the received data
    for( Int k=0; k<recvTotal/2; ++k )
    {
        const Int iDest = recvBuf[2*k+0];
        const Int i     = recvBuf[2*k+1];

        const Int iDestLoc = pInv.LocalRow(iDest);
        pInv.SetLocal( iDestLoc, 0, i );
    }
}

} // namespace El
