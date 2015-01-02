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

    const mpi::Comm colComm = p.ColComm();
    const Int commSize = mpi::Size( colComm );
    std::vector<int> sendCounts(commSize,0), sendDispls(commSize),
                     recvCounts(commSize,0), recvDispls(commSize);

    // Compute the send counts
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
    {
        const Int iDest = p.GetLocal(iLoc,0);
        const Int owner = pInv.RowOwner(iDest);
        sendCounts[owner] += 2; // we'll send the global index and the value
    }
    // Perform a small AllToAll to get the receive counts
    mpi::AllToAll( sendCounts.data(), 1, recvCounts.data(), 1, colComm );

    // Compute the displacements
    int sendTotal=0, recvTotal=0;
    for( Int q=0; q<commSize; ++q )
    {
        sendDispls[q] = sendTotal;
        recvDispls[q] = recvTotal;
        sendTotal += sendCounts[q];
        recvTotal += recvCounts[q];
    }

    // Pack the send data
    std::vector<Int> sendBuf(sendTotal);
    auto offsets = sendDispls;
    for( Int iLoc=0; iLoc<p.LocalHeight(); ++iLoc )
    {
        const Int i     = p.GlobalRow(iLoc);
        const Int iDest = p.GetLocal(iLoc,0);
        const Int owner = pInv.RowOwner(iDest);
        sendBuf[offsets[owner]++] = iDest;
        sendBuf[offsets[owner]++] = i;
    }

    // Perform the actual exchange
    std::vector<Int> recvBuf(recvTotal);
    mpi::AllToAll
    ( sendBuf.data(), sendCounts.data(), sendDispls.data(),
      recvBuf.data(), recvCounts.data(), recvDispls.data(), colComm );
    SwapClear( sendBuf );
    SwapClear( sendCounts );
    SwapClear( sendDispls );

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
