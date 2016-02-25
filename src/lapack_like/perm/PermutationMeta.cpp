/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

PermutationMeta::PermutationMeta
( const DistMatrix<Int,STAR,STAR>& perm,
  const DistMatrix<Int,STAR,STAR>& invPerm,
        Int permAlign,
        mpi::Comm permComm )
{
    DEBUG_ONLY(
      CSE cse("PermutationMeta::PermutationMeta");
      AssertSameGrids( perm, invPerm );
    )
    comm = permComm;
    align = permAlign;
    const Int permStride = mpi::Size( permComm );
    const Int permShift = Shift( mpi::Rank(permComm), permAlign, permStride );
    const Int b = perm.Height();
    const Int* permBuf = perm.LockedBuffer();
    const Int* invPermBuf = invPerm.LockedBuffer();

    // Form the metadata
    //
    // Extract the send and recv counts from the image and preimage.
    // There are three different types of exchanges:
    //   1. [0,b) -> [0,b)
    //   2. [0,b) -> [b,n)
    //   3. [b,n) -> [0,b)
    // The fourth possibility, [b,n) -> [b,n), is impossible due to the 
    // fact that indices pulled in from [b,n) are stuck in [0,b) due to the
    // fact that the i'th pivot exchanges index i with some index k >= i.
    // 
    sendCounts.resize( permStride, 0 );
    recvCounts.resize( permStride, 0 );
    FastResize( sendIdx, b );
    FastResize( recvIdx, b );
    FastResize( sendRanks, b );
    FastResize( recvRanks, b );

    sendIdx.resize( 0 );
    recvIdx.resize( 0 );
    sendRanks.resize( 0 );
    recvRanks.resize( 0 );
    for( Int i=0; i<b; ++i )
    {
        const Int preVal = permBuf[i];
        const Int postVal = invPermBuf[i];        

        // Handle sends 
        if( Mod(i,permStride) == permShift )
        {
            const Int iLoc = (i-permShift) / permStride;
            const Int sendTo = Mod(postVal+permAlign,permStride);
            sendIdx.push_back( iLoc );
            sendRanks.push_back( sendTo );
            ++sendCounts[sendTo];
        }
        if( preVal >= b && Mod(preVal,permStride) == permShift )
        {
            const Int iLoc = (preVal-permShift) / permStride;
            const Int sendTo = Mod(i+permAlign,permStride);
            sendIdx.push_back( iLoc );
            sendRanks.push_back( sendTo );
            ++sendCounts[sendTo];
        }

        // Handle recvs
        if( Mod(postVal,permStride) == permShift )
        {
            const Int iLoc = (postVal-permShift) / permStride;
            const Int recvFrom = Mod(i+permAlign,permStride);
            recvIdx.push_back( iLoc );
            recvRanks.push_back( recvFrom );
            ++recvCounts[recvFrom];
        }
        if( preVal >= b && Mod(i,permStride) == permShift )
        {
            const Int iLoc = (i-permShift) / permStride;
            const Int recvFrom = Mod(preVal+permAlign,permStride);
            recvIdx.push_back( iLoc );
            recvRanks.push_back( recvFrom );
            ++recvCounts[recvFrom];
        }
    }

    // Construct the send and recv displacements from the counts
    const Int totalSend = Scan( sendCounts, sendDispls );
    const Int totalRecv = Scan( recvCounts, recvDispls );
    DEBUG_ONLY(
      if( totalSend != totalRecv )
          LogicError
          ("Send and recv counts do not match: send=",totalSend,", recv=",
           totalRecv);
    )
}

} // namespace El
