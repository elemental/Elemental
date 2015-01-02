/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

PermutationMeta::PermutationMeta
( const AbstractDistMatrix<Int>& perm,
  const AbstractDistMatrix<Int>& invPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("PermutationMeta::PermutationMeta");
        AssertSameGrids( perm, invPerm );
        if( perm.ColDist() != invPerm.ColDist() )
            LogicError("perm and invPerm must have the same column dist");
        if( perm.RowDist() != invPerm.RowDist() )
            LogicError("perm and invPerm must have the same row dist");
        if( perm.ColAlign() != invPerm.ColAlign() )
            LogicError("perm and invPerm must align");
        if( perm.RowDist() != Collect( perm.ColDist() ) )
            LogicError("permutations must have column distributions");
    )
    comm = perm.ColComm();
    align = perm.ColAlign();
    const Int b = perm.Height();
    const Int stride = perm.ColStride();

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
    sendCounts.resize( stride, 0 );
    recvCounts.resize( stride, 0 );

    for( Int i=0; i<b; ++i )
    {
        const Int preVal = perm.Get(i,0);
        const Int postVal = invPerm.Get(i,0);        

        // Handle sends 
        if( invPerm.IsLocalRow(i) )
        {
            const Int iLoc = invPerm.LocalRow( i );
            const Int sendTo = invPerm.RowOwner( postVal );
            sendIdx.push_back( iLoc );
            sendRanks.push_back( sendTo );
            ++sendCounts[sendTo];
        }
        if( preVal >= b && invPerm.IsLocalRow(preVal) )
        {
            const Int iLoc = invPerm.LocalRow( preVal );
            const Int sendTo = invPerm.RowOwner( i) ;
            sendIdx.push_back( iLoc );
            sendRanks.push_back( sendTo );
            ++sendCounts[sendTo];
        }

        // Handle recvs
        if( invPerm.IsLocalRow(postVal) )
        {
            const Int iLoc = invPerm.LocalRow( postVal );
            const Int recvFrom = invPerm.RowOwner( i );
            recvIdx.push_back( iLoc );
            recvRanks.push_back( recvFrom );
            ++recvCounts[recvFrom];
        }
        if( preVal >= b && invPerm.IsLocalRow(i) )
        {
            const Int iLoc = invPerm.LocalRow( i );
            const Int recvFrom = invPerm.RowOwner( preVal );
            recvIdx.push_back( iLoc );
            recvRanks.push_back( recvFrom );
            ++recvCounts[recvFrom];
        }
    }

    // Construct the send and recv displacements from the counts
    sendDispls.resize( stride );
    recvDispls.resize( stride );
    Int totalSend=0, totalRecv=0;
    for( Int i=0; i<stride; ++i )
    {
        sendDispls[i] = totalSend;
        recvDispls[i] = totalRecv;
        totalSend += sendCounts[i];
        totalRecv += recvCounts[i];
    }
    DEBUG_ONLY(
        if( totalSend != totalRecv )
            LogicError
            ("Send and recv counts do not match: send=",totalSend,", recv=",
             totalRecv);
    )
}

} // namespace El
