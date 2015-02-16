/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace copy {

template<typename T>
void Exchange
( const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& B, 
  int sendRank, int recvRank, mpi::Comm comm )
{
    DEBUG_ONLY(
      CallStackEntry cse("copy::Exchange");
      AssertSameGrids( A, B );
    )
    const int myRank = mpi::Rank( comm );
    DEBUG_ONLY(
      if( myRank == sendRank && myRank != recvRank )
          LogicError("Sending to self but receiving from someone else");
      if( myRank != sendRank && myRank == recvRank )
          LogicError("Receiving from self but sending to someone else");
    )
    B.Resize( A.Height(), A.Width() );
    if( myRank == sendRank )
    {
        B.Matrix() = A.LockedMatrix();
        return;
    }

    const Int localHeightA = A.LocalHeight();
    const Int localHeightB = B.LocalHeight();
    const Int localWidthA = A.LocalWidth();
    const Int localWidthB = B.LocalWidth();
    const Int sendSize = localHeightA*localWidthA;
    const Int recvSize = localHeightB*localWidthB;

    const bool contigA = ( A.LocalHeight() == A.LDim() );
    const bool contigB = ( B.LocalHeight() == B.LDim() );

    if( contigA && contigB )
    {
        mpi::SendRecv
        ( A.LockedBuffer(), sendSize, sendRank,
          B.Buffer(),       recvSize, recvRank, comm );
    }
    else if( contigB )
    {
        // Pack A's data
        vector<T> buf( sendSize );
        copy::util::InterleaveMatrix
        ( localHeightA, localWidthA,
          A.LockedBuffer(), 1, A.LDim(),
          buf.data(),       1, localHeightA );

        // Exchange with the partner
        mpi::SendRecv
        ( buf.data(), sendSize, sendRank,
          B.Buffer(), recvSize, recvRank, comm );
    }
    else if( contigA )
    {
        // Exchange with the partner
        vector<T> buf( recvSize );
        mpi::SendRecv
        ( A.LockedBuffer(), sendSize, sendRank,
          buf.data(),       recvSize, recvRank, comm );

        // Unpack
        copy::util::InterleaveMatrix
        ( localHeightB, localWidthB,
          buf.data(), 1, localHeightB,
          B.Buffer(), 1, B.LDim() );
    }
    else
    {
        // Pack A's data
        vector<T> sendBuf( sendSize );
        copy::util::InterleaveMatrix
        ( localHeightA, localWidthA,
          A.LockedBuffer(), 1, A.LDim(),
          sendBuf.data(),   1, localHeightA );

        // Exchange with the partner
        vector<T> recvBuf( recvSize );
        mpi::SendRecv
        ( sendBuf.data(), sendSize, sendRank,
          recvBuf.data(), recvSize, recvRank, comm );

        // Unpack
        copy::util::InterleaveMatrix
        ( localHeightB, localWidthB,
          recvBuf.data(), 1, localHeightB,
          B.Buffer(),     1, B.LDim() );
    }
}

#define PROTO(T) \
  template void Exchange \
  ( const AbstractDistMatrix<T>& A, AbstractDistMatrix<T>& B, \
    int sendRank, int recvRank, mpi::Comm comm );

#include "El/macros/Instantiate.h"

} // namespace copy
} // namespace El
