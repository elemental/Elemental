/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_BLAS_COPY_EXCHANGE_HPP
#define EL_BLAS_COPY_EXCHANGE_HPP

namespace El {
namespace copy {

template<typename T>
void Exchange
( const ElementalMatrix<T>& A, 
        ElementalMatrix<T>& B, 
  int sendRank, int recvRank, mpi::Comm comm )
{
    DEBUG_CSE
    DEBUG_ONLY(AssertSameGrids( A, B ))
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
        Copy( A.LockedMatrix(), B.Matrix() );
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
        vector<T> buf;
        FastResize( buf, sendSize );
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
        vector<T> buf;
        FastResize( buf, recvSize );
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
        vector<T> sendBuf;
        FastResize( sendBuf, sendSize );
        copy::util::InterleaveMatrix
        ( localHeightA, localWidthA,
          A.LockedBuffer(), 1, A.LDim(),
          sendBuf.data(),   1, localHeightA );

        // Exchange with the partner
        vector<T> recvBuf;
        FastResize( recvBuf, recvSize );
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

template<typename T,Dist U,Dist V>
void ColwiseVectorExchange
( const DistMatrix<T,ProductDist<U,V>(),STAR>& A,
        DistMatrix<T,ProductDist<V,U>(),STAR>& B )
{
    DEBUG_CSE
    AssertSameGrids( A, B );
    if( !B.Participating() )
        return;

    const Int distSize = A.DistSize();
    const Int colDiff = A.ColShift() - B.ColShift();
    const Int sendRankB = Mod( B.DistRank()+colDiff, distSize );
    const Int recvRankA = Mod( A.DistRank()-colDiff, distSize );
    const Int recvRankB =
      (recvRankA/A.PartialColStride())+
      (recvRankA%A.PartialColStride())*A.PartialUnionColStride();
    copy::Exchange( A, B, sendRankB, recvRankB, B.DistComm() );
}

template<typename T,Dist U,Dist V>
void RowwiseVectorExchange
( const DistMatrix<T,STAR,ProductDist<U,V>()>& A,
        DistMatrix<T,STAR,ProductDist<V,U>()>& B )
{
    DEBUG_CSE
    AssertSameGrids( A, B );
    if( !B.Participating() )
        return;

    const Int distSize = A.DistSize();
    const Int rowDiff = A.RowShift() - B.RowShift();
    const Int sendRankB = Mod( B.DistRank()+rowDiff, distSize );
    const Int recvRankA = Mod( A.DistRank()-rowDiff, distSize );
    const Int recvRankB =
      (recvRankA/A.PartialRowStride())+
      (recvRankA%A.PartialRowStride())*A.PartialUnionRowStride();
    copy::Exchange( A, B, sendRankB, recvRankB, B.DistComm() );
}

} // namespace copy
} // namespace El

#endif // ifndef EL_BLAS_COPY_EXCHANGE_HPP
