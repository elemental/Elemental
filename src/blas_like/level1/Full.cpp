/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename T>
void Full( const SparseMatrix<T>& A, Matrix<T>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Full"))
    const Int m = A.Height();
    const Int n = A.Width();
    const Int numEntries = A.NumEntries();
    Zeros( B, m, n );
    for( Int e=0; e<numEntries; ++e )
        B.Set( A.Row(e), A.Col(e), A.Value(e) );
}

template<typename T>
void Full( const DistSparseMatrix<T>& A, AbstractDistMatrix<T>& BPre )
{
    DEBUG_ONLY(CallStackEntry cse("Full"))

    auto BPtr = WriteProxy<T,MC,MR>(&BPre);
    auto& B = *BPtr;
    const Int m = A.Height();
    const Int n = A.Width();
    Zeros( B, m, n );

    // Determine the metadata
    // ----------------------
    mpi::Comm comm = A.Comm();
    const int commSize = mpi::Size( comm );
    const Int numLocalEntries = A.NumLocalEntries();
    vector<int> sendSizes( commSize, 0 );
    for( Int e=0; e<numLocalEntries; ++e )
        ++sendSizes[ B.Owner( A.Row(e), A.Col(e) ) ];
    vector<int> recvSizes(commSize);
    mpi::AllToAll( sendSizes.data(), 1, recvSizes.data(), 1, comm );
    vector<int> sendOffs, recvOffs;
    const int totalSend = Scan( sendSizes, sendOffs );
    const int totalRecv = Scan( recvSizes, recvOffs );

    // Pack the data
    // -------------
    vector<Int> sSendBuf(totalSend), tSendBuf(totalSend);
    vector<T> vSendBuf(totalSend);
    auto offs = sendOffs;
    for( Int e=0; e<numLocalEntries; ++e )
    {
        const Int i = A.Row(e); 
        const Int j = A.Col(e);
        const T value = A.Value(e);
        const int owner = B.Owner( i, j );
        sSendBuf[offs[owner]] = i;
        tSendBuf[offs[owner]] = j;
        vSendBuf[offs[owner]] = value;
        ++offs[owner];
    }

    // Exchange the data
    // -----------------
    vector<Int> sRecvBuf(totalRecv), tRecvBuf(totalRecv);
    vector<T> vRecvBuf(totalRecv);
    mpi::AllToAll
    ( sSendBuf.data(), sendSizes.data(), sendOffs.data(),
      sRecvBuf.data(), recvSizes.data(), recvOffs.data(), comm );
    mpi::AllToAll
    ( tSendBuf.data(), sendSizes.data(), sendOffs.data(),
      tRecvBuf.data(), recvSizes.data(), recvOffs.data(), comm );
    mpi::AllToAll
    ( vSendBuf.data(), sendSizes.data(), sendOffs.data(),
      vRecvBuf.data(), recvSizes.data(), recvOffs.data(), comm );
    
    // Unpack the data
    // ---------------
    for( Int e=0; e<totalRecv; ++e )
        B.Update( sRecvBuf[e], tRecvBuf[e], vRecvBuf[e] );
}

template<typename T>
Matrix<T> Full( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("Full"))
    Matrix<T> B;
    Full( A, B );
    return B;
}

// NOTE: A DistSparseMatrix version does not exist since it is not yet clear
//       whether Elemental can currently handle creating a grid in such a 
//       routine without a memory leak


#define PROTO(T) \
  template void Full( const SparseMatrix<T>& A, Matrix<T>& B ); \
  template void Full \
  ( const DistSparseMatrix<T>& A, AbstractDistMatrix<T>& B ); \
  template Matrix<T> Full( const SparseMatrix<T>& A );

#include "El/macros/Instantiate.h"

} // namespace El
