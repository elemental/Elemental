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
    DEBUG_ONLY(CSE cse("Full"))
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
    DEBUG_ONLY(CSE cse("Full"))

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

    // Pack the data
    // -------------
    vector<int> sendOffs;
    const int totalSend = Scan( sendSizes, sendOffs );
    vector<Entry<T>> sendBuf(totalSend);
    auto offs = sendOffs;
    for( Int e=0; e<numLocalEntries; ++e )
    {
        const Int i = A.Row(e); 
        const Int j = A.Col(e);
        const T value = A.Value(e);
        const int owner = B.Owner( i, j );
        sendBuf[offs[owner]++] = Entry<T>{ i, j, value };
    }

    // Exchange and unpack
    // -------------------
    auto recvBuf = mpi::AllToAll( sendBuf, sendSizes, sendOffs, comm );
    for( auto& entry : recvBuf )
        B.Update( entry );
}

template<typename T>
Matrix<T> Full( const SparseMatrix<T>& A )
{
    DEBUG_ONLY(CSE cse("Full"))
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

#define EL_ENABLE_QUAD
#include "El/macros/Instantiate.h"

} // namespace El
