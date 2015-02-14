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
void PermuteCols
( Matrix<T>& A, const Matrix<Int>& perm, const Matrix<Int>& invPerm )
{
    const Int b = perm.Height();
    DEBUG_ONLY(
        CallStackEntry cse("PermuteCols");
        if( A.Width() < b || b != invPerm.Height() )
            LogicError
            ("perm and invPerm must be vectors of equal length that are not "
             "wider than A.");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    // TODO: Optimize this routine

    // Make a copy of the first b columns
    auto AColPanView = A( IR(0,m), IR(0,b) );
    auto AColPanCopy( AColPanView );

    // Make a copy of the preimage columns
    Matrix<T> APreimageCopy( m, b );
    for( Int j=0; j<b; ++j )
    {
        const Int jPre = perm.Get(j,0);
        if( jPre >= b )
            MemCopy( APreimageCopy.Buffer(0,j), A.LockedBuffer(0,jPre), m );
    }

    // Apply the permutations
    for( Int j=0; j<b; ++j )
    {
        const Int jPre = perm.Get(j,0);
        const Int jPost = invPerm.Get(j,0);
        // Move row[j] into row[jPost]
        MemCopy( A.Buffer(0,jPost), AColPanCopy.LockedBuffer(0,j), m );
        // Move row[jPre] into row[j]
        if( jPre >= b )
            MemCopy( A.Buffer(0,j), APreimageCopy.LockedBuffer(0,j), m );
    }
}

template<typename T> 
void PermuteCols( Matrix<T>& A, const Matrix<Int>& perm )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteCols"))
    Matrix<Int> invPerm;
    InvertPermutation( perm, invPerm );
    PermuteCols( A, perm, invPerm );
}

template<typename T> 
void InversePermuteCols( Matrix<T>& A, const Matrix<Int>& invPerm )
{
    DEBUG_ONLY(CallStackEntry cse("InversePermuteCols"))
    Matrix<Int> perm;
    InvertPermutation( invPerm, perm );
    PermuteCols( A, perm, invPerm );
}

template<typename T> 
void PermuteCols( AbstractDistMatrix<T>& A, const PermutationMeta& oldMeta )
{
    DEBUG_ONLY(
        CallStackEntry cse("PermuteCols");
        if( A.RowComm() != oldMeta.comm )
            LogicError("Invalid communicator in metadata");
        if( A.RowAlign() != oldMeta.align )
            LogicError("Invalid alignment in metadata");
    )
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;

    const Int localHeight = A.LocalHeight();
    PermutationMeta meta = oldMeta;
    meta.ScaleUp( localHeight ); 
    // Fill vectors with the send data
    auto offsets = meta.sendDispls;
    const int totalSend = meta.TotalSend();
    vector<T> sendData( mpi::Pad(totalSend) );
    const int numSends = meta.sendIdx.size();
    for( int send=0; send<numSends; ++send )
    {
        const int jLoc = meta.sendIdx[send];
        const int rank = meta.sendRanks[send];
        MemCopy
        ( &sendData[offsets[rank]], A.LockedBuffer(0,jLoc), localHeight );
        offsets[rank] += localHeight;
    }

    // Communicate all pivot rows
    const int totalRecv = meta.TotalRecv();
    vector<T> recvData( mpi::Pad(totalRecv) );
    mpi::AllToAll
    ( sendData.data(), meta.sendCounts.data(), meta.sendDispls.data(),
      recvData.data(), meta.recvCounts.data(), meta.recvDispls.data(), 
      meta.comm );

    // Unpack the recv data
    offsets = meta.recvDispls;
    const int numRecvs = meta.recvIdx.size();
    for( int recv=0; recv<numRecvs; ++recv )
    {
        const int jLoc = meta.recvIdx[recv];
        const int rank = meta.recvRanks[recv];
        MemCopy
        ( A.Buffer(0,jLoc), &recvData[offsets[rank]], localHeight );
        offsets[rank] += localHeight;
    }
}

template<typename T,Dist U,Dist V>
void PermuteCols
( DistMatrix<T,U,V>& A, 
  const AbstractDistMatrix<Int>& permPre, 
  const AbstractDistMatrix<Int>& invPermPre )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteCols"))

    ProxyCtrl ctrl;
    ctrl.rootConstrain = true;
    ctrl.colConstrain = true;
    ctrl.root = A.Root();
    ctrl.colAlign = A.RowAlign();

    auto permPtr    = ReadProxy<Int,V,Collect<U>()>( &permPre,    ctrl );
    auto invPermPtr = ReadProxy<Int,V,Collect<U>()>( &invPermPre, ctrl );

    if( A.Participating() )
    {
        PermutationMeta meta( *permPtr, *invPermPtr );
        PermuteCols( A, meta );
    }
}

template<typename T>
void PermuteCols
( AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<Int>& perm,
  const AbstractDistMatrix<Int>& invPerm )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteCols"))
    const Dist U = APre.ColDist();
    const Dist V = APre.RowDist();
    #define GUARD(CDIST,RDIST) U == CDIST && V == RDIST
    #define PAYLOAD(CDIST,RDIST) \
        auto& A = dynamic_cast<DistMatrix<T,CDIST,RDIST>&>(APre); \
        PermuteCols( A, perm, invPerm );
    #include "El/macros/GuardAndPayload.h"
}

template<typename T>
void PermuteCols
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& perm )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteCols"))
    const Grid& g = A.Grid();
    DistMatrix<Int,VC,STAR> invPerm(g);
    InvertPermutation( perm, invPerm );
    PermuteCols( A, perm, invPerm );
}

template<typename T>
void InversePermuteCols
( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& invPerm )
{
    DEBUG_ONLY(CallStackEntry cse("InversePermuteCols"))
    const Grid& g = A.Grid();
    DistMatrix<Int,VC,STAR> perm(g);
    InvertPermutation( invPerm, perm );
    PermuteCols( A, perm, invPerm );
}

#define PROTO(T) \
  template void PermuteCols( Matrix<T>& A, const Matrix<Int>& perm ); \
  template void PermuteCols \
  ( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& perm ); \
  template void InversePermuteCols( Matrix<T>& A, const Matrix<Int>& perm ); \
  template void InversePermuteCols \
  ( AbstractDistMatrix<T>& A, const AbstractDistMatrix<Int>& perm ); \
  template void PermuteCols \
  ( Matrix<T>& A, const Matrix<Int>& perm, const Matrix<Int>& invPerm ); \
  template void PermuteCols \
  ( AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<Int>& perm, \
    const AbstractDistMatrix<Int>& invPerm ); \
  template void PermuteCols \
  ( AbstractDistMatrix<T>& A, const PermutationMeta& oldMeta );

#include "El/macros/Instantiate.h"

} // namespace El
