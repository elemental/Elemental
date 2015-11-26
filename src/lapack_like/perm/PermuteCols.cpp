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
(       Matrix<T>& A,
  const Matrix<Int>& perm,
  const Matrix<Int>& invPerm )
{
    const Int b = perm.Height();
    DEBUG_ONLY(
      CSE cse("PermuteCols");
      if( A.Width() < b || b != invPerm.Height() )
          LogicError
          ("perm and invPerm must be vectors of equal length that are not "
           "wider than A.");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    // Make a copy of the first b columns
    auto AColPanView = A( IR(0,m), IR(0,b) );
    auto AColPanCopy( AColPanView );

    // Make a copy of the preimage columns
    Matrix<T> APreimageCopy( m, b );

    const Int* permBuf = perm.LockedBuffer();
    const Int* invPermBuf = invPerm.LockedBuffer();
          T* ABuf = A.Buffer();
          T* APreBuf = APreimageCopy.Buffer();
          T* AColPanBuf = AColPanCopy.Buffer();
    const Int ALDim = A.LDim();
    const Int APreLDim = APreimageCopy.LDim();
    const Int AColPanLDim = AColPanCopy.LDim();

    for( Int j=0; j<b; ++j )
    {
        const Int jPre = permBuf[j];
        if( jPre >= b )
            MemCopy( &APreBuf[j*APreLDim], &ABuf[jPre*ALDim], m );
    }

    // Apply the permutations
    for( Int j=0; j<b; ++j )
    {
        const Int jPre = permBuf[j];
        const Int jPost = invPermBuf[j];
        // Move row[j] into row[jPost]
        MemCopy( &ABuf[jPost*ALDim], &AColPanBuf[j*AColPanLDim], m );
        // Move row[jPre] into row[j]
        if( jPre >= b )
            MemCopy( &ABuf[j*ALDim], &APreBuf[j*APreLDim], m );
    }
}

template<typename T> 
void PermuteCols( Matrix<T>& A, const Matrix<Int>& perm )
{
    DEBUG_ONLY(CSE cse("PermuteCols"))
    Matrix<Int> invPerm;
    InvertPermutation( perm, invPerm );
    PermuteCols( A, perm, invPerm );
}

template<typename T> 
void InversePermuteCols( Matrix<T>& A, const Matrix<Int>& invPerm )
{
    DEBUG_ONLY(CSE cse("InversePermuteCols"))
    Matrix<Int> perm;
    InvertPermutation( invPerm, perm );
    PermuteCols( A, perm, invPerm );
}

template<typename T> 
void PermuteCols( ElementalMatrix<T>& A, const PermutationMeta& oldMeta )
{
    DEBUG_ONLY(
      CSE cse("PermuteCols");
      if( A.RowComm() != oldMeta.comm )
          LogicError("Invalid communicator in metadata");
      if( A.RowAlign() != oldMeta.align )
          LogicError("Invalid alignment in metadata");
    )
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;

    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    const Int localHeight = A.LocalHeight();
    PermutationMeta meta = oldMeta;
    meta.ScaleUp( localHeight ); 
    // Fill vectors with the send data
    auto offsets = meta.sendDispls;
    const int totalSend = meta.TotalSend();
    vector<T> sendData;
    FastResize( sendData, mpi::Pad(totalSend) );
    const int numSends = meta.sendIdx.size();
    for( int send=0; send<numSends; ++send )
    {
        const int jLoc = meta.sendIdx[send];
        const int rank = meta.sendRanks[send];
        MemCopy( &sendData[offsets[rank]], &ABuf[jLoc*ALDim], localHeight );
        offsets[rank] += localHeight;
    }

    // Communicate all pivot rows
    const int totalRecv = meta.TotalRecv();
    vector<T> recvData;
    FastResize( recvData, mpi::Pad(totalRecv) );
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
        MemCopy( &ABuf[jLoc*ALDim], &recvData[offsets[rank]], localHeight );
        offsets[rank] += localHeight;
    }
}

template<typename T>
void PermuteCols
(       ElementalMatrix<T>& A, 
  const ElementalMatrix<Int>& permPre, 
  const ElementalMatrix<Int>& invPermPre )
{
    DEBUG_ONLY(CSE cse("PermuteCols"))

    DistMatrixReadProxy<Int,Int,STAR,STAR>
      permProx( permPre ),
      invPermProx( invPermPre );
    auto& perm = permProx.GetLocked();
    auto& invPerm = invPermProx.GetLocked();

    if( A.Participating() )
    {
        PermutationMeta meta( perm, invPerm, A.RowAlign(), A.RowComm() );
        PermuteCols( A, meta );
    }
}

template<typename T>
void PermuteCols
(       ElementalMatrix<T>& A,
  const ElementalMatrix<Int>& perm )
{
    DEBUG_ONLY(CSE cse("PermuteCols"))
    const Grid& g = A.Grid();
    DistMatrix<Int,STAR,STAR> invPerm(g);
    InvertPermutation( perm, invPerm );
    PermuteCols( A, perm, invPerm );
}

template<typename T>
void InversePermuteCols
(       ElementalMatrix<T>& A,
  const ElementalMatrix<Int>& invPerm )
{
    DEBUG_ONLY(CSE cse("InversePermuteCols"))
    const Grid& g = A.Grid();
    DistMatrix<Int,STAR,STAR> perm(g);
    InvertPermutation( invPerm, perm );
    PermuteCols( A, perm, invPerm );
}

#define PROTO(T) \
  template void PermuteCols \
  (       Matrix<T>& A, \
    const Matrix<Int>& perm ); \
  template void PermuteCols \
  (       ElementalMatrix<T>& A, \
    const ElementalMatrix<Int>& perm ); \
  template void InversePermuteCols \
  (       Matrix<T>& A, \
    const Matrix<Int>& perm ); \
  template void InversePermuteCols \
  (       ElementalMatrix<T>& A, \
    const ElementalMatrix<Int>& perm ); \
  template void PermuteCols \
  (       Matrix<T>& A, \
    const Matrix<Int>& perm, \
    const Matrix<Int>& invPerm ); \
  template void PermuteCols \
  (       ElementalMatrix<T>& A, \
    const ElementalMatrix<Int>& perm, \
    const ElementalMatrix<Int>& invPerm ); \
  template void PermuteCols \
  (       ElementalMatrix<T>& A, \
    const PermutationMeta& oldMeta );

#include "El/macros/Instantiate.h"

} // namespace El
