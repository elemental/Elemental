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
void PermuteRows
(       Matrix<T>& A,
  const Matrix<Int>& perm,
  const Matrix<Int>& invPerm )
{
    const Int b = perm.Height();
    DEBUG_ONLY(
      CSE cse("PermuteRows");
      if( A.Height() < b || b != invPerm.Height() )
          LogicError
          ("perm and invPerm must be vectors of equal length that are not "
           "taller than A.");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    const Int* permBuf = perm.LockedBuffer();
    const Int* invPermBuf = invPerm.LockedBuffer();

    // Make a copy of the first b rows
    auto ARowPanView = A( IR(0,b), IR(0,n) );
    Matrix<T> ARowPanTrans;
    Transpose( ARowPanView, ARowPanTrans );

    // Make a copy of the preimage rows
    Matrix<T> APreimageTrans( n, b );

    T* ABuf = A.Buffer();
    T* APreBuf = APreimageTrans.Buffer(); 
    T* ARowPanBuf = ARowPanTrans.Buffer();
    const Int ALDim = A.LDim();
    const Int APreLDim = APreimageTrans.LDim();
    const Int ARowPanLDim = ARowPanTrans.LDim();

    for( Int i=0; i<b; ++i ) 
    {
        const Int iPre = permBuf[i];
        if( iPre >= b )
            for( Int j=0; j<n; ++j )
                APreBuf[j+i*APreLDim] = ABuf[iPre+j*ALDim];
    }

    // Apply the permutations
    for( Int i=0; i<b; ++i )
    {
        const Int iPre = permBuf[i];
        const Int iPost = invPermBuf[i];
        // Move row[i] into row[image[i]]
        for( Int j=0; j<n; ++j )
            ABuf[iPost+j*ALDim] = ARowPanBuf[j+i*ARowPanLDim];
        if( iPre >= b )
        {
            // Move row[preimage[i]] into row[i]
            for( Int j=0; j<n; ++j )
                ABuf[i+j*ALDim] = APreBuf[j+i*APreLDim];
        }
    }
}

template<typename T> 
void PermuteRows( Matrix<T>& A, const Matrix<Int>& perm )
{
    DEBUG_ONLY(CSE cse("PermuteRows"))
    Matrix<Int> invPerm;
    InvertPermutation( perm, invPerm );
    PermuteRows( A, perm, invPerm );
}

template<typename T> 
void InversePermuteRows( Matrix<T>& A, const Matrix<Int>& invPerm )
{
    DEBUG_ONLY(CSE cse("InversePermuteRows"))
    Matrix<Int> perm;
    InvertPermutation( invPerm, perm );
    PermuteRows( A, perm, invPerm );
}

template<typename T>
void PermuteRows
(       AbstractDistMatrix<T>& A,
  const PermutationMeta& oldMeta,
  bool inverse )
{
    DEBUG_ONLY(
      CSE cse("PermuteRows");
      if( A.ColComm() != oldMeta.comm )
          LogicError("Invalid communicator in metadata");
      if( A.ColAlign() != oldMeta.align )
          LogicError("Invalid alignment in metadata");
    )
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;

    T* ABuf = A.Buffer();
    const Int ALDim = A.LDim();

    const Int localWidth = A.LocalWidth();
    PermutationMeta meta = oldMeta;
    meta.ScaleUp( localWidth );

    if( inverse )
    {
        // Fill vectors with the send data
        auto offsets = meta.recvDispls;
        const int totalSend = meta.TotalRecv();
        vector<T> sendData;
        FastResize( sendData, mpi::Pad(totalSend) );
        const int numSends = meta.recvIdx.size();
        for( int send=0; send<numSends; ++send )
        {
            const int iLoc = meta.recvIdx[send];
            const int rank = meta.recvRanks[send];
        
            StridedMemCopy
            ( &sendData[offsets[rank]], 1, &ABuf[iLoc], ALDim, localWidth );
            offsets[rank] += localWidth;
        }

        // Communicate all pivot rows
        const int totalRecv = meta.TotalSend();
        vector<T> recvData;
        FastResize( recvData, mpi::Pad(totalRecv) );
        mpi::AllToAll
        ( sendData.data(), meta.recvCounts.data(), meta.recvDispls.data(),
          recvData.data(), meta.sendCounts.data(), meta.sendDispls.data(),
          meta.comm );

        // Unpack the recv data
        offsets = meta.sendDispls;
        const int numRecvs = meta.sendIdx.size();
        for( int recv=0; recv<numRecvs; ++recv )
        {
            const int iLoc = meta.sendIdx[recv];
            const int rank = meta.sendRanks[recv];
            StridedMemCopy
            ( &ABuf[iLoc], ALDim, &recvData[offsets[rank]], 1,localWidth );
            offsets[rank] += localWidth;
        }
    }
    else
    {
        // Fill vectors with the send data
        auto offsets = meta.sendDispls;
        const int totalSend = meta.TotalSend();
        vector<T> sendData;
        FastResize( sendData, mpi::Pad(totalSend) );
        const int numSends = meta.sendIdx.size();
        for( int send=0; send<numSends; ++send )
        {
            const int iLoc = meta.sendIdx[send];
            const int rank = meta.sendRanks[send];
        
            StridedMemCopy
            ( &sendData[offsets[rank]], 1, &ABuf[iLoc], ALDim, localWidth );
            offsets[rank] += localWidth;
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
            const int iLoc = meta.recvIdx[recv];
            const int rank = meta.recvRanks[recv];
            StridedMemCopy
            ( &ABuf[iLoc], ALDim, &recvData[offsets[rank]], 1,localWidth );
            offsets[rank] += localWidth;
        }
    }
}

template<typename T>
void PermuteRows
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& permPre,
  const AbstractDistMatrix<Int>& invPermPre )
{
    DEBUG_ONLY(CSE cse("PermuteRows"))

    DistMatrixReadProxy<Int,Int,STAR,STAR>
      permProx( permPre ),
      invPermProx( invPermPre );
    auto& perm = permProx.GetLocked();
    auto& invPerm = invPermProx.GetLocked();
   
    if( A.Participating() )
    {
        PermutationMeta meta( perm, invPerm, A.ColAlign(), A.ColComm() );
        PermuteRows( A, meta );
    }
}

template<typename T>
void PermuteRows
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& perm )
{
    DEBUG_ONLY(CSE cse("PermuteRows"))
    const Grid& g = A.Grid();
    DistMatrix<Int,STAR,STAR> invPerm(g);
    InvertPermutation( perm, invPerm );
    PermuteRows( A, perm, invPerm );
}

template<typename T>
void InversePermuteRows
(       AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<Int>& invPerm )
{
    DEBUG_ONLY(CSE cse("InversePermuteRows"))
    const Grid& g = A.Grid();
    DistMatrix<Int,STAR,STAR> perm(g);
    InvertPermutation( invPerm, perm );
    PermuteRows( A, perm, invPerm );
}

#define PROTO(T) \
  template void PermuteRows \
  (       Matrix<T>& A, \
    const Matrix<Int>& perm ); \
  template void PermuteRows \
  (       AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<Int>& perm ); \
  template void InversePermuteRows \
   (       Matrix<T>& A, \
     const Matrix<Int>& perm ); \
  template void InversePermuteRows \
  (       AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<Int>& perm ); \
  template void PermuteRows \
  (       Matrix<T>& A, \
    const Matrix<Int>& perm, \
    const Matrix<Int>& invPerm ); \
  template void PermuteRows \
  (       AbstractDistMatrix<T>& A, \
    const AbstractDistMatrix<Int>& perm, \
    const AbstractDistMatrix<Int>& invPerm ); \
  template void PermuteRows \
  (       AbstractDistMatrix<T>& A, \
    const PermutationMeta& oldMeta, \
    bool inverse );

#include "El/macros/Instantiate.h"

} // namespace El
