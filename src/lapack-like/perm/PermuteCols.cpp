/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

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
    auto AColPanView = LockedView( A, 0, 0, m, b );
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

template<typename T,Dist U,Dist V> 
void PermuteCols( DistMatrix<T,U,V>& A, const PermutationMeta& oldMeta )
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
    std::vector<T> sendData( mpi::Pad(totalSend) );
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
    std::vector<T> recvData( mpi::Pad(totalRecv) );
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

template<typename T,Dist U,Dist V,Dist UPerm>
void PermuteCols
( DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& invPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("PermuteCols");
        if( perm.ColAlign() != invPerm.ColAlign() )
            LogicError("misaligned perm and invPerm");
    )
    const Grid& g = A.Grid();
    const Dist UGath = GatheredDist<U>();
    DistMatrix<Int,V,UGath> perm_V_UGath(g), invPerm_V_UGath(g);
    if( V == UPerm && A.RowAlign() == perm.ColAlign() )
    {
        perm_V_UGath = LockedView( perm );
        invPerm_V_UGath = LockedView( invPerm );
    }
    else
    {
        perm_V_UGath.AlignWith( A );
        perm_V_UGath = perm;
        invPerm_V_UGath.AlignWith( A );
        invPerm_V_UGath = invPerm;
    }

    if( A.Participating() )
    {
        PermutationMeta meta( perm_V_UGath, invPerm_V_UGath );
        PermuteCols( A, meta );
    }
}

template<typename T,Dist U,Dist V,Dist UPerm>
void PermuteCols
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteCols"))
    const Grid& g = A.Grid();
    const Dist UGath = GatheredDist<U>();
    DistMatrix<Int,V,UGath> perm_V_UGath(g), invPerm_V_UGath(g);
    perm_V_UGath.AlignWith( A );
    perm_V_UGath = perm;
    InvertPermutation( perm_V_UGath, invPerm_V_UGath );
    PermuteCols( A, perm_V_UGath, invPerm_V_UGath );
}

template<typename T,Dist U,Dist V,Dist UPerm>
void InversePermuteCols
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& invPerm )
{
    DEBUG_ONLY(CallStackEntry cse("InversePermuteCols"))
    const Grid& g = A.Grid();
    const Dist UGath = GatheredDist<U>();
    DistMatrix<Int,V,UGath> perm_V_UGath(g), invPerm_V_UGath(g);
    invPerm_V_UGath.AlignWith( A );
    invPerm_V_UGath = invPerm;
    InvertPermutation( invPerm_V_UGath, perm_V_UGath );
    PermuteCols( A, perm_V_UGath, invPerm_V_UGath );
}

#define PROTO_DIST_INTERNAL(T,U,V,UPERM) \
  template void PermuteCols \
  (       DistMatrix<T,U,V>& A, \
    const DistMatrix<Int,UPERM,GatheredDist<UPERM>()>& perm ); \
  template void InversePermuteCols \
  (       DistMatrix<T,U,V>& A, \
    const DistMatrix<Int,UPERM,GatheredDist<UPERM>()>& perm ); \
  template void PermuteCols \
  ( DistMatrix<T,U,V>& A, \
    const DistMatrix<Int,UPERM,GatheredDist<UPERM>()>& perm, \
    const DistMatrix<Int,UPERM,GatheredDist<UPERM>()>& invPerm );

#define PROTO_DIST(T,U,V) \
  PROTO_DIST_INTERNAL(T,U,V,CIRC) \
  PROTO_DIST_INTERNAL(T,U,V,MC  ) \
  PROTO_DIST_INTERNAL(T,U,V,MD  ) \
  PROTO_DIST_INTERNAL(T,U,V,MR  ) \
  PROTO_DIST_INTERNAL(T,U,V,STAR) \
  PROTO_DIST_INTERNAL(T,U,V,VC  ) \
  PROTO_DIST_INTERNAL(T,U,V,VR  ) \
  template void PermuteCols \
  ( DistMatrix<T,U,V>& A, const PermutationMeta& oldMeta );

#define PROTO(T) \
  template void PermuteCols( Matrix<T>& A, const Matrix<Int>& perm ); \
  template void InversePermuteCols( Matrix<T>& A, const Matrix<Int>& perm ); \
  template void PermuteCols \
  ( Matrix<T>& A, const Matrix<Int>& perm, const Matrix<Int>& invPerm ); \
  PROTO_DIST(T,CIRC,CIRC) \
  PROTO_DIST(T,MC,  MR  ) \
  PROTO_DIST(T,MC,  STAR) \
  PROTO_DIST(T,MD,  STAR) \
  PROTO_DIST(T,MR,  MC  ) \
  PROTO_DIST(T,MR,  STAR) \
  PROTO_DIST(T,STAR,MC  ) \
  PROTO_DIST(T,STAR,MD  ) \
  PROTO_DIST(T,STAR,MR  ) \
  PROTO_DIST(T,STAR,STAR) \
  PROTO_DIST(T,STAR,VC  ) \
  PROTO_DIST(T,STAR,VR  ) \
  PROTO_DIST(T,VC,  STAR) \
  PROTO_DIST(T,VR,  STAR)

PROTO(Int)
PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
