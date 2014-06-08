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
void PermuteRows
( Matrix<T>& A, const Matrix<Int>& perm, const Matrix<Int>& invPerm )
{
    const Int b = perm.Height();
    DEBUG_ONLY(
        CallStackEntry cse("PermuteRows");
        if( A.Height() < b || b != invPerm.Height() )
            LogicError
            ("perm and invPerm must be vectors of equal length that are not "
             "taller than A.");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    if( m == 0 || n == 0 )
        return;

    // TODO: Optimize this routine

    // Make a copy of the first b rows
    auto ARowPanView = LockedView( A, 0, 0, b, n );
    auto ARowPanCopy( ARowPanView );

    // Make a copy of the preimage rows
    Matrix<T> APreimageCopy( b, n );
    for( Int i=0; i<b; ++i ) 
    {
        const Int iPre = perm.Get(i,0);
        if( iPre >= b )
            for( Int j=0; j<n; ++j )
                APreimageCopy.Set(i,j,A.Get(iPre,j));
    }

    // Apply the permutations
    for( Int i=0; i<b; ++i )
    {
        const Int iPre = perm.Get(i,0);
        const Int iPost = invPerm.Get(i,0);
        // Move row[i] into row[image[i]]
        for( Int j=0; j<n; ++j )
            A.Set(iPost,j,ARowPanCopy.Get(i,j));
        if( iPre >= b )
        {
            // Move row[preimage[i]] into row[i]
            for( Int j=0; j<n; ++j )
                A.Set(i,j,APreimageCopy.Get(i,j));
        }
    }
}

template<typename T> 
void PermuteRows( Matrix<T>& A, const Matrix<Int>& perm )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteRows"))
    Matrix<Int> invPerm;
    InvertPermutation( perm, invPerm );
    PermuteRows( A, perm, invPerm );
}

template<typename T> 
void InversePermuteRows( Matrix<T>& A, const Matrix<Int>& invPerm )
{
    DEBUG_ONLY(CallStackEntry cse("InversePermuteRows"))
    Matrix<Int> perm;
    InvertPermutation( invPerm, perm );
    PermuteRows( A, perm, invPerm );
}

template<typename T,Dist U,Dist V>
void PermuteRows( DistMatrix<T,U,V>& A, const PermutationMeta& oldMeta )
{
    DEBUG_ONLY(
        CallStackEntry cse("PermuteRows");
        if( A.ColComm() != oldMeta.comm )
            LogicError("Invalid communicator in metadata");
        if( A.ColAlign() != oldMeta.align )
            LogicError("Invalid alignment in metadata");
    )
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;

    const Int localWidth = A.LocalWidth();
    const Int ldim = A.LDim();
    PermutationMeta meta = oldMeta;
    meta.ScaleUp( localWidth );

    // Fill vectors with the send data
    auto offsets = meta.sendDispls;
    const int totalSend = meta.TotalSend();
    std::vector<T> sendData( mpi::Pad(totalSend) );
    const int numSends = meta.sendIdx.size();
    for( int send=0; send<numSends; ++send )
    {
        const int iLoc = meta.sendIdx[send];
        const int rank = meta.sendRanks[send];
        
        StridedMemCopy
        ( &sendData[offsets[rank]], 1,
          A.LockedBuffer(iLoc,0),   ldim, localWidth );
        offsets[rank] += localWidth;
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
        const int iLoc = meta.recvIdx[recv];
        const int rank = meta.recvRanks[recv];
        StridedMemCopy
        ( A.Buffer(iLoc,0),         ldim,
          &recvData[offsets[rank]], 1,    localWidth );
        offsets[rank] += localWidth;
    }
}

template<typename T,Dist U,Dist V,Dist UPerm>
void PermuteRows
( DistMatrix<T,U,V>& A,
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm,
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& invPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("PermuteRows");
        if( perm.ColAlign() != invPerm.ColAlign() )
            LogicError("misaligned perm and invPerm");
    )
    const Grid& g = A.Grid();
    const Dist VGath = GatheredDist<V>();
    DistMatrix<Int,U,VGath> perm_U_VGath(g), invPerm_U_VGath(g);
    if( U == UPerm && A.ColAlign() == perm.ColAlign() )
    {
        perm_U_VGath = LockedView( perm );
        invPerm_U_VGath = LockedView( invPerm );
    }
    else
    {
        perm_U_VGath.AlignWith( A );
        perm_U_VGath = perm;
        invPerm_U_VGath.AlignWith( A );
        invPerm_U_VGath = invPerm;
    }
        
    if( A.Participating() )
    {
        PermutationMeta meta( perm_U_VGath, invPerm_U_VGath );
        PermuteRows( A, meta );
    }
}

template<typename T,Dist U,Dist V,Dist UPerm>
void PermuteRows
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& perm )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteRows"))
    const Grid& g = A.Grid();
    const Dist VGath = GatheredDist<V>();
    DistMatrix<Int,U,VGath> perm_U_VGath(g), invPerm_U_VGath(g);
    perm_U_VGath.AlignWith( A );
    perm_U_VGath = perm;
    InvertPermutation( perm_U_VGath, invPerm_U_VGath );
    PermuteRows( A, perm_U_VGath, invPerm_U_VGath );
}

template<typename T,Dist U,Dist V,Dist UPerm>
void InversePermuteRows
(       DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,GatheredDist<UPerm>()>& invPerm )
{
    DEBUG_ONLY(CallStackEntry cse("InversePermuteRows"))
    const Grid& g = A.Grid();
    const Dist VGath = GatheredDist<V>();
    DistMatrix<Int,U,VGath> perm_U_VGath(g), invPerm_U_VGath(g);
    invPerm_U_VGath.AlignWith( A );
    invPerm_U_VGath = invPerm;
    InvertPermutation( invPerm_U_VGath, perm_U_VGath );
    PermuteRows( A, perm_U_VGath, invPerm_U_VGath );
}

#define PROTO_DIST_INTERNAL(T,U,V,UPERM) \
  template void PermuteRows \
  (       DistMatrix<T,U,V>& A, \
    const DistMatrix<Int,UPERM,GatheredDist<UPERM>()>& perm ); \
  template void InversePermuteRows \
  (       DistMatrix<T,U,V>& A, \
    const DistMatrix<Int,UPERM,GatheredDist<UPERM>()>& perm ); \
  template void PermuteRows \
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
  template void PermuteRows \
  ( DistMatrix<T,U,V>& A, const PermutationMeta& oldMeta );

#define PROTO(T) \
  template void PermuteRows( Matrix<T>& A, const Matrix<Int>& perm ); \
  template void InversePermuteRows( Matrix<T>& A, const Matrix<Int>& perm ); \
  template void PermuteRows \
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
