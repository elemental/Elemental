/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_PERMUTECOLS_HPP
#define ELEM_LAPACK_PERMUTECOLS_HPP

#include "./PermutationMeta.hpp"
#include "./InvertPermutation.hpp"

namespace elem {

template<typename T> 
inline void
PermuteCols( Matrix<T>& A, const Matrix<Int>& perm, const Matrix<Int>& invPerm )
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
inline void
PermuteCols( Matrix<T>& A, const Matrix<Int>& perm )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteCols"))
    Matrix<Int> invPerm;
    InvertPermutation( perm, invPerm );
    PermuteCols( A, perm, invPerm );
}

template<typename T> 
inline void
InversePermuteCols( Matrix<T>& A, const Matrix<Int>& invPerm )
{
    DEBUG_ONLY(CallStackEntry cse("InversePermuteCols"))
    Matrix<Int> perm;
    InvertPermutation( invPerm, perm );
    PermuteCols( A, perm, invPerm );
}

template<typename T,Dist U,Dist V> 
inline void
PermuteCols( DistMatrix<T,U,V>& A, const PermutationMeta& oldMeta )
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
inline void
PermuteCols
( DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,STAR>& perm, 
  const DistMatrix<Int,UPerm,STAR>& invPerm )
{
    DEBUG_ONLY(
        CallStackEntry cse("PermuteCols");
        if( perm.ColAlign() != invPerm.ColAlign() )
            LogicError("misaligned perm and invPerm");
    )
    const Grid& g = A.Grid();
    DistMatrix<Int,V,STAR> perm_V_STAR(g), invPerm_V_STAR(g);
    if( V == UPerm && A.RowAlign() == perm.ColAlign() )
    {
        perm_V_STAR = LockedView( perm );
        invPerm_V_STAR = LockedView( invPerm );
    }
    else
    {
        perm_V_STAR.AlignWith( A );
        perm_V_STAR = perm;
        invPerm_V_STAR.AlignWith( A );
        invPerm_V_STAR = invPerm;
    }

    if( A.Participating() )
    {
        PermutationMeta meta( perm_V_STAR, invPerm_V_STAR );
        PermuteCols( A, meta );
    }
}

template<typename T,Dist U,Dist V,Dist UPerm>
inline void
PermuteCols
( DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,STAR>& perm )
{
    DEBUG_ONLY(CallStackEntry cse("PermuteCols"))
    const Grid& g = A.Grid();
    DistMatrix<Int,V,STAR> perm_V_STAR(g), invPerm_V_STAR(g);
    perm_V_STAR.AlignWith( A );
    perm_V_STAR = perm;
    InvertPermutation( perm_V_STAR, invPerm_V_STAR );
    PermuteCols( A, perm_V_STAR, invPerm_V_STAR );
}

template<typename T,Dist U,Dist V,Dist UPerm>
inline void
InversePermuteCols
( DistMatrix<T,U,V>& A, 
  const DistMatrix<Int,UPerm,STAR>& invPerm )
{
    DEBUG_ONLY(CallStackEntry cse("InversePermuteCols"))
    const Grid& g = A.Grid();
    DistMatrix<Int,V,STAR> perm_V_STAR(g), invPerm_V_STAR(g);
    invPerm_V_STAR.AlignWith( A );
    invPerm_V_STAR = invPerm;
    InvertPermutation( invPerm_V_STAR, perm_V_STAR );
    PermuteCols( A, perm_V_STAR, invPerm_V_STAR );
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_PERMUTECOLS_HPP
