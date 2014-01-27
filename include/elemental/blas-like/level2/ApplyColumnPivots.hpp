/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_APPLYCOLUMNPIVOTS_HPP
#define ELEM_APPLYCOLUMNPIVOTS_HPP

#include "./ComposePivots.hpp"

namespace elem {

template<typename F>
inline void
ApplyColumnPivots( Matrix<F>& A, const Matrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyColumnPivots");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( p.Height() > A.Width() )
            LogicError("p cannot be longer than width of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = p.Height();
    for( Int j=0; j<numPivots; ++j )
    {
        const Int k = p.Get(j,0);
        F* Aj = A.Buffer(0,j);
        F* Ak = A.Buffer(0,k);
        for( Int i=0; i<height; ++i )
        {
            F temp = Aj[i];
            Aj[i] = Ak[i];
            Ak[i] = temp;
        }
    }
}

template<typename F>
inline void
ApplyInverseColumnPivots( Matrix<F>& A, const Matrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyInverseColumnPivots");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( p.Height() > A.Width() )
            LogicError("p cannot be larger than width of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = p.Height();
    for( Int j=numPivots-1; j>=0; --j )
    {
        const Int k = p.Get(j,0);
        F* Aj = A.Buffer(0,j);
        F* Ak = A.Buffer(0,k);
        for( Int i=0; i<height; ++i )
        {
            F temp = Aj[i];
            Aj[i] = Ak[i];
            Ak[i] = temp;
        }
    }
}

template<typename F,Dist U1,Dist V1,
                    Dist U2,Dist V2>
inline void
ApplyColumnPivots( DistMatrix<F,U1,V1>& A, const DistMatrix<Int,U2,V2>& p )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyColumnPivots"))
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyColumnPivots( A, p_STAR_STAR );
}

template<typename F,Dist U1,Dist V1,
                    Dist U2,Dist V2>
inline void
ApplyInverseColumnPivots
( DistMatrix<F,U1,V1>& A, const DistMatrix<Int,U2,V2>& p )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyInverseColumnPivots"))
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyInverseColumnPivots( A, p_STAR_STAR );
}

template<typename F,Dist U,Dist V>
inline void
ApplyColumnPivots( DistMatrix<F,U,V>& A, const DistMatrix<Int,STAR,STAR>& p )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyColumnPivots"))
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyColumnPivots( A, image, preimage );
}

template<typename F,Dist U,Dist V>
inline void
ApplyInverseColumnPivots
( DistMatrix<F,U,V>& A, const DistMatrix<Int,STAR,STAR>& p )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyInverseColumnPivots"))
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyColumnPivots( A, preimage, image );
}

template<typename F> 
inline void
ApplyColumnPivots
( Matrix<F>& A, 
  const std::vector<Int>& image,
  const std::vector<Int>& preimage )
{
    const Int b = image.size();
    DEBUG_ONLY(
        CallStackEntry cse("ApplyColumnPivots");
        if( A.Width() < b || b != int(preimage.size()) )
            LogicError
            ("image and preimage must be vectors of equal length that are not "
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
    Matrix<F> APreimageCopy( m, b );
    for( Int j=0; j<b; ++j )
    {
        const Int jPre = preimage[j];
        if( jPre >= b )
            MemCopy( APreimageCopy.Buffer(0,j), A.LockedBuffer(0,jPre), m );
    }

    // Apply the permutations
    for( Int j=0; j<b; ++j )
    {
        const Int jPre = preimage[j];
        const Int jPost = image[j];
        // Move row[i] into row[image[i]]
        MemCopy( A.Buffer(0,jPost), AColPanCopy.LockedBuffer(0,j), m );
        // Move row[preimage[i]] into row[i]
        if( jPre >= b )
            MemCopy( A.Buffer(0,j), APreimageCopy.LockedBuffer(0,j), m );
    }
}

template<typename F,Dist U,Dist V> 
inline void
ApplyColumnPivots( DistMatrix<F,U,V>& A, const PivotMeta& oldMeta )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyColumnPivots");
        if( A.RowComm() != oldMeta.comm )
            LogicError("Invalid communicator in metadata");
        if( A.RowAlign() != oldMeta.align )
            LogicError("Invalid alignment in metadata");
    )
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;

    const Int localHeight = A.LocalHeight();
    PivotMeta meta = oldMeta;
    meta.ScaleUp( localHeight );

    // Fill vectors with the send data
    auto offsets = meta.sendDispls;
    const int totalSend = meta.TotalSend();
    std::vector<F> sendData( mpi::Pad(totalSend) );
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
    std::vector<F> recvData( mpi::Pad(totalRecv) );
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

template<typename F,Dist U,Dist V> 
inline void
ApplyColumnPivots
( DistMatrix<F,U,V>& A, 
  const std::vector<Int>& image, const std::vector<Int>& preimage )
{
    if( A.Participating() )
    {
        auto meta = FormPivotMeta( A.RowComm(), A.RowAlign(), image, preimage );
        ApplyColumnPivots( A, meta );
    }
}

} // namespace elem

#endif // ifndef ELEM_APPLYCOLUMNPIVOTS_HPP
