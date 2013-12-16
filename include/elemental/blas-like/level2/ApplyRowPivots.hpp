/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_APPLYROWPIVOTS_HPP
#define ELEM_LAPACK_APPLYROWPIVOTS_HPP

#include "elemental/blas-like/level2/ComposePivots.hpp"

namespace elem {

template<typename F>
inline void
ApplyRowPivots( Matrix<F>& A, const Matrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyRowPivots");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( p.Height() > A.Height() )
            LogicError("p cannot be larger than height of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = p.Height();
    const Int ldim = A.LDim();
    for( Int i=0; i<numPivots; ++i )
    {
        const Int k = p.Get(i,0);
        F* Ai = A.Buffer(i,0);
        F* Ak = A.Buffer(k,0);
        for( Int j=0; j<width; ++j )
        {
            F temp = Ai[j*ldim];
            Ai[j*ldim] = Ak[j*ldim];
            Ak[j*ldim] = temp;
        }
    }
}

template<typename F>
inline void
ApplyInverseRowPivots( Matrix<F>& A, const Matrix<Int>& p )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyInverseRowPivots");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( p.Height() > A.Height() )
            LogicError("p cannot be larger than height of A");
    )
    const Int height = A.Height();
    const Int width = A.Width();
    if( height == 0 || width == 0 )
        return;

    const Int numPivots = p.Height();
    const Int ldim = A.LDim();
    for( Int i=numPivots-1; i>=0; --i )
    {
        const Int k = p.Get(i,0);
        F* Ai = A.Buffer(i,0);
        F* Ak = A.Buffer(k,0);
        for( Int j=0; j<width; ++j )
        {
            F temp = Ai[j*ldim];
            Ai[j*ldim] = Ak[j*ldim];
            Ak[j*ldim] = temp;
        }
    }
}

template<typename F,Distribution U1,Distribution V1,
                    Distribution U2,Distribution V2>
inline void
ApplyRowPivots( DistMatrix<F,U1,V1>& A, const DistMatrix<Int,U2,V2>& p )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyRowPivots"))
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyRowPivots( A, p_STAR_STAR );
}

template<typename F,Distribution U1,Distribution V1,
                    Distribution U2,Distribution V2>
inline void
ApplyInverseRowPivots
( DistMatrix<F,U1,V1>& A, const DistMatrix<Int,U2,V2>& p )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyInverseRowPivots"))
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ApplyInverseRowPivots( A, p_STAR_STAR );
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyRowPivots( DistMatrix<F,U,V>& A, const DistMatrix<Int,STAR,STAR>& p )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyRowPivots"))
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyRowPivots( A, image, preimage );
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyInverseRowPivots
( DistMatrix<F,U,V>& A, const DistMatrix<Int,STAR,STAR>& p )
{
    DEBUG_ONLY(CallStackEntry cse("ApplyInverseRowPivots"))
    std::vector<Int> image, preimage;
    ComposePivots( p, image, preimage );
    ApplyRowPivots( A, preimage, image );
}

template<typename F> 
inline void
ApplyRowPivots
( Matrix<F>& A, 
  const std::vector<Int>& image,
  const std::vector<Int>& preimage )
{
    const Int b = image.size();
    DEBUG_ONLY(
        CallStackEntry cse("ApplyRowPivots");
        if( A.Height() < b || b != (int)preimage.size() )
            LogicError
            ("image and preimage must be vectors of equal length that are not "
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
    Matrix<F> APreimageCopy( b, n );
    for( Int i=0; i<b; ++i ) 
    {
        const Int iPre = preimage[i];
        if( iPre >= b )
            for( Int j=0; j<n; ++j )
                APreimageCopy.Set(i,j,A.Get(iPre,j));
    }

    // Apply the permutations
    for( Int i=0; i<b; ++i )
    {
        const Int iPre = preimage[i];
        const Int iPost = image[i];
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

template<typename F,Distribution U,Distribution V>
inline void
ApplyRowPivots( DistMatrix<F,U,V>& A, const PivotMeta& oldMeta )
{
    DEBUG_ONLY(
        CallStackEntry cse("ApplyRowPivots");
        if( A.ColComm() != oldMeta.comm )
            LogicError("Invalid communicator in metadata");
        if( A.ColAlign() != oldMeta.align )
            LogicError("Invalid alignment in metadata");
    )
    if( A.Height() == 0 || A.Width() == 0 || !A.Participating() )
        return;

    const Int localWidth = A.LocalWidth();
    const Int ldim = A.LDim();
    PivotMeta meta = oldMeta;
    meta.ScaleUp( localWidth );

    // Fill vectors with the send data
    auto offsets = meta.sendDispls;
    const int totalSend = meta.TotalSend();
    std::vector<F> sendData( mpi::Pad(totalSend) );
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
        const int iLoc = meta.recvIdx[recv];
        const int rank = meta.recvRanks[recv];
        StridedMemCopy
        ( A.Buffer(iLoc,0),         ldim,
          &recvData[offsets[rank]], 1,    localWidth );
        offsets[rank] += localWidth;
    }
}

template<typename F,Distribution U,Distribution V>
inline void
ApplyRowPivots
( DistMatrix<F,U,V>& A,
  const std::vector<Int>& image, const std::vector<Int>& preimage )
{
    if( A.Participating() )
    {
        auto meta = FormPivotMeta( A.ColComm(), A.ColAlign(), image, preimage );
        ApplyRowPivots( A, meta );
    }
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_APPLYROWPIVOTS_HPP
