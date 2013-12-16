/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_COMPOSEPIVOTS_HPP
#define ELEM_LAPACK_COMPOSEPIVOTS_HPP

#include <algorithm>

namespace elem {

// Meant for composing an entire pivot vector for an n x n matrix.
// Requires O(n) work for an n x n matrix.

inline void
ComposePivots
( const Matrix<Int>& p,
  std::vector<Int>& image, std::vector<Int>& preimage )
{
    DEBUG_ONLY(
        CallStackEntry cse("ComposePivots");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
    )
    const Int n = p.Height();
    if( n == 0 )
    {
        image.resize( 0 );
        preimage.resize( 0 );
        return;
    }

    const Int* pBuffer = p.LockedBuffer();
    const Int range = *std::max_element( pBuffer, pBuffer+n ) + 1;

    // Construct the preimage of {0,...,range-1} under the permutation in 
    // O(range) work
    preimage.resize( range );
    for( Int i=0; i<range; ++i )
        preimage[i] = i;
    for( Int i=0; i<n; ++i )
    {
        const Int j = pBuffer[i];
        std::swap( preimage[i], preimage[j] );
    }

    // Construct the image of {0,...,range-1} under the permutation in 
    // O(range) work
    image.resize( range );
    for( Int i=0; i<range; ++i )
        image[preimage[i]] = i;
}

inline void
ComposePivots
( const DistMatrix<Int,VC,STAR>& p,
  std::vector<Int>& image, std::vector<Int>& preimage )
{
    DEBUG_ONLY(CallStackEntry cse("ComposePivots"))
    DistMatrix<Int,STAR,STAR> p_STAR_STAR( p );
    ComposePivots( p_STAR_STAR, image, preimage );
}

inline void
ComposePivots
( const DistMatrix<Int,STAR,STAR>& p, 
  std::vector<Int>& image, std::vector<Int>& preimage )
{ ComposePivots( p.LockedMatrix(), image, preimage ); }

// Meant for composing the pivots from a panel factorization, where b pivots
// were performed in an n x n matrix. 
// Requires O(b^2) work.

inline void
ComposePivots
( const Matrix<Int>& p, Int pivotOffset,
  std::vector<Int>& image, std::vector<Int>& preimage )
{
    DEBUG_ONLY(
        CallStackEntry cse("ComposePivots");
        if( p.Width() != 1 )
            LogicError("p must be a column vector");
        if( pivotOffset < 0 )
            LogicError("pivotOffset must be non-negative");
    )
    const Int b = p.Height();
    const Int* pBuffer = p.LockedBuffer();

    // Construct the preimage of {0,...,b-1} under permutation in O(b^2) work
    preimage.resize( b );
    for( Int i=0; i<b; ++i )
    {
        Int k = pBuffer[i]-pivotOffset;
        for( Int j=i-1; j>=0; --j )
        {
            if( pBuffer[j]-pivotOffset == k )
                k = j;
            else if( j == k )
                k = pBuffer[j]-pivotOffset;
        }
        preimage[i] = k;
    }
    
    // Construct the image of {0,...,b-1} under the permutation in O(b^2) work
    image.resize( b );
    for( Int i=0; i<b; ++i )
    {
        Int k = i;
        for( Int j=0; j<Min(k+1,b); ++j )
        {
            if( pBuffer[j]-pivotOffset == k )
                k = j;
            else if( j == k )
                k = pBuffer[j]-pivotOffset;
        }
        image[i] = k;
    }
}

inline void
ComposePivots
( const DistMatrix<Int,STAR,STAR>& p, Int pivotOffset,
  std::vector<Int>& image, std::vector<Int>& preimage )
{ ComposePivots( p.LockedMatrix(), pivotOffset, image, preimage ); }

inline PivotMeta
FormPivotMeta
( mpi::Comm comm, Int align,
  const std::vector<Int>& image,
  const std::vector<Int>& preimage )
{
    DEBUG_ONLY(CallStackEntry cse("FormPivotMeta"))
    PivotMeta meta;
    meta.comm = comm;
    meta.align = align;

    const Int b = image.size();
    const Int rank = mpi::CommRank( comm );
    const Int stride = mpi::CommSize( comm );
    const Int shift = Shift( rank, align, stride );

    // Form the metadata
    //
    // Extract the send and recv counts from the image and preimage.
    // There are three different types of exchanges:
    //   1. [0,b) -> [0,b)
    //   2. [0,b) -> [b,n)
    //   3. [b,n) -> [0,b)
    // The fourth possibility, [b,n) -> [b,n), is impossible due to the 
    // fact that indices pulled in from [b,n) are stuck in [0,b) due to the
    // fact that the i'th pivot exchanges index i with some index k >= i.
    // 
    meta.sendCounts.resize( stride, 0 );
    meta.recvCounts.resize( stride, 0 );

    for( Int j=0; j<b; ++j )
    {
        // Handle send 
        if( rank == ((align+j)%stride) )
        {
            const Int jLoc = (j-shift) / stride;
            const Int sendTo = (align+image[j]) % stride;
            meta.sendIdx.push_back( jLoc );
            meta.sendRanks.push_back( sendTo );
            ++meta.sendCounts[sendTo];
        }
        if( preimage[j] >= b && rank == ((align+preimage[j])%stride) )
        {
            const Int jLoc = (preimage[j]-shift) / stride;
            const Int sendTo = (align+j) % stride;
            meta.sendIdx.push_back( jLoc );
            meta.sendRanks.push_back( sendTo );
            ++meta.sendCounts[sendTo];
        }
        // Handle recv
        if( rank == ((align+image[j])%stride) )
        {
            const Int jLoc = (image[j]-shift) / stride;     
            const Int recvFrom = (align+j) % stride;
            meta.recvIdx.push_back( jLoc );
            meta.recvRanks.push_back( recvFrom );
            ++meta.recvCounts[recvFrom];
        }
        if( preimage[j] >= b && rank == ((align+j)%stride) )
        {
            const Int jLoc = (j-shift) / stride;
            const Int recvFrom = (align+preimage[j]) % stride;
            meta.recvIdx.push_back( jLoc );
            meta.recvRanks.push_back( recvFrom );
            ++meta.recvCounts[recvFrom];
        }
    }

    // Construct the send and recv displacements from the counts
    meta.sendDispls.resize( stride );
    meta.recvDispls.resize( stride );
    Int totalSend=0, totalRecv=0;
    for( Int i=0; i<stride; ++i )
    {
        meta.sendDispls[i] = totalSend;
        meta.recvDispls[i] = totalRecv;
        totalSend += meta.sendCounts[i];
        totalRecv += meta.recvCounts[i];
    }
    DEBUG_ONLY(
        if( totalSend != totalRecv )
            LogicError
            ("Send and recv counts do not match: send=",totalSend,", recv=",
             totalRecv);
    )

    return meta;
}

} // namespace elem

#endif // ifndef ELEM_LAPACK_COMPOSEPIVOTS_HPP
