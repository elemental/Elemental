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
#ifndef RELEASE
    CallStackEntry entry("ComposePivots");
    if( p.Width() != 1 )
        LogicError("p must be a column vector");
#endif
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
#ifndef RELEASE    
    CallStackEntry entry("ComposePivots");
#endif
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
#ifndef RELEASE
    CallStackEntry entry("ComposePivots");
    if( p.Width() != 1 )
        LogicError("p must be a column vector");
    if( pivotOffset < 0 )
        LogicError("pivotOffset must be non-negative");
#endif
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

} // namespace elem

#endif // ifndef ELEM_LAPACK_COMPOSEPIVOTS_HPP
