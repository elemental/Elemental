/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_COMPOSEPIVOTS_HPP
#define LAPACK_COMPOSEPIVOTS_HPP

namespace elem {

// Meant for composing an entire pivot vector for an n x n matrix.
// Requires O(n) work for an n x n matrix.

inline void
ComposePivots
( const Matrix<int>& p,
  std::vector<int>& image, std::vector<int>& preimage )
{
#ifndef RELEASE
    PushCallStack("ComposePivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
#endif
    const int n = p.Height();
    const int* pBuffer = p.LockedBuffer();

    // Construct the preimage of {0,...,n-1} under the permutation in O(n) work
    preimage.resize( n );
    for( int i=0; i<n; ++i )
        preimage[i] = i;
    for( int i=0; i<n; ++i )
    {
        const int j = pBuffer[i];
        const int k = preimage[j];
        preimage[j] = preimage[i];
        preimage[i] = k;
    }

    // Construct the image of {0,...,n-1} under the permutation in O(n) work
    image.resize( n );
    for( int i=0; i<n; ++i )
        image[preimage[i]] = i;
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
ComposePivots
( const DistMatrix<int,VC,STAR>& p,
  std::vector<int>& image, std::vector<int>& preimage )
{
#ifndef RELEASE    
    PushCallStack("ComposePivots");
#endif
    DistMatrix<int,STAR,STAR> p_STAR_STAR( p );
    ComposePivots( p_STAR_STAR, image, preimage );
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
ComposePivots
( const DistMatrix<int,STAR,STAR>& p, 
  std::vector<int>& image, std::vector<int>& preimage )
{ ComposePivots( p.LockedMatrix(), image, preimage ); }

// Meant for composing the pivots from a panel factorization, where b pivots
// were performed in an n x n matrix. 
// Requires O(b^2) work.

inline void
ComposePivots
( const Matrix<int>& p, int pivotOffset,
  std::vector<int>& image, std::vector<int>& preimage )
{
#ifndef RELEASE
    PushCallStack("ComposePivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( pivotOffset < 0 )
        throw std::logic_error("pivotOffset must be non-negative");
#endif
    const int b = p.Height();
    const int* pBuffer = p.LockedBuffer();

    // Construct the preimage of {0,...,b-1} under permutation in O(b^2) work
    preimage.resize( b );
    for( int i=0; i<b; ++i )
    {
        int k = pBuffer[i]-pivotOffset;
        for( int j=i-1; j>=0; --j )
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
    for( int i=0; i<b; ++i )
    {
        int k = i;
        for( int j=0; j<std::min(k+1,b); ++j )
        {
            if( pBuffer[j]-pivotOffset == k )
                k = j;
            else if( j == k )
                k = pBuffer[j]-pivotOffset;
        }
        image[i] = k;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
ComposePivots
( const DistMatrix<int,STAR,STAR>& p, int pivotOffset,
  std::vector<int>& image, std::vector<int>& preimage )
{ ComposePivots( p.LockedMatrix(), pivotOffset, image, preimage ); }

} // namespace elem

#endif // ifndef LAPACK_COMPOSEPIVOTS_HPP
