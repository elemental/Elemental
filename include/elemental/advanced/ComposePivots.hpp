/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elemental {

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

    // Construct the image of {0,...,n-1} under the permutation in O(n) work
    image.resize( n );
    for( int i=0; i<n; ++i )
        image[i] = i;
    for( int i=0; i<n; ++i )
    {
        const int j = pBuffer[i];
        const int k = image[j];
        image[j] = image[i];
        image[i] = k;
    }

    // Construct the preimage of {0,...,n-1} under the permutation in O(n) work
    preimage.resize( n );
    for( int i=0; i<n; ++i )
        preimage[image[i]] = i;
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
{ ComposePivots( p.LockedLocalMatrix(), image, preimage ); }

// Meant for composing the pivots from a panel factorization, where b pivots
// were performed in an n x n matrix. 
// Requires O(b^2) work.

inline void
internal::ComposePanelPivots
( const Matrix<int>& p, int pivotOffset,
  std::vector<int>& image, std::vector<int>& preimage )
{
#ifndef RELEASE
    PushCallStack("internal::ComposePanelPivots");
    if( p.Width() != 1 )
        throw std::logic_error("p must be a column vector");
    if( pivotOffset < 0 )
        throw std::logic_error("pivotOffset must be non-negative");
#endif
    const int b = p.Height();
    const int* pBuffer = p.LockedBuffer();

    // Construct the image of {0,...,b-1} under the permutation in O(b^2) work
    image.resize( b );
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
        image[i] = k;
    }
    
    // Construct the preimage of {0,...,b-1} under the permutation in 
    // O(b^2) work
    preimage.resize( b );
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
        preimage[i] = k;
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

inline void
internal::ComposePanelPivots
( const DistMatrix<int,STAR,STAR>& p, int pivotOffset,
  std::vector<int>& image, std::vector<int>& preimage )
{ ComposePanelPivots( p.LockedLocalMatrix(), pivotOffset, image, preimage ); }

} // namespace elemental
