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

namespace elem {

template<typename T> 
inline void
Walsh( int k, Matrix<T>& A, bool binary )
{
#ifndef RELEASE
    PushCallStack("Walsh");
#endif
    if( k < 1 )
        throw std::logic_error("Walsh matrices are only defined for k>=1");

    const unsigned n = 1u<<k;
    A.ResizeTo( n, n );

    // Run a simple O(n^2 log n) algorithm for computing the entries
    // based upon successive sign flips
    const T onValue = 1;
    const T offValue = ( binary ? 0 : -1 );
    for( unsigned j=0; j<n; ++j )
    {
        for( unsigned i=0; i<n; ++i )
        {
            // Recurse on the quadtree, flipping the sign of the entry each
            // time we are in the bottom-right quadrant
            unsigned r = i;     
            unsigned s = j;
            unsigned t = n;
            bool on = true;
            while( t != 1u )
            {
                t >>= 1;
                if( r >= t && s >= t )
                    on = !on;
                r %= t;
                s %= t;
            }

            if( on )
                A.Set( i, j, onValue );
            else
                A.Set( i, j, offValue );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T,Distribution U,Distribution V>
inline void
Walsh( int k, DistMatrix<T,U,V>& A, bool binary )
{
#ifndef RELEASE
    PushCallStack("Walsh");
#endif
    if( k < 1 )
        throw std::logic_error("Walsh matrices are only defined for k>=1");

    const unsigned n = 1u<<k;
    A.ResizeTo( n, n );

    // Run an O(n^2 log n / p) algorithm based upon successive sign flips
    const T onValue = 1;
    const T offValue = ( binary ? 0 : -1 );
    const unsigned localHeight = A.LocalHeight();
    const unsigned localWidth = A.LocalWidth();
    const unsigned colShift = A.ColShift();
    const unsigned rowShift = A.RowShift();
    const unsigned colStride = A.ColStride();
    const unsigned rowStride = A.RowStride();
    for( unsigned jLocal=0; jLocal<localWidth; ++jLocal )
    {
        const unsigned j = rowShift + jLocal*rowStride;
        for( unsigned iLocal=0; iLocal<localHeight; ++iLocal )
        {
            const unsigned i = colShift + iLocal*colStride;

            // Recurse on the quadtree, flipping the sign of the entry each
            // time we are in the bottom-right quadrant
            unsigned r = i;     
            unsigned s = j;
            unsigned t = n;
            bool on = true;
            while( t != 1u )
            {
                t >>= 1;
                if( r >= t && s >= t )
                    on = !on;
                r %= t;
                s %= t;
            }
            if( on )
                A.SetLocalEntry( iLocal, jLocal, onValue );
            else
                A.SetLocalEntry( iLocal, jLocal, offValue );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
