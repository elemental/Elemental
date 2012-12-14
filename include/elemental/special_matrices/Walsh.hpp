/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
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
                A.SetLocal( iLocal, jLocal, onValue );
            else
                A.SetLocal( iLocal, jLocal, offValue );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem
