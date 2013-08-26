/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_WALSH_HPP
#define ELEM_MATRICES_WALSH_HPP

namespace elem {

template<typename T> 
inline void
MakeWalsh( Matrix<T>& A, Int k, bool binary=false )
{
#ifndef RELEASE
    CallStackEntry cse("MakeWalsh");
#endif
    if( k < 1 )
        LogicError("Walsh matrices are only defined for k>=1");

    const Unsigned n = 1u<<k;
    if( A.Height() != Int(n) || A.Width() != Int(n) )
        LogicError("Invalid input matrix size");

    // Run a simple O(n^2 log n) algorithm for computing the entries
    // based upon successive sign flips
    const T onValue = 1;
    const T offValue = ( binary ? 0 : -1 );
    for( Unsigned j=0; j<n; ++j )
    {
        for( Unsigned i=0; i<n; ++i )
        {
            // Recurse on the quadtree, flipping the sign of the entry each
            // time we are in the bottom-right quadrant
            Unsigned r = i;     
            Unsigned s = j;
            Unsigned t = n;
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
}

template<typename T,Distribution U,Distribution V>
inline void
MakeWalsh( DistMatrix<T,U,V>& A, Int k, bool binary=false )
{
#ifndef RELEASE
    CallStackEntry cse("MakeWalsh");
#endif
    if( k < 1 )
        LogicError("Walsh matrices are only defined for k>=1");

    const Unsigned n = 1u<<k;
    if( A.Height() != Int(n) || A.Width() != Int(n) )
        LogicError("Invalid input matrix size");

    // Run an O(n^2 log n / p) algorithm based upon successive sign flips
    const T onValue = 1;
    const T offValue = ( binary ? 0 : -1 );
    const Unsigned localHeight = A.LocalHeight();
    const Unsigned localWidth = A.LocalWidth();
    const Unsigned colShift = A.ColShift();
    const Unsigned rowShift = A.RowShift();
    const Unsigned colStride = A.ColStride();
    const Unsigned rowStride = A.RowStride();
    for( Unsigned jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Unsigned j = rowShift + jLoc*rowStride;
        for( Unsigned iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Unsigned i = colShift + iLoc*colStride;

            // Recurse on the quadtree, flipping the sign of the entry each
            // time we are in the bottom-right quadrant
            Unsigned r = i;     
            Unsigned s = j;
            Unsigned t = n;
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
                A.SetLocal( iLoc, jLoc, onValue );
            else
                A.SetLocal( iLoc, jLoc, offValue );
        }
    }
}

template<typename T> 
inline void
Walsh( Matrix<T>& A, Int k, bool binary=false )
{
#ifndef RELEASE
    CallStackEntry cse("Walsh");
#endif
    if( k < 1 )
        LogicError("Walsh matrices are only defined for k>=1");
    const Unsigned n = 1u<<k;
    A.ResizeTo( n, n );
    MakeWalsh( A, k, binary );
}

template<typename T> 
inline Matrix<T>
Walsh( Int k, bool binary=false )
{
    Matrix<T> A;
    Walsh( A, k, binary );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void
Walsh( DistMatrix<T,U,V>& A, Int k, bool binary=false )
{
#ifndef RELEASE
    CallStackEntry cse("Walsh");
#endif
    if( k < 1 )
        LogicError("Walsh matrices are only defined for k>=1");
    const Unsigned n = 1u<<k;
    A.ResizeTo( n, n );
    MakeWalsh( A, k, binary );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Walsh( const Grid& g, Int k, bool binary=false )
{
    DistMatrix<T,U,V> A(g);
    Walsh( A, k, binary ); 
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_WALSH_HPP
