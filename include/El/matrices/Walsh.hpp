/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_WALSH_HPP
#define EL_WALSH_HPP

namespace El {

// TODO: Get rid of MakeWalsh routine? It doesn't seem useful.

template<typename T> 
inline void
MakeWalsh( Matrix<T>& A, Int k, bool binary=false )
{
    DEBUG_ONLY(CallStackEntry cse("MakeWalsh"))
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

template<typename T>
inline void
MakeWalsh( AbstractDistMatrix<T>& A, Int k, bool binary=false )
{
    DEBUG_ONLY(CallStackEntry cse("MakeWalsh"))
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
    for( Unsigned jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Unsigned j = A.GlobalCol(jLoc);
        for( Unsigned iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Unsigned i = A.GlobalRow(iLoc);

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
    DEBUG_ONLY(CallStackEntry cse("Walsh"))
    if( k < 1 )
        LogicError("Walsh matrices are only defined for k>=1");
    const Unsigned n = 1u<<k;
    A.Resize( n, n );
    MakeWalsh( A, k, binary );
}

template<typename T>
inline void
Walsh( AbstractDistMatrix<T>& A, Int k, bool binary=false )
{
    DEBUG_ONLY(CallStackEntry cse("Walsh"))
    if( k < 1 )
        LogicError("Walsh matrices are only defined for k>=1");
    const Unsigned n = 1u<<k;
    A.Resize( n, n );
    MakeWalsh( A, k, binary );
}

} // namespace El

#endif // ifndef EL_WALSH_HPP
