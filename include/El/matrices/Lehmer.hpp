/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LEHMER_HPP
#define EL_LEHMER_HPP

namespace El {

template<typename F> 
inline void
Lehmer( Matrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Lehmer"))
    L.Resize( n, n );
    for( Int j=0; j<n; ++j )
    {
        for( Int i=0; i<j; ++i )
            L.Set( i, j, F(i+1)/F(j+1) );
        for( Int i=j; i<n; ++i )
            L.Set( i, j, F(j+1)/F(i+1) );
    }
}

template<typename F>
inline void
Lehmer( AbstractDistMatrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Lehmer"))
    L.Resize( n, n );
    const Int localHeight = L.LocalHeight();
    const Int localWidth = L.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = L.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = L.GlobalRow(iLoc);
            if( i < j )
                L.SetLocal( iLoc, jLoc, F(i+1)/F(j+1) );
            else
                L.SetLocal( iLoc, jLoc, F(j+1)/F(i+1) );
        }
    }
}

template<typename F>
inline void
Lehmer( AbstractBlockDistMatrix<F>& L, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Lehmer"))
    L.Resize( n, n );
    const Int localHeight = L.LocalHeight();
    const Int localWidth = L.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = L.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = L.GlobalRow(iLoc);
            if( i < j )
                L.SetLocal( iLoc, jLoc, F(i+1)/F(j+1) );
            else
                L.SetLocal( iLoc, jLoc, F(j+1)/F(i+1) );
        }
    }
}

} // namespace El

#endif // ifndef EL_LEHMER_HPP
