/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PARTER_HPP
#define EL_PARTER_HPP

namespace El {

template<typename F> 
inline void
Parter( Matrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    const F oneHalf = F(1)/F(2);
    P.Resize( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            P.Set( i, j, F(1)/(F(i)-F(j)+oneHalf) );
}

template<typename F>
inline void
Parter( AbstractDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    const F oneHalf = F(1)/F(2);
    P.Resize( n, n );
    const Int localHeight = P.LocalHeight();
    const Int localWidth = P.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = P.GlobalRow(iLoc);
            P.SetLocal( iLoc, jLoc, F(1)/(F(i)-F(j)+oneHalf) );
        }
    }
}

template<typename F>
inline void
Parter( AbstractBlockDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    const F oneHalf = F(1)/F(2);
    P.Resize( n, n );
    const Int localHeight = P.LocalHeight();
    const Int localWidth = P.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = P.GlobalRow(iLoc);
            P.SetLocal( iLoc, jLoc, F(1)/(F(i)-F(j)+oneHalf) );
        }
    }
}

} // namespace El

#endif // ifndef EL_PARTER_HPP
