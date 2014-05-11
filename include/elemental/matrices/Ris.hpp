/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_RIS_HPP
#define ELEM_RIS_HPP

namespace elem {

template<typename F> 
inline void
Ris( Matrix<F>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    const F oneHalf = F(1)/F(2);
    R.Resize( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            R.Set( i, j, oneHalf/(F(n-i-j)-oneHalf) );
}

template<typename F,Dist U,Dist V>
inline void
Ris( DistMatrix<F,U,V>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    const F oneHalf = F(1)/F(2);
    R.Resize( n, n );
    const Int localHeight = R.LocalHeight();
    const Int localWidth = R.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = R.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = R.GlobalRow(iLoc);
            R.SetLocal( iLoc, jLoc, oneHalf/(F(n-i-j)-oneHalf) );
        }
    }
}

template<typename F,Dist U,Dist V>
inline void
Ris( BlockDistMatrix<F,U,V>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    const F oneHalf = F(1)/F(2);
    R.Resize( n, n );
    const Int localHeight = R.LocalHeight();
    const Int localWidth = R.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = R.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = R.GlobalRow(iLoc);
            R.SetLocal( iLoc, jLoc, oneHalf/(F(n-i-j)-oneHalf) );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_RIS_HPP
