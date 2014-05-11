/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_PEI_HPP
#define ELEM_PEI_HPP

namespace elem {

template<typename T> 
inline void
Pei( Matrix<T>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    P.Resize( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            P.Set( i, j, T(1) );
    for( Int j=0; j<n; ++j )
        P.Update( j, j, alpha );
}

template<typename T,Dist U,Dist V>
inline void
Pei( DistMatrix<T,U,V>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    P.Resize( n, n );
    const Int localHeight = P.LocalHeight();
    const Int localWidth = P.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = P.GlobalRow(iLoc);
            P.SetLocal( iLoc, jLoc, T(1) );
            if( i == j )
                P.UpdateLocal( iLoc, jLoc, alpha );
        }
    }
}

template<typename T,Dist U,Dist V>
inline void
Pei( BlockDistMatrix<T,U,V>& P, Int n, T alpha )
{
    DEBUG_ONLY(CallStackEntry cse("Pei"))
    P.Resize( n, n );
    const Int localHeight = P.LocalHeight();
    const Int localWidth = P.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = P.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = P.GlobalRow(iLoc);
            P.SetLocal( iLoc, jLoc, T(1) );
            if( i == j )
                P.UpdateLocal( iLoc, jLoc, alpha );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_PEI_HPP
