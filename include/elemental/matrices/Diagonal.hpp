/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DIAGONAL_HPP
#define ELEM_DIAGONAL_HPP

#include "./Zeros.hpp"

namespace elem {

template<typename S,typename T> 
inline void
Diagonal( Matrix<S>& D, const std::vector<T>& d )
{
    DEBUG_ONLY(CallStackEntry cse("Diagonal"))
    const Int n = d.size();
    Zeros( D, n, n );

    for( Int j=0; j<n; ++j )
        D.Set( j, j, d[j] );
}

template<typename S,typename T,Dist U,Dist V>
inline void
Diagonal( DistMatrix<S,U,V>& D, const std::vector<T>& d )
{
    DEBUG_ONLY(CallStackEntry cse("Diagonal"))
    const Int n = d.size();
    Zeros( D, n, n );

    const Int localWidth = D.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = D.GlobalCol(jLoc);
        D.Set( j, j, d[j] );
    }
}

template<typename S,typename T,Dist U,Dist V>
inline void
Diagonal( BlockDistMatrix<S,U,V>& D, const std::vector<T>& d )
{
    DEBUG_ONLY(CallStackEntry cse("Diagonal"))
    const Int n = d.size();
    Zeros( D, n, n );

    const Int localWidth = D.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = D.GlobalCol(jLoc);
        D.Set( j, j, d[j] );
    }
}

} // namespace elem

#endif // ifndef ELEM_DIAGONAL_HPP
