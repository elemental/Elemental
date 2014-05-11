/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TOEPLITZ_HPP
#define ELEM_TOEPLITZ_HPP

namespace elem {

template<typename S,typename T> 
inline void
Toeplitz( Matrix<S>& A, Int m, Int n, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Toeplitz"))
    const Int length = m+n-1;
    if( a.size() != Unsigned(length) )
        LogicError("a was the wrong size");
    A.Resize( m, n );

    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, a[i-j+(n-1)] );
}

template<typename S,typename T,Dist U,Dist V>
inline void
Toeplitz( DistMatrix<S,U,V>& A, Int m, Int n, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Toeplitz"))
    const Int length = m+n-1;
    if( a.size() != Unsigned(length) )
        LogicError("a was the wrong size");
    A.Resize( m, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, a[i-j+(n-1)] );
        }
    }
}

template<typename S,typename T,Dist U,Dist V>
inline void
Toeplitz( BlockDistMatrix<S,U,V>& A, Int m, Int n, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Toeplitz"))
    const Int length = m+n-1;
    if( a.size() != Unsigned(length) )
        LogicError("a was the wrong size");
    A.Resize( m, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, a[i-j+(n-1)] );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_TOEPLITZ_HPP
