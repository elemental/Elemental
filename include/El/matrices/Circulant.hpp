/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CIRCULANT_HPP
#define EL_CIRCULANT_HPP

namespace El {

template<typename T> 
inline void
Circulant( Matrix<T>& A, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Circulant"))
    const Int n = a.size();
    A.Resize( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            A.Set( i, j, a[Mod(i-j,n)] );
}

template<typename T>
inline void
Circulant( AbstractDistMatrix<T>& A, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Circulant"))
    const Int n = a.size();
    A.Resize( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, a[Mod(i-j,n)] );
        }
    }
}

template<typename T>
inline void
Circulant( AbstractBlockDistMatrix<T>& A, const std::vector<T>& a )
{
    DEBUG_ONLY(CallStackEntry cse("Circulant"))
    const Int n = a.size();
    A.Resize( n, n );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, a[Mod(i-j,n)] );
        }
    }
}

} // namespace El

#endif // ifndef EL_CIRCULANT_HPP
