/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MINIJ_HPP
#define EL_MINIJ_HPP

namespace El {

template<typename T> 
inline void
MinIJ( Matrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            M.Set( i, j, std::min(i+1,j+1) );
}

template<typename T>
inline void
MinIJ( AbstractDistMatrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    const Int localHeight = M.LocalHeight();
    const Int localWidth = M.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = M.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = M.GlobalRow(iLoc);
            M.SetLocal( iLoc, jLoc, std::min(i+1,j+1) );
        }
    }
}

template<typename T>
inline void
MinIJ( AbstractBlockDistMatrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    const Int localHeight = M.LocalHeight();
    const Int localWidth = M.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = M.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = M.GlobalRow(iLoc);
            M.SetLocal( iLoc, jLoc, std::min(i+1,j+1) );
        }
    }
}

} // namespace El

#endif // ifndef EL_MINIJ_HPP
