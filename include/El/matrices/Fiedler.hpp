/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_FIEDLER_HPP
#define EL_FIEDLER_HPP

namespace El {

template<typename F> 
inline void
Fiedler( Matrix<F>& A, const std::vector<F>& c )
{
    DEBUG_ONLY(CallStackEntry cse("Fiedler"))
    const Int n = c.size();
    A.Resize( n, n );
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<n; ++i )
            A.Set( i, j, Abs(c[i]-c[j]) );
}

template<typename F>
inline void
Fiedler( AbstractDistMatrix<F>& A, const std::vector<F>& c )
{
    DEBUG_ONLY(CallStackEntry cse("Fiedler"))
    const Int n = c.size();
    A.Resize( n, n );
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, Abs(c[i]-c[j]) );
        }
    }
}

template<typename F>
inline void
Fiedler( AbstractBlockDistMatrix<F>& A, const std::vector<F>& c )
{
    DEBUG_ONLY(CallStackEntry cse("Fiedler"))
    const Int n = c.size();
    A.Resize( n, n );
    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            A.SetLocal( iLoc, jLoc, Abs(c[i]-c[j]) );
        }
    }
}

} // namespace El

#endif // ifndef EL_FIEDLER_HPP
