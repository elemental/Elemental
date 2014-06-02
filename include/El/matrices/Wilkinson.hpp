/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_WILKINSON_HPP
#define EL_WILKINSON_HPP

#include EL_ZEROS_INC

namespace El {

template<typename T> 
inline void
Wilkinson( Matrix<T>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("Wilkinson"))
    const Int n = 2*k+1;
    A.Resize( n, n );
    MakeZeros( A );

    for( Int j=0; j<n; ++j )
    {
        if( j <= k )
            A.Set( j, j, T(k-j) );
        else
            A.Set( j, j, T(j-k) );

        if( j > 0 )
            A.Set( j-1, j, T(1) );
        if( j < n-1 )
            A.Set( j+1, j, T(1) );
    }
}

template<typename T>
inline void
Wilkinson( AbstractDistMatrix<T>& A, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("Wilkinson"))
    const Int n = 2*k+1;
    A.Resize( n, n );
    MakeZeros( A );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = A.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = A.GlobalRow(iLoc);
            if( i == j )
            {
                if( j <= k )
                    A.SetLocal( iLoc, jLoc, T(k-j) );
                else
                    A.SetLocal( iLoc, jLoc, T(j-k) );
            }
            else if( i == j-1 || i == j+1 )
                A.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

} // namespace El

#endif // ifndef EL_WILKINSON_HPP
