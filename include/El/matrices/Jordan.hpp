/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_JORDAN_HPP
#define EL_JORDAN_HPP

namespace El {

template<typename T> 
inline void
MakeJordan( Matrix<T>& J, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("MakeJordan"))
    Zero( J );
    const Int m = J.Height();
    const Int n = J.Width();
    for( Int j=0; j<std::min(m,n); ++j )
    {
        J.Set( j, j, lambda );
        if( j != 0 )
            J.Set( j-1, j, T(1) );
    }
}

template<typename T>
inline void
MakeJordan( AbstractDistMatrix<T>& J, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("MakeJordan"))
    Zero( J.Matrix() );

    const Int localHeight = J.LocalHeight();
    const Int localWidth = J.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = J.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i == j )
                J.SetLocal( iLoc, jLoc, lambda );
            else if( i == j-1 )
                J.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

template<typename T>
inline void
MakeJordan( AbstractBlockDistMatrix<T>& J, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("MakeJordan"))
    Zero( J.Matrix() );

    const Int localHeight = J.LocalHeight();
    const Int localWidth = J.LocalWidth();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = J.GlobalCol(jLoc);
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = J.GlobalRow(iLoc);
            if( i == j )
                J.SetLocal( iLoc, jLoc, lambda );
            else if( i == j-1 )
                J.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

template<typename T>
inline void
Jordan( Matrix<T>& J, Int n, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Jordan"))
    J.Resize( n, n );
    MakeJordan( J, lambda );
}

template<typename T>
inline void
Jordan( AbstractDistMatrix<T>& J, Int n, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Jordan"))
    J.Resize( n, n );
    MakeJordan( J, lambda );
}

template<typename T>
inline void
Jordan( AbstractBlockDistMatrix<T>& J, Int n, T lambda )
{
    DEBUG_ONLY(CallStackEntry cse("Jordan"))
    J.Resize( n, n );
    MakeJordan( J, lambda );
}

} // namespace El

#endif // ifndef EL_JORDAN_HPP
