/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_JORDAN_HPP
#define ELEM_MATRICES_JORDAN_HPP

#include "elemental/blas-like/level1/Zero.hpp"

namespace elem {

template<typename T> 
inline void
MakeJordan( Matrix<T>& J, T lambda )
{
#ifndef RELEASE
    CallStackEntry cse("MakeJordan");
#endif
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

template<typename T,Distribution U,Distribution V>
inline void
MakeJordan( DistMatrix<T,U,V>& J, T lambda )
{
#ifndef RELEASE
    CallStackEntry cse("MakeJordan");
#endif
    Zero( J.Matrix() );

    const Int localHeight = J.LocalHeight();
    const Int localWidth = J.LocalWidth();
    const Int colShift = J.ColShift();
    const Int rowShift = J.RowShift();
    const Int colStride = J.ColStride();
    const Int rowStride = J.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
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
#ifndef RELEASE
    CallStackEntry cse("Jordan");
#endif
    J.ResizeTo( n, n );
    MakeJordan( J, lambda );
}

template<typename T>
inline Matrix<T>
Jordan( Int n, T lambda )
{
    Matrix<T> J( n, n );
    MakeJordan( J, lambda );
    return J;
}

template<typename T,Distribution U,Distribution V>
inline void
Jordan( DistMatrix<T,U,V>& J, Int n, T lambda )
{
#ifndef RELEASE
    CallStackEntry cse("Jordan");
#endif
    J.ResizeTo( n, n );
    MakeJordan( J, lambda );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Jordan( const Grid& g, Int n, T lambda )
{
    DistMatrix<T,U,V> J( n, n, g );
    MakeJordan( J, lambda );
    return J;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_JORDAN_HPP
