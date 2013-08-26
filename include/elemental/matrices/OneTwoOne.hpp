/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_ONETWOONE_HPP
#define ELEM_MATRICES_ONETWOONE_HPP

#include "elemental/matrices/Zeros.hpp"

namespace elem {

template<typename T> 
inline void
MakeOneTwoOne( Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeOneTwoOne");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix 1-2-1");
    MakeZeros( A );

    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
    {
        A.Set( j, j, T(2) );
        if( j < n-1 )
        {
            A.Set( j+1, j, T(1) );
            A.Set( j, j+1, T(1) );
        }
    }
}

template<typename T,Distribution U,Distribution V>
inline void
MakeOneTwoOne( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeOneTwoOne");
#endif
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix 1-2-1");
    MakeZeros( A );

    const Int localHeight = A.LocalHeight();
    const Int localWidth = A.LocalWidth();
    const Int colShift = A.ColShift();
    const Int rowShift = A.RowShift();
    const Int colStride = A.ColStride();
    const Int rowStride = A.RowStride();
    for( Int jLoc=0; jLoc<localWidth; ++jLoc )
    {
        const Int j = rowShift + jLoc*rowStride;
        for( Int iLoc=0; iLoc<localHeight; ++iLoc )
        {
            const Int i = colShift + iLoc*colStride;
            if( i == j )
                A.SetLocal( iLoc, jLoc, T(2) );
            else if( i == j-1 || i == j+1 )
                A.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

template<typename T> 
inline void
OneTwoOne( Matrix<T>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("OneTwoOne");
#endif
    A.ResizeTo( n, n );
    MakeOneTwoOne( A );
}

template<typename T> 
inline Matrix<T>
OneTwoOne( Int n )
{
    Matrix<T> A( n, n );
    MakeOneTwoOne( A );
    return A;
}

template<typename T,Distribution U,Distribution V> 
inline void
OneTwoOne( DistMatrix<T,U,V>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("OneTwoOne");
#endif
    A.ResizeTo( n, n );
    MakeOneTwoOne( A );
}

template<typename T,Distribution U=MC,Distribution V=MR> 
inline DistMatrix<T,U,V>
OneTwoOne( const Grid& g, Int n )
{
    DistMatrix<T,U,V> A( n, n, g );
    MakeOneTwoOne( A );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_ONETWOONE_HPP
