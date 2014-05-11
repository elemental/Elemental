/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_ONETWOONE_HPP
#define ELEM_ONETWOONE_HPP

#include ELEM_ZEROS_INC

namespace elem {

template<typename T> 
inline void
MakeOneTwoOne( Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOneTwoOne"))
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

template<typename T,Dist U,Dist V>
inline void
MakeOneTwoOne( DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOneTwoOne"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix 1-2-1");
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
                A.SetLocal( iLoc, jLoc, T(2) );
            else if( i == j-1 || i == j+1 )
                A.SetLocal( iLoc, jLoc, T(1) );
        }
    }
}

template<typename T,Dist U,Dist V>
inline void
MakeOneTwoOne( BlockDistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOneTwoOne"))
    if( A.Height() != A.Width() )
        LogicError("Cannot make a non-square matrix 1-2-1");
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
    DEBUG_ONLY(CallStackEntry cse("OneTwoOne"))
    A.Resize( n, n );
    MakeOneTwoOne( A );
}

template<typename T,Dist U,Dist V> 
inline void
OneTwoOne( DistMatrix<T,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("OneTwoOne"))
    A.Resize( n, n );
    MakeOneTwoOne( A );
}

template<typename T,Dist U,Dist V> 
inline void
OneTwoOne( BlockDistMatrix<T,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("OneTwoOne"))
    A.Resize( n, n );
    MakeOneTwoOne( A );
}

} // namespace elem

#endif // ifndef ELEM_ONETWOONE_HPP
