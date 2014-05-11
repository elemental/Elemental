/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRIW_HPP
#define ELEM_TRIW_HPP

#include ELEM_SETDIAGONAL_INC
#include ELEM_ZEROS_INC

namespace elem {

template<typename T> 
inline void
TriW( Matrix<T>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals must be non-negative");
    Zeros( A, n, n );
    SetDiagonal( A, 1 );
    for( Int j=0; j<Min(n-1,k); ++j ) 
        SetDiagonal( A, alpha, j+1 );
}

template<typename T,Dist U,Dist V>
inline void
TriW( DistMatrix<T,U,V>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals must be non-negative");
    Zeros( A, n, n );
    SetDiagonal( A, 1 );
    for( Int j=0; j<Min(n-1,k); ++j )
        SetDiagonal( A, alpha, j+1 );
}

template<typename T,Dist U,Dist V>
inline void
TriW( BlockDistMatrix<T,U,V>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals must be non-negative");
    Zeros( A, n, n );
    SetDiagonal( A, 1 );
    for( Int j=0; j<Min(n-1,k); ++j )
        SetDiagonal( A, alpha, j+1 );
}

} // namespace elem

#endif // ifndef ELEM_TRIW_HPP
