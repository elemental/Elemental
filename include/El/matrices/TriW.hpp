/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_TRIW_HPP
#define EL_TRIW_HPP

#include EL_ZEROS_INC

namespace El {

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

template<typename T>
inline void
TriW( AbstractDistMatrix<T>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals must be non-negative");
    Zeros( A, n, n );
    SetDiagonal( A, 1 );
    for( Int j=0; j<Min(n-1,k); ++j )
        SetDiagonal( A, alpha, j+1 );
}

template<typename T>
inline void
TriW( AbstractBlockDistMatrix<T>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals must be non-negative");
    Zeros( A, n, n );
    SetDiagonal( A, 1 );
    for( Int j=0; j<Min(n-1,k); ++j )
        SetDiagonal( A, alpha, j+1 );
}

} // namespace El

#endif // ifndef EL_TRIW_HPP
