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

#include ELEM_TOEPLITZ_INC

namespace elem {

template<typename T> 
inline void
TriW( Matrix<T>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals of ones must be non-negative");
    const Int numDiags = ( (n>0)&&(m>0) ? m+n-1 : 0 );
    std::vector<T> a( numDiags, 0 );
    if( n > 0 )
        a[n-1] = 1;
    for( Int j=0; j<std::min(n-1,k); ++j )
        a[n-2-j] = alpha;
    Toeplitz( A, m, n, a );
}

#ifndef SWIG
template<typename T> 
inline Matrix<T>
TriW( Int m, Int n, T alpha, Int k )
{
    Matrix<T> A;
    TriW( A, m, n, alpha, k );
    return A;
}
#endif

template<typename T,Dist U,Dist V>
inline void
TriW( DistMatrix<T,U,V>& A, Int m, Int n, T alpha, Int k )
{
    DEBUG_ONLY(CallStackEntry cse("TriW"))
    if( k < 0 )
        LogicError("Number of superdiagonals of ones must be non-negative");
    const Int numDiags = ( (n>0)&&(m>0) ? m+n-1 : 0 );
    std::vector<T> a( numDiags, 0 );
    if( n > 0 )
        a[n-1] = 1;
    for( Int j=0; j<std::min(n-1,k); ++j )
        a[n-2-j] = alpha;
    Toeplitz( A, m, n, a );
}

#ifndef SWIG
template<typename T,Dist U=MC,Dist V=MR>
inline DistMatrix<T,U,V>
TriW( const Grid& g, Int m, Int n, T alpha, Int k )
{
    DistMatrix<T,U,V> A(g);
    TriW( A, m, n, alpha, k );
    return A;
}
#endif

// TODO: MakeTriW?

} // namespace elem

#endif // ifndef ELEM_TRIW_HPP
