/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_TRIW_HPP
#define MATRICES_TRIW_HPP

#include "elemental/matrices/Toeplitz.hpp"

namespace elem {

template<typename T> 
inline void
TriW( Matrix<T>& A, int m, int n, T alpha, int k )
{
#ifndef RELEASE
    CallStackEntry entry("TriW");
#endif
    if( k < 0 )
        throw std::logic_error
        ("Number of superdiagonals of ones must be non-negative");
    const int numDiags = ( (n>0)&&(m>0) ? m+n-1 : 0 );
    std::vector<T> a( numDiags, 0 );
    if( n > 0 )
        a[n-1] = 1;
    for( int j=0; j<std::min(n-1,k); ++j )
        a[n-2-j] = alpha;
    Toeplitz( A, m, n, a );
}

template<typename T,Distribution U,Distribution V>
inline void
TriW( DistMatrix<T,U,V>& A, int m, int n, T alpha, int k )
{
#ifndef RELEASE
    CallStackEntry entry("TriW");
#endif
    if( k < 0 )
        throw std::logic_error
        ("Number of superdiagonals of ones must be non-negative");
    const int numDiags = ( (n>0)&&(m>0) ? m+n-1 : 0 );
    std::vector<T> a( numDiags, 0 );
    if( n > 0 )
        a[n-1] = 1;
    for( int j=0; j<std::min(n-1,k); ++j )
        a[n-2-j] = alpha;
    Toeplitz( A, m, n, a );
}

// TODO: MakeTriW?

} // namespace elem

#endif // ifndef MATRICES_TRIW_HPP
