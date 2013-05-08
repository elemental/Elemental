/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_GRCAR_HPP
#define MATRICES_GRCAR_HPP

#include "elemental/matrices/Toeplitz.hpp"

namespace elem {

template<typename T> 
inline void
Grcar( Matrix<T>& A, int n, int k=3 )
{
#ifndef RELEASE
    CallStackEntry entry("Grcar");
#endif
    if( k < 0 )
        throw std::logic_error
        ("Number of superdiagonals of ones must be non-negative");
    const int numDiags = ( n>0 ? 2*n-1 : 0 );
    std::vector<T> a( numDiags, 0 );
    if( n > 1 )
        a[n-2] = -1;
    if( n > 0 )
        a[n-1] = 1;
    for( int j=0; j<std::min(n-1,k); ++j )
        a[j+n] = 1;
    Toeplitz( A, n, n, a );
}

template<typename T,Distribution U,Distribution V>
inline void
Grcar( DistMatrix<T,U,V>& A, int n, int k=3 )
{
#ifndef RELEASE
    CallStackEntry entry("Grcar");
#endif
    if( k < 0 )
        throw std::logic_error
        ("Number of superdiagonals of ones must be non-negative");
    const int numDiags = ( n>0 ? 2*n-1 : 0 );
    std::vector<T> a( numDiags, 0 );
    if( n > 1 )
        a[n-2] = -1;
    if( n > 0 )
        a[n-1] = 1;
    for( int j=0; j<std::min(n-1,k); ++j )
        a[j+n] = 1;
    Toeplitz( A, n, n, a );
}

// TODO: MakeGrcar?

} // namespace elem

#endif // ifndef MATRICES_GRCAR_HPP
