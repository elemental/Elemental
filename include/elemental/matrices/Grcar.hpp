/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_GRCAR_HPP
#define ELEM_MATRICES_GRCAR_HPP

#include "elemental/matrices/Toeplitz.hpp"

namespace elem {

template<typename T> 
inline void
Grcar( Matrix<T>& A, Int n, Int k=3 )
{
#ifndef RELEASE
    CallStackEntry cse("Grcar");
#endif
    if( k < 0 )
        LogicError("Number of superdiagonals of ones must be non-negative");
    const Int numDiags = ( n>0 ? 2*n-1 : 0 );
    std::vector<T> a( numDiags, 0 );
    if( n > 1 )
        a[n] = -1;
    if( n > 0 )
        a[n-1] = 1;
    for( Int j=0; j<std::min(n-1,k); ++j )
        a[n-2-j] = 1;
    Toeplitz( A, n, n, a );
}

template<typename T> 
inline Matrix<T>
Grcar( Int n, Int k=3 )
{
    Matrix<T> A;
    Grcar( A, n, k );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void
Grcar( DistMatrix<T,U,V>& A, Int n, Int k=3 )
{
#ifndef RELEASE
    CallStackEntry cse("Grcar");
#endif
    if( k < 0 )
        LogicError("Number of superdiagonals of ones must be non-negative");
    const Int numDiags = ( n>0 ? 2*n-1 : 0 );
    std::vector<T> a( numDiags, 0 );
    if( n > 1 )
        a[n] = -1;
    if( n > 0 )
        a[n-1] = 1;
    for( Int j=0; j<std::min(n-1,k); ++j )
        a[n-2-j] = 1;
    Toeplitz( A, n, n, a );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Grcar( const Grid& g, Int n, Int k=3 )
{
    DistMatrix<T,U,V> A(g);
    Grcar( A, n, k );
    return A;
}

// TODO: MakeGrcar?

} // namespace elem

#endif // ifndef ELEM_MATRICES_GRCAR_HPP
