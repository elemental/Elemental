/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_ZEROS_HPP
#define ELEM_MATRICES_ZEROS_HPP

#include "elemental/blas-like/level1/Zero.hpp"

namespace elem {

template<typename T> 
inline void
MakeZeros( Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeZeros");
#endif
    Zero( A );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeZeros( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry cse("MakeZeros");
#endif
    Zero( A.Matrix() );
}

template<typename T>
inline void
Zeros( Matrix<T>& A, Int m, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Zeros");
#endif
    A.ResizeTo( m, n );
    MakeZeros( A );
}

template<typename T>
inline Matrix<T>
Zeros( Int m, Int n )
{
    Matrix<T> A( m, n );
    MakeZeros( A ); 
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void
Zeros( DistMatrix<T,U,V>& A, Int m, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Zeros");
#endif
    A.ResizeTo( m, n );
    MakeZeros( A );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Zeros( const Grid& g, Int m, Int n )
{
    DistMatrix<T,U,V> A( m, n, g );
    MakeZeros( A );
    return A;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_ZEROS_HPP
