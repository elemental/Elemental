/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef MATRICES_ZEROS_HPP
#define MATRICES_ZEROS_HPP

#include "elemental/blas-like/level1/Zero.hpp"

namespace elem {

template<typename T>
inline void
Zeros( Matrix<T>& A, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Zeros");
#endif
    A.ResizeTo( m, n );
    MakeZeros( A );
}

template<typename T,Distribution U,Distribution V>
inline void
Zeros( DistMatrix<T,U,V>& A, int m, int n )
{
#ifndef RELEASE
    CallStackEntry entry("Zeros");
#endif
    A.ResizeTo( m, n );
    MakeZeros( A );
}

template<typename T> 
inline void
MakeZeros( Matrix<T>& A )
{
#ifndef RELEASE
    CallStackEntry entry("MakeZeros");
#endif
    Zero( A );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeZeros( DistMatrix<T,U,V>& A )
{
#ifndef RELEASE
    CallStackEntry entry("MakeZeros");
#endif
    Zero( A.Matrix() );
}

} // namespace elem

#endif // ifndef MATRICES_ZEROS_HPP
