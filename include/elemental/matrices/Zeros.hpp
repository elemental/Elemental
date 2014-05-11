/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_ZEROS_HPP
#define ELEM_ZEROS_HPP

#include ELEM_ZERO_INC

namespace elem {

template<typename T> 
inline void
MakeZeros( Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeZeros"))
    Zero( A );
}

template<typename T,Dist U,Dist V>
inline void
MakeZeros( DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeZeros"))
    Zero( A.Matrix() );
}

template<typename T,Dist U,Dist V>
inline void
MakeZeros( BlockDistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeZeros"))
    Zero( A.Matrix() );
}

template<typename T>
inline void
Zeros( Matrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Zeros"))
    A.Resize( m, n );
    MakeZeros( A );
}

template<typename T,Dist U,Dist V>
inline void
Zeros( DistMatrix<T,U,V>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Zeros"))
    A.Resize( m, n );
    MakeZeros( A );
}

template<typename T,Dist U,Dist V>
inline void
Zeros( BlockDistMatrix<T,U,V>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Zeros"))
    A.Resize( m, n );
    MakeZeros( A );
}

} // namespace elem

#endif // ifndef ELEM_ZEROS_HPP
