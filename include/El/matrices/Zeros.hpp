/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_ZEROS_HPP
#define EL_ZEROS_HPP

namespace El {

template<typename T> 
inline void
MakeZeros( Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeZeros"))
    Zero( A );
}

template<typename T>
inline void
MakeZeros( AbstractDistMatrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeZeros"))
    Zero( A.Matrix() );
}

template<typename T>
inline void
MakeZeros( AbstractBlockDistMatrix<T>& A )
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

template<typename T>
inline void
Zeros( AbstractDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Zeros"))
    A.Resize( m, n );
    MakeZeros( A );
}

template<typename T>
inline void
Zeros( AbstractBlockDistMatrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Zeros"))
    A.Resize( m, n );
    MakeZeros( A );
}

} // namespace El

#endif // ifndef EL_ZEROS_HPP
