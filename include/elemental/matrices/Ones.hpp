/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_ONES_HPP
#define ELEM_ONES_HPP

namespace elem {

template<typename T> 
inline void
MakeOnes( Matrix<T>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    const Int m = A.Height();
    const Int n = A.Width();
    for( Int j=0; j<n; ++j )
        for( Int i=0; i<m; ++i )
            A.Set( i, j, T(1) );
}

template<typename T,Dist U,Dist V>
inline void
MakeOnes( DistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    MakeOnes( A.Matrix() );
}

template<typename T,Dist U,Dist V>
inline void
MakeOnes( BlockDistMatrix<T,U,V>& A )
{
    DEBUG_ONLY(CallStackEntry cse("MakeOnes"))
    MakeOnes( A.Matrix() );
}

template<typename T>
inline void
Ones( Matrix<T>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

template<typename T,Dist U,Dist V>
inline void
Ones( DistMatrix<T,U,V>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

template<typename T,Dist U,Dist V>
inline void
Ones( BlockDistMatrix<T,U,V>& A, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ones"))
    A.Resize( m, n );
    MakeOnes( A );
}

} // namespace elem

#endif // ifndef ELEM_ONES_HPP
