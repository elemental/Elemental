/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LOTKIN_HPP
#define ELEM_LOTKIN_HPP

#include ELEM_HILBERT_INC

namespace elem {

template<typename F>
inline void
Lotkin( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Lotkin"))
    Hilbert( A, n );
    // Set first row to all ones
    for( Int j=0; j<n; ++j )
        A.Set( 0, j, F(1) );
}

#ifndef SWIG
template<typename F>
inline Matrix<F>
Lotkin( Int n )
{
    Matrix<F> A;
    Lotkin( A, n );
    return A;
}
#endif

template<typename F,Dist U,Dist V>
inline void
Lotkin( DistMatrix<F,U,V>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Lotkin"))
    Hilbert( A, n );
    // Set first row to all ones
    if( A.ColShift() == 0 )
    {
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            A.SetLocal( 0, jLoc, F(1) );
    } 
}

#ifndef SWIG
template<typename F,Dist U=MC,Dist V=MR>
inline DistMatrix<F,U,V>
Lotkin( const Grid& g, Int n )
{
    DistMatrix<F,U,V> A(g);
    Lotkin( A, n );
    return A;
}
#endif

// TODO: MakeLotkin?

} // namespace elem

#endif // ifndef ELEM_LOTKIN_HPP
