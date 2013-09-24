/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_LOTKIN_HPP
#define ELEM_MATRICES_LOTKIN_HPP

#include "elemental/matrices/Hilbert.hpp"

namespace elem {

template<typename F>
inline void
Lotkin( Matrix<F>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Lotkin");
#endif
    Hilbert( A, n );
    // Set first row to all ones
    for( Int j=0; j<n; ++j )
        A.Set( 0, j, F(1) );
}

template<typename F>
inline Matrix<F>
Lotkin( Int n )
{
    Matrix<F> A;
    Lotkin( A, n );
    return A;
}

template<typename F,Distribution U,Distribution V>
inline void
Lotkin( DistMatrix<F,U,V>& A, Int n )
{
#ifndef RELEASE
    CallStackEntry cse("Lotkin");
#endif
    Hilbert( A, n );
    // Set first row to all ones
    if( A.ColShift() == 0 )
    {
        const Int localWidth = A.LocalWidth();
        for( Int jLoc=0; jLoc<localWidth; ++jLoc )
            A.SetLocal( 0, jLoc, F(1) );
    } 
}

template<typename F,Distribution U=MC,Distribution V=MR>
inline DistMatrix<F,U,V>
Lotkin( const Grid& g, Int n )
{
    DistMatrix<F,U,V> A(g);
    Lotkin( A, n );
    return A;
}

// TODO: MakeLotkin?

} // namespace elem

#endif // ifndef ELEM_MATRICES_LOTKIN_HPP
