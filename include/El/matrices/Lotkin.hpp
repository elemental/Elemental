/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LOTKIN_HPP
#define EL_LOTKIN_HPP

#include EL_HILBERT_INC

namespace El {

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

template<typename F>
inline void
Lotkin( AbstractDistMatrix<F>& A, Int n )
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

template<typename F>
inline void
Lotkin( AbstractBlockDistMatrix<F>& A, Int n )
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

} // namespace El

#endif // ifndef EL_LOTKIN_HPP
