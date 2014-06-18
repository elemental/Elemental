/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_GCDMATRIX_HPP
#define EL_GCDMATRIX_HPP

namespace El {

template<typename T>
inline void
GCDMatrix( Matrix<T>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    IndexDependentFill( G, []( Int i, Int j ) { return T(GCD(i+1,j+1)); } );
}

template<typename T>
inline void
GCDMatrix( AbstractDistMatrix<T>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    IndexDependentFill( G, []( Int i, Int j ) { return T(GCD(i+1,j+1)); } );
}

template<typename T>
inline void
GCDMatrix( AbstractBlockDistMatrix<T>& G, Int m, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("GCDMatrix"))
    G.Resize( m, n );
    IndexDependentFill( G, []( Int i, Int j ) { return T(GCD(i+1,j+1)); } );
}

} // namespace El

#endif // ifndef EL_GCDMATRIX_HPP
