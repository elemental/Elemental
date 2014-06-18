/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_MINIJ_HPP
#define EL_MINIJ_HPP

namespace El {

template<typename T> 
inline void
MinIJ( Matrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    IndexDependentFill( M, []( Int i, Int j ) { return T(Min(i+1,j+1)); } );
}

template<typename T>
inline void
MinIJ( AbstractDistMatrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    IndexDependentFill( M, []( Int i, Int j ) { return T(Min(i+1,j+1)); } );
}

template<typename T>
inline void
MinIJ( AbstractBlockDistMatrix<T>& M, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("MinIJ"))
    M.Resize( n, n );
    IndexDependentFill( M, []( Int i, Int j ) { return T(Min(i+1,j+1)); } );
}

} // namespace El

#endif // ifndef EL_MINIJ_HPP
