/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_HILBERT_HPP
#define EL_HILBERT_HPP

namespace El {

template<typename F>
inline void
Hilbert( Matrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Hilbert"))
    A.Resize( n, n );
    IndexDependentFill( A, []( Int i, Int j ) { return F(1)/F(i+j+1); } );
}

template<typename F>
inline void
Hilbert( AbstractDistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Hilbert"))
    A.Resize( n, n );
    IndexDependentFill( A, []( Int i, Int j ) { return F(1)/F(i+j+1); } );
}

template<typename F>
inline void
Hilbert( AbstractBlockDistMatrix<F>& A, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Hilbert"))
    A.Resize( n, n );
    IndexDependentFill( A, []( Int i, Int j ) { return F(1)/F(i+j+1); } );
}

} // namespace El

#endif // ifndef EL_HILBERT_HPP
