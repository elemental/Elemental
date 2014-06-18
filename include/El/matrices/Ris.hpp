/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_RIS_HPP
#define EL_RIS_HPP

namespace El {

template<typename F> 
inline void
Ris( Matrix<F>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    R.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( R, [=]( Int i, Int j ) { return oneHalf/(F(n-i-j)-oneHalf); } );
}

template<typename F>
inline void
Ris( AbstractDistMatrix<F>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    R.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( R, [=]( Int i, Int j ) { return oneHalf/(F(n-i-j)-oneHalf); } );
}

template<typename F>
inline void
Ris( AbstractBlockDistMatrix<F>& R, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Ris"))
    R.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( R, [=]( Int i, Int j ) { return oneHalf/(F(n-i-j)-oneHalf); } );
}

} // namespace El

#endif // ifndef EL_RIS_HPP
