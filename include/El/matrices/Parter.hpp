/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_PARTER_HPP
#define EL_PARTER_HPP

namespace El {

template<typename F> 
inline void
Parter( Matrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    P.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( P, [=]( Int i, Int j ) { return F(1)/(F(i)-F(j)+oneHalf); } );
}

template<typename F>
inline void
Parter( AbstractDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    P.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( P, [=]( Int i, Int j ) { return F(1)/(F(i)-F(j)+oneHalf); } );
}

template<typename F>
inline void
Parter( AbstractBlockDistMatrix<F>& P, Int n )
{
    DEBUG_ONLY(CallStackEntry cse("Parter"))
    P.Resize( n, n );
    const F oneHalf = F(1)/F(2);
    IndexDependentFill
    ( P, [=]( Int i, Int j ) { return F(1)/(F(i)-F(j)+oneHalf); } );
}

} // namespace El

#endif // ifndef EL_PARTER_HPP
