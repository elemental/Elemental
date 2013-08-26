/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_FORSYTHE_HPP
#define ELEM_MATRICES_FORSYTHE_HPP

#include "elemental/matrices/Jordan.hpp"

namespace elem {

template<typename T> 
inline void
MakeForsythe( Matrix<T>& J, T alpha, T lambda )
{
#ifndef RELEASE
    CallStackEntry cse("MakeForsythe");
#endif
    MakeJordan( J, lambda );
    const Int m = J.Height();
    const Int n = J.Width();
    if( m > 0 && n > 0 )
        J.Set( m-1, 0, alpha );
}

template<typename T,Distribution U,Distribution V>
inline void
MakeForsythe( DistMatrix<T,U,V>& J, T alpha, T lambda )
{
#ifndef RELEASE
    CallStackEntry cse("MakeForsythe");
#endif
    MakeJordan( J, lambda );
    const Int m = J.Height();
    const Int n = J.Width();
    if( m > 0 && n > 0 )
        J.Set( m-1, 0, alpha );
}

template<typename T>
inline Matrix<T>
Forsythe( Int n, T alpha, T lambda )
{
    Matrix<T> J( n, n );
    MakeForsythe( J, alpha, lambda );
    return J;
}

template<typename T,Distribution U,Distribution V>
inline void
Forsythe( DistMatrix<T,U,V>& J, Int n, T alpha, T lambda )
{
#ifndef RELEASE
    CallStackEntry cse("Forsythe");
#endif
    J.ResizeTo( n, n );
    MakeForsythe( J, alpha, lambda );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Forsythe( const Grid& g, Int n, T alpha, T lambda )
{
    DistMatrix<T,U,V> J( n, n, g );
    MakeForsythe( J, alpha, lambda );
    return J;
}

} // namespace elem

#endif // ifndef ELEM_MATRICES_FORSYTHE_HPP
