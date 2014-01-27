/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_DOT_HPP
#define ELEM_DOT_HPP

#include "./HilbertSchmidt.hpp"

namespace elem {

template<typename F> 
inline F
Dot( const Matrix<F>& A, const Matrix<F>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dot"))
    return HilbertSchmidt( A, B );
}

template<typename F,Dist U,Dist V> 
inline F
Dot( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B )
{
    DEBUG_ONLY(CallStackEntry cse("Dot"))
    return HilbertSchmidt( A, B );
}

} // namespace elem

#endif // ifndef ELEM_DOT_HPP
