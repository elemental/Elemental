/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_DOT_HPP
#define LAPACK_DOT_HPP

#include "elemental/lapack-like/HilbertSchmidt.hpp"

namespace elem {

template<typename F> 
inline F
Dot( const Matrix<F>& A, const Matrix<F>& B )
{
#ifndef RELEASE
    PushCallStack("Dot");
#endif
    const F innerProd = HilbertSchmidt( A, B );
#ifndef RELEASE
    PopCallStack();
#endif
    return innerProd;
}

template<typename F,Distribution U,Distribution V> 
inline F
Dot( const DistMatrix<F,U,V>& A, const DistMatrix<F,U,V>& B )
{
#ifndef RELEASE
    PushCallStack("Dot");
#endif
    const F innerProd = HilbertSchmidt( A, B );
#ifndef RELEASE
    PopCallStack();
#endif
    return innerProd;
}

} // namespace elem

#endif // ifndef LAPACK_DOT_HPP
