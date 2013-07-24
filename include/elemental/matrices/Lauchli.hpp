/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_MATRICES_LAUCHLI_HPP
#define ELEM_MATRICES_LAUCHLI_HPP

#include "elemental/matrices/Diagonal.hpp"
#include "elemental/matrices/Ones.hpp"

namespace elem {

template<typename T>
inline void
Lauchli( Matrix<T>& A, int n, T mu )
{
#ifndef RELEASE
    CallStackEntry entry("Lauchli");
#endif
    A.ResizeTo( n+1, n );
    Matrix<T> ABlock;

    View( ABlock, A, 0, 0, 1, n );
    MakeOnes( ABlock );

    std::vector<T> d(n,mu);
    View( ABlock, A, 1, 0, n, n );
    Diagonal( ABlock, d );
}

template<typename T,Distribution U,Distribution V>
inline void
Lauchli( DistMatrix<T,U,V>& A, int n, T mu )
{
#ifndef RELEASE
    CallStackEntry entry("Lauchli");
#endif
    A.ResizeTo( n+1, n );
    DistMatrix<T,U,V> ABlock( A.Grid() );

    View( ABlock, A, 0, 0, 1, n );
    MakeOnes( ABlock );

    std::vector<T> d(n,mu);
    View( ABlock, A, 1, 0, n, n );
    Diagonal( ABlock, d );
}

// TODO: MakeLauchli?

} // namespace elem

#endif // ifndef ELEM_MATRICES_LAUCHLI_HPP
