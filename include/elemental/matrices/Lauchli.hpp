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
Lauchli( Matrix<T>& A, Int n, T mu )
{
#ifndef RELEASE
    CallStackEntry cse("Lauchli");
#endif
    A.ResizeTo( n+1, n );

    auto ABlock = View( A, 0, 0, 1, n );
    MakeOnes( ABlock );

    std::vector<T> d(n,mu);
    ABlock = View( A, 1, 0, n, n );
    Diagonal( ABlock, d );
}

template<typename T>
inline Matrix<T>
Lauchli( Int n, T mu )
{
    Matrix<T> A;
    Lauchli( A, n, mu );
    return A;
}

template<typename T,Distribution U,Distribution V>
inline void
Lauchli( DistMatrix<T,U,V>& A, Int n, T mu )
{
#ifndef RELEASE
    CallStackEntry cse("Lauchli");
#endif
    A.ResizeTo( n+1, n );

    auto ABlock = View( A, 0, 0, 1, n );
    MakeOnes( ABlock );

    std::vector<T> d(n,mu);
    ABlock = View( A, 1, 0, n, n );
    Diagonal( ABlock, d );
}

template<typename T,Distribution U=MC,Distribution V=MR>
inline DistMatrix<T,U,V>
Lauchli( const Grid& g, Int n, T mu )
{
    DistMatrix<T,U,V> A(g);
    Lauchli( A, n, mu );
    return A;
}

// TODO: MakeLauchli?

} // namespace elem

#endif // ifndef ELEM_MATRICES_LAUCHLI_HPP
