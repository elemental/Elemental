/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LAUCHLI_HPP
#define EL_LAUCHLI_HPP

#include EL_ONES_INC

namespace El {

template<typename T>
inline void
Lauchli( Matrix<T>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Lauchli"))
    A.Resize( n+1, n );

    auto ABlock = View( A, 0, 0, 1, n );
    MakeOnes( ABlock );

    std::vector<T> d(n,mu);
    ABlock = View( A, 1, 0, n, n );
    Diagonal( ABlock, d );
}

template<typename T,Dist U,Dist V>
inline void
Lauchli( DistMatrix<T,U,V>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Lauchli"))
    A.Resize( n+1, n );

    auto ABlock = View( A, 0, 0, 1, n );
    MakeOnes( ABlock );

    std::vector<T> d(n,mu);
    ABlock = View( A, 1, 0, n, n );
    Diagonal( ABlock, d );
}

template<typename T,Dist U,Dist V>
inline void
Lauchli( BlockDistMatrix<T,U,V>& A, Int n, T mu )
{
    DEBUG_ONLY(CallStackEntry cse("Lauchli"))
    A.Resize( n+1, n );

    auto ABlock = View( A, 0, 0, 1, n );
    MakeOnes( ABlock );

    std::vector<T> d(n,mu);
    ABlock = View( A, 1, 0, n, n );
    Diagonal( ABlock, d );
}

} // namespace El

#endif // ifndef EL_LAUCHLI_HPP
