/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_WIGNER_HPP
#define ELEM_WIGNER_HPP

#include ELEM_MAKEHERMITIAN_INC
#include ELEM_GAUSSIAN_INC

namespace elem {

template<typename T>
inline void
Wigner( Matrix<T>& A, Int n, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Wigner"))
    Gaussian( A, n, n, mean, stddev );
    MakeHermitian( LOWER, A );
}

template<typename T,Dist U,Dist V>
inline void
Wigner( DistMatrix<T,U,V>& A, Int n, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Wigner"))
    Gaussian( A, n, n, mean, stddev );
    MakeHermitian( LOWER, A );
}

} // namespace elem

#endif // ifndef ELEM_WIGNER_HPP
