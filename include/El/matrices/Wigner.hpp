/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_WIGNER_HPP
#define EL_WIGNER_HPP

#include EL_GAUSSIAN_INC

namespace El {

template<typename T>
inline void
Wigner( Matrix<T>& A, Int n, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Wigner"))
    Gaussian( A, n, n, mean, stddev );
    MakeHermitian( LOWER, A );
}

template<typename T>
inline void
Wigner( AbstractDistMatrix<T>& A, Int n, T mean=0, Base<T> stddev=1 )
{
    DEBUG_ONLY(CallStackEntry cse("Wigner"))
    Gaussian( A, n, n, mean, stddev );
    MakeHermitian( LOWER, A );
}

} // namespace El

#endif // ifndef EL_WIGNER_HPP
