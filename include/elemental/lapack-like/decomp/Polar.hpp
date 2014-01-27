/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_POLAR_HPP
#define ELEM_POLAR_HPP

#include ELEM_SIGN_INC

#include "./Polar/QDWH.hpp"
#include "./Polar/SVD.hpp"

namespace elem {

// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is 
// Hermitian positive semi-definite. On exit, A is overwritten with Q.

template<typename F>
inline void
Polar( Matrix<F>& A, Matrix<F>& P )
{
    DEBUG_ONLY(CallStackEntry cse("Polar"))
    polar::SVD( A, P );
}

template<typename F>
inline void
Polar( DistMatrix<F>& A, DistMatrix<F>& P )
{
    DEBUG_ONLY(CallStackEntry cse("Polar"))
    polar::SVD( A, P );
}

template<typename F>
inline void
HermitianPolar( UpperOrLower uplo, Matrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    HermitianSign( uplo, A );
}

template<typename F>
inline void
HermitianPolar( UpperOrLower uplo, DistMatrix<F>& A )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    HermitianSign( uplo, A );
}

template<typename F>
inline void
HermitianPolar( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    HermitianSign( uplo, A, P );
}

template<typename F>
inline void
HermitianPolar( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& P )
{
    DEBUG_ONLY(CallStackEntry cse("HermitianPolar"))
    HermitianSign( uplo, A, P );
}

} // namespace elem

#endif // ifndef ELEM_POLAR_HPP
