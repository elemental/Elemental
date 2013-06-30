/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_POLAR_HPP
#define LAPACK_POLAR_HPP

#include "elemental/lapack-like/Polar/SVD.hpp"
#include "elemental/lapack-like/Polar/QDWH.hpp"
#include "elemental/lapack-like/Sign.hpp"

namespace elem {

//
// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is 
// Hermitian positive semi-definite. On exit, A is overwritten with Q.
//

template<typename F>
inline void
Polar( Matrix<F>& A, Matrix<F>& P )
{
#ifndef RELEASE
    CallStackEntry entry("Polar");
#endif
    polar::SVD( A, P );
}

template<typename F>
inline void
Polar( DistMatrix<F>& A, DistMatrix<F>& P )
{
#ifndef RELEASE
    CallStackEntry entry("Polar");
#endif
    polar::SVD( A, P );
}

template<typename F>
inline void
HermitianPolar( UpperOrLower uplo, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianPolar");
#endif
    HermitianSign( uplo, A );
}

template<typename F>
inline void
HermitianPolar( UpperOrLower uplo, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianPolar");
#endif
    HermitianSign( uplo, A );
}

template<typename F>
inline void
HermitianPolar( UpperOrLower uplo, Matrix<F>& A, Matrix<F>& P )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianPolar");
#endif
    HermitianSign( uplo, A, P );
}

template<typename F>
inline void
HermitianPolar( UpperOrLower uplo, DistMatrix<F>& A, DistMatrix<F>& P )
{
#ifndef RELEASE
    CallStackEntry entry("HermitianPolar");
#endif
    HermitianSign( uplo, A, P );
}

} // namespace elem

#endif // ifndef LAPACK_POLAR_HPP
