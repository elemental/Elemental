/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_POLAR_SVD_HPP
#define LAPACK_POLAR_SVD_HPP

#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/lapack-like/HermitianFunction.hpp"
#include "elemental/lapack-like/SVD.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace polar {

//
// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is 
// Hermitian positive semi-definite. On exit, A is overwritten with Q.
//

template<typename F>
inline void
SVD( Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("polar::SVD");
#endif
    // Get the SVD of A
    typedef BASE(F) R;
    Matrix<R> s;
    Matrix<F> U, V;
    U = A;
    elem::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );
}

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<F>& P )
{
#ifndef RELEASE
    CallStackEntry entry("polar::SVD");
#endif
    // Get the SVD of A
    typedef BASE(F) R;
    Matrix<R> s;
    Matrix<F> U, V;
    U = A;
    elem::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );

    // Form P := V Sigma V^H in P
    HermitianFromEVD( LOWER, P, s, V );
}

template<typename F>
inline void
SVD( DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("polar::SVD");
#endif
    // Get the SVD of A
    typedef BASE(F) R;
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    elem::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );
}

template<typename F>
inline void
SVD( DistMatrix<F>& A, DistMatrix<F>& P )
{
#ifndef RELEASE
    CallStackEntry entry("polar::SVD");
#endif
    // Get the SVD of A
    typedef BASE(F) R;
    const Grid& g = A.Grid();
    DistMatrix<R,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    elem::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );

    // Form P := V Sigma V^H in P
    HermitianFromEVD( LOWER, P, s, V );
}

} // namespace polar
} // namespace elem

#endif // ifndef LAPACK_POLAR_SVD_HPP
