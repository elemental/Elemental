/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_POLAR_SVD_HPP
#define EL_POLAR_SVD_HPP

namespace El {
namespace polar {

// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is 
// Hermitian positive semi-definite. On exit, A is overwritten with Q.

template<typename F>
inline void
SVD( Matrix<F>& A )
{
    DEBUG_ONLY(CSE cse("polar::SVD"))
    // Get the SVD of A
    typedef Base<F> Real;
    Matrix<Real> s;
    Matrix<F> U, V;
    U = A;
    El::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );
}

template<typename F>
inline void
SVD( Matrix<F>& A, Matrix<F>& P )
{
    DEBUG_ONLY(CSE cse("polar::SVD"))
    // Get the SVD of A
    typedef Base<F> Real;
    Matrix<Real> s;
    Matrix<F> U, V;
    U = A;
    El::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );

    // Form P := V Sigma V^H in P
    HermitianFromEVD( LOWER, P, s, V );
}

template<typename F>
inline void
SVD( ElementalMatrix<F>& APre )
{
    DEBUG_ONLY(CSE cse("polar::SVD"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    // Get the SVD of A
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    El::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );
}

template<typename F>
inline void
SVD( ElementalMatrix<F>& APre, ElementalMatrix<F>& PPre )
{
    DEBUG_ONLY(CSE cse("polar::SVD"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> PProx( PPre );
    auto& A = AProx.Get();
    auto& P = PProx.Get();

    // Get the SVD of A
    typedef Base<F> Real;
    const Grid& g = A.Grid();
    DistMatrix<Real,VR,STAR> s(g);
    DistMatrix<F> U(g), V(g);
    U = A;
    El::SVD( U, s, V );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, F(1), U, V, A );

    // Form P := V Sigma V^H in P
    HermitianFromEVD( LOWER, P, s, V );
}

} // namespace polar
} // namespace El

#endif // ifndef EL_POLAR_SVD_HPP
