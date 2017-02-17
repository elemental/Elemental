/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_POLAR_SVD_HPP
#define EL_POLAR_SVD_HPP

namespace El {
namespace polar {

// Compute the polar decomposition of A, A = Q P, where Q is unitary and P is
// Hermitian positive semi-definite. On exit, A is overwritten with Q.

template<typename Field>
void SVD( Matrix<Field>& A )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    // Get the SVD of A
    Matrix<Real> s;
    Matrix<Field> U, V;
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    El::SVD( A, U, s, V, ctrl );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, Field(1), U, V, A );
}

template<typename Field>
void SVD( Matrix<Field>& A, Matrix<Field>& P )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    // Get the SVD of A
    Matrix<Real> s;
    Matrix<Field> U, V;
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    El::SVD( A, U, s, V, ctrl );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, Field(1), U, V, A );

    // Form P := V Sigma V^H in P
    HermitianFromEVD( LOWER, P, s, V );
}

template<typename Field>
void SVD( AbstractDistMatrix<Field>& APre )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();
    const Grid& g = A.Grid();

    // Get the SVD of A
    DistMatrix<Real,VR,STAR> s(g);
    DistMatrix<Field> U(g), V(g);
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    El::SVD( A, U, s, V, ctrl );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, Field(1), U, V, A );
}

template<typename Field>
void SVD( AbstractDistMatrix<Field>& APre, AbstractDistMatrix<Field>& PPre )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    DistMatrixWriteProxy<Field,Field,MC,MR> PProx( PPre );
    auto& A = AProx.Get();
    auto& P = PProx.Get();
    const Grid& g = A.Grid();

    // Get the SVD of A
    DistMatrix<Real,VR,STAR> s(g);
    DistMatrix<Field> U(g), V(g);
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    El::SVD( A, U, s, V, ctrl );

    // Form Q := U V^H in A
    Gemm( NORMAL, ADJOINT, Field(1), U, V, A );

    // Form P := V Sigma V^H in P
    HermitianFromEVD( LOWER, P, s, V );
}

} // namespace polar
} // namespace El

#endif // ifndef EL_POLAR_SVD_HPP
