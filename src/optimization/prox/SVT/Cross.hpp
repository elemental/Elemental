/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVT_CROSS_HPP
#define EL_SVT_CROSS_HPP

namespace El {

namespace svt {

template<typename Field>
Int Cross( Matrix<Field>& A, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    Matrix<Field> U;
    Matrix<Real> s;
    Matrix<Field> V;
    SVDCtrl<Real> ctrl;
    // It is perhaps misleading to have 'approach' stored within 'bidiagSVDCtrl'
    // when the PRODUCT_SVD approach reduces to tridiagonal form instead; we
    // could think of this as implicitly forming the Grammian of the bidiagonal
    // matrix
    //
    // TODO(poulson): A more aggressive tolerance? tol = tau is too strong.
    ctrl.bidiagSVDCtrl.approach = PRODUCT_SVD;
    SVD( A, U, s, V, ctrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, Field(1), U, V, Field(0), A );

    return ZeroNorm( s );
}

template<typename Field>
Int Cross
( AbstractDistMatrix<Field>& APre, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    typedef Base<Field> Real;

    auto& A = AProx.Get();

    DistMatrix<Real,VR,STAR> s( A.Grid() );
    DistMatrix<Field> U( A.Grid() ), V( A.Grid() );
    SVDCtrl<Real> ctrl;
    // See the equivalent note above
    ctrl.bidiagSVDCtrl.approach = PRODUCT_SVD;
    SVD( A, U, s, V, ctrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, Field(1), U, V, Field(0), A );

    return ZeroNorm( s );
}

template<typename Field>
Int Cross
( DistMatrix<Field,VC,STAR>& A, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrix<Field,VC,STAR> U( A.Grid() );
    DistMatrix<Real,STAR,STAR> s( A.Grid() );
    DistMatrix<Field,STAR,STAR> V( A.Grid() );

    SVDCtrl<Real> ctrl;
    // See the equivalent note above
    ctrl.bidiagSVDCtrl.approach = PRODUCT_SVD;
    SVD( A, U, s, V, ctrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    LocalGemm( NORMAL, ADJOINT, Field(1), U, V, Field(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace El

#endif // ifndef EL_SVT_CROSS_HPP
