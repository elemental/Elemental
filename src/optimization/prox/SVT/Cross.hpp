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

template<typename F>
Int Cross( Matrix<F>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CSE cse("svt::Cross"))
    typedef Base<F> Real;

    Matrix<F> U;
    Matrix<Real> s;
    Matrix<F> V;
    SVDCtrl<Real> ctrl;
    ctrl.approach = PRODUCT_SVD;
    ctrl.tol = tau;
    ctrl.relative = relative;
    SVD( A, U, s, V, ctrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

template<typename F>
Int Cross( ElementalMatrix<F>& APre, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CSE cse("svt::Cross"))

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    typedef Base<F> Real;

    auto& A = AProx.Get();

    DistMatrix<Real,VR,STAR> s( A.Grid() );
    DistMatrix<F> U( A.Grid() ), V( A.Grid() );
    SVDCtrl<Real> ctrl;
    ctrl.approach = PRODUCT_SVD;
    ctrl.tol = tau;
    ctrl.relative = relative;
    SVD( A, U, s, V, ctrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

template<typename F>
Int Cross( DistMatrix<F,VC,STAR>& A, Base<F> tau, bool relative )
{
    DEBUG_ONLY(CSE cse("svt::Cross"))
    typedef Base<F> Real;

    DistMatrix<F,VC,STAR> U( A.Grid() );
    DistMatrix<Real,STAR,STAR> s( A.Grid() );
    DistMatrix<F,STAR,STAR> V( A.Grid() );

    SVDCtrl<Real> ctrl;
    ctrl.approach = PRODUCT_SVD;
    ctrl.tol = tau;
    ctrl.relative = relative;
    SVD( A, U, s, V, ctrl );

    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    LocalGemm( NORMAL, ADJOINT, F(1), U, V, F(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace El

#endif // ifndef EL_SVT_CROSS_HPP
