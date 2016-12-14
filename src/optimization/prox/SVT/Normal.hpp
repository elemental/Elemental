/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SVT_NORMAL_HPP
#define EL_SVT_NORMAL_HPP

namespace El {

namespace svt {

template<typename Field>
Int Normal( Matrix<Field>& A, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    Matrix<Field> U, V;
    Matrix<Real> s;
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    SVD( A, U, s, V, ctrl );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, Field(1), U, V, Field(0), A );

    return ZeroNorm( s );
}

template<typename Field>
Int Normal
( AbstractDistMatrix<Field>& APre, const Base<Field>& tau, bool relative )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrix<Real,VR,STAR> s( A.Grid() );
    DistMatrix<Field> U( A.Grid() ), V( A.Grid() );
    SVDCtrl<Real> ctrl;
    ctrl.overwrite = true;
    SVD( A, U, s, V );
    SoftThreshold( s, tau, relative );
    DiagonalScale( RIGHT, NORMAL, s, U );
    Gemm( NORMAL, ADJOINT, Field(1), U, V, Field(0), A );

    return ZeroNorm( s );
}

} // namespace svt
} // namespace El

#endif // ifndef EL_SVT_NORMAL_HPP
