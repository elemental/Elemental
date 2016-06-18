/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_QR_HPP
#define EL_CHOLESKY_QR_HPP

namespace El {
namespace qr {

// NOTE: This version is designed for tall-skinny matrices and is much less
//       numerically stable than Householder-based QR factorizations
//
// Computes the QR factorization of full-rank tall-skinny matrix A and 
// overwrites A with Q
//

template<typename F> 
void Cholesky( Matrix<F>& A, Matrix<F>& R )
{
    DEBUG_CSE
    if( A.Height() < A.Width() )
        LogicError("A^H A will be singular");
    Herk( UPPER, ADJOINT, Base<F>(1), A, R );
    El::Cholesky( UPPER, R );
    Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), R, A );
}

template<typename F> 
void Cholesky( ElementalMatrix<F>& APre, ElementalMatrix<F>& RPre )
{
    DEBUG_CSE
    const Int m = APre.Height();
    const Int n = APre.Width();
    if( m < n )
        LogicError("A^H A will be singular");

    DistMatrixReadWriteProxy<F,F,VC,STAR> AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR> RProx( RPre );
    auto& A = AProx.Get();
    auto& R = RProx.Get();

    Zeros( R, n, n );
    Herk( UPPER, ADJOINT, Base<F>(1), A.Matrix(), Base<F>(0), R.Matrix() );
    El::AllReduce( R, A.ColComm() );
    El::Cholesky( UPPER, R.Matrix() );
    Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), R.Matrix(), A.Matrix() );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_CHOLESKY_HPP
