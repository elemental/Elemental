/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CHOLESKY_QR_HPP
#define EL_CHOLESKY_QR_HPP

#include EL_ZEROS_INC

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
    DEBUG_ONLY(CallStackEntry cse("qr::Cholesky"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( height < width )
        LogicError("A^H A will be singular");
    Herk( UPPER, ADJOINT, F(1), A, R );
    El::Cholesky( UPPER, R );
    Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), R, A );
}

template<typename F> 
void Cholesky( DistMatrix<F,VC,STAR>& A, DistMatrix<F,STAR,STAR>& R )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Cholesky"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( height < width )
        LogicError("A^H A will be singular");
    Zeros( R, width, width );
    Herk( UPPER, ADJOINT, F(1), A.Matrix(), F(0), R.Matrix() );
    R.SumOver( A.ColComm() );
    El::Cholesky( UPPER, R.Matrix() );
    Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), R.Matrix(), A.Matrix() );
}

} // namespace qr
} // namespace El

#endif // ifndef EL_QR_CHOLESKY_HPP
