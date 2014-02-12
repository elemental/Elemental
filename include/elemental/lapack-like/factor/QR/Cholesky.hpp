/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CHOLESKY_QR_HPP
#define ELEM_CHOLESKY_QR_HPP

#include ELEM_HERK_INC
#include ELEM_TRSM_INC

#include ELEM_CHOLESKY_INC

#include ELEM_ZEROS_INC

namespace elem {
namespace qr {

// NOTE: This version is designed for tall-skinny matrices and is much less
//       numerically stable than Householder-based QR factorizations
//
// Computes the QR factorization of full-rank tall-skinny matrix A and 
// overwrites A with Q
//

template<typename F> 
inline void
Cholesky( Matrix<F>& A, Matrix<F>& R )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Cholesky"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( height < width )
        LogicError("A^H A will be singular");
    Herk( UPPER, ADJOINT, F(1), A, R );
    elem::Cholesky( UPPER, R );
    Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), R, A );
}

template<typename F> 
inline void
Cholesky( DistMatrix<F,VC,STAR>& A, DistMatrix<F,STAR,STAR>& R )
{
    DEBUG_ONLY(CallStackEntry cse("qr::Cholesky"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( height < width )
        LogicError("A^H A will be singular");
    Zeros( R, width, width );
    Herk( UPPER, ADJOINT, F(1), A.Matrix(), F(0), R.Matrix() );
    R.SumOver( A.ColComm() );
    elem::Cholesky( UPPER, R.Matrix() );
    Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), R.Matrix(), A.Matrix() );
}

} // namespace qr
} // namespace elem

#endif // ifndef ELEM_QR_CHOLESKY_HPP
