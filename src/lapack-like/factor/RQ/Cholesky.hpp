/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CHOLESKY_RQ_HPP
#define EL_CHOLESKY_RQ_HPP



namespace El {
namespace rq {

// NOTE: This version is designed for short-fat matrices and is much less
//       numerically stable than Householder-based RQ factorizations
//

template<typename F> 
void Cholesky( Matrix<F>& A, Matrix<F>& R )
{
    DEBUG_ONLY(CallStackEntry cse("rq::Cholesky"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( height > width )
        LogicError("A A^H will be singular");
    Herk( UPPER, NORMAL, F(1), A, R );
    El::ReverseCholesky( UPPER, R );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R, A );
}

template<typename F> 
void Cholesky( DistMatrix<F,STAR,VR>& A, DistMatrix<F,STAR,STAR>& R )
{
    DEBUG_ONLY(CallStackEntry cse("rq::Cholesky"))
    const Int height = A.Height();
    const Int width = A.Width();
    if( height > width )
        LogicError("A A^H will be singular");
    Zeros( R, height, height );
    Herk( UPPER, NORMAL, F(1), A.Matrix(), F(0), R.Matrix() );
    R.SumOver( A.RowComm() );
    El::ReverseCholesky( UPPER, R.Matrix() );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R.Matrix(), A.Matrix() );
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_CHOLESKY_HPP
