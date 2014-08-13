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
    if( A.Height() > A.Width() )
        LogicError("A A^H will be singular");
    Herk( UPPER, NORMAL, F(1), A, R );
    El::ReverseCholesky( UPPER, R );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R, A );
}

template<typename F> 
void Cholesky( AbstractDistMatrix<F>& APre, AbstractDistMatrix<F>& RPre )
{
    DEBUG_ONLY(CallStackEntry cse("rq::Cholesky"))
    const Int m = APre.Height();
    const Int n = APre.Width();
    if( m > n )
        LogicError("A A^H will be singular");

    // Proxies cannot be resized since they might be views
    RPre.Resize( m, m );

    const Grid& g = APre.Grid();
    DistMatrix<F,VR,STAR> A(g);
    DistMatrix<F,STAR,STAR> R(g);
    Copy( APre, A, READ_WRITE_PROXY );
    Copy( RPre, R, WRITE_PROXY );

    Zero( R );
    Herk( UPPER, NORMAL, F(1), A.Matrix(), F(0), R.Matrix() );
    R.SumOver( A.RowComm() );
    El::ReverseCholesky( UPPER, R.Matrix() );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R.Matrix(), A.Matrix() );

    Copy( A, APre, RESTORE_READ_WRITE_PROXY );
    Copy( R, RPre, RESTORE_WRITE_PROXY );
}

} // namespace rq
} // namespace El

#endif // ifndef EL_RQ_CHOLESKY_HPP
