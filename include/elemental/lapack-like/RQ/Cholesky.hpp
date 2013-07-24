/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_CHOLESKY_RQ_HPP
#define ELEM_LAPACK_CHOLESKY_RQ_HPP

#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/lapack-like/Cholesky.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace rq {

// NOTE: This version is designed for short-fat matrices and is much less
//       numerically stable than Householder-based RQ factorizations
//

template<typename F> 
inline void
Cholesky( Matrix<F>& A, Matrix<F>& R )
{
#ifndef RELEASE
    CallStackEntry entry("rq::Cholesky");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( height > width )
        throw std::logic_error("A A^H will be singular");
    Herk( UPPER, NORMAL, F(1), A, R );
    elem::ReverseCholesky( UPPER, R );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R, A );
}

template<typename F> 
inline void
Cholesky( DistMatrix<F,STAR,VR>& A, DistMatrix<F,STAR,STAR>& R )
{
#ifndef RELEASE
    CallStackEntry entry("rq::Cholesky");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( height > width )
        throw std::logic_error("A A^H will be singular");
    Zeros( R, height, height );
    Herk( UPPER, NORMAL, F(1), A.Matrix(), F(0), R.Matrix() );
    R.SumOverGrid();
    elem::ReverseCholesky( UPPER, R.Matrix() );
    Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), R.Matrix(), A.Matrix() );
}

} // namespace rq
} // namespace elem

#endif // ifndef ELEM_LAPACK_RQ_CHOLESKY_HPP
