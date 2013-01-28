/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_CHOLESKY_QR_HPP
#define LAPACK_CHOLESKY_QR_HPP

#include "elemental/lapack-like/Cholesky.hpp"

namespace elem {

// NOTE: This version is designed for tall-skinny matrices and is much less
//       numerically stable than Householder-based QR factorizations
//
// Computes the QR factorization of full-rank tall-skinny matrix A and 
// overwrites A with Q
//

template<typename F> 
inline void
CholeskyQR( Matrix<F>& A, Matrix<F>& R )
{
#ifndef RELEASE
    PushCallStack("CholeskyQR");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( height < width )
        throw std::logic_error("A^H A will be singular");
    Zeros( width, width, R );
    Herk( UPPER, ADJOINT, F(1), A, F(0), R );
    Cholesky( UPPER, R );
    Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), R, A );
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename F> 
inline void
CholeskyQR( DistMatrix<F,VC,STAR>& A, DistMatrix<F,STAR,STAR>& R )
{
#ifndef RELEASE
    PushCallStack("CholeskyQR");
#endif
    const int height = A.Height();
    const int width = A.Width();
    if( height < width )
        throw std::logic_error("A^H A will be singular");
    Zeros( width, width, R );
    Herk( UPPER, ADJOINT, F(1), A.LocalMatrix(), F(0), R.LocalMatrix() );
    R.SumOverGrid();
    Cholesky( UPPER, R.LocalMatrix() );
    Trsm
    ( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), R.LocalMatrix(), A.LocalMatrix() );
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace elem

#endif // ifndef LAPACK_CHOLESKY_QR_HPP
