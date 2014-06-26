/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El-lite.hpp"

namespace El {

template<typename F> 
void Ridge
( const Matrix<F>& A, const Matrix<F>& B, 
  Base<F> alpha, Matrix<F>& X, RidgeAlg alg )
{
    DEBUG_ONLY(CallStackEntry cse("Ridge"))
    const Int m = A.Height();
    const Int n = A.Width();

    if( m >= n )
    {
        Matrix<F> Z;
        if( alg == RIDGE_CHOLESKY )
        {
            Herk( LOWER, ADJOINT, F(1), A, Z );
            UpdateDiagonal( Z, F(alpha*alpha) );
            Cholesky( LOWER, Z );
            Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else if( alg == RIDGE_QR )
        {
            Zeros( Z, m+n, n );
            auto ZT = View( Z, 0, 0, m, n );
            auto ZB = View( Z, m, 0, n, n );
            ZT = A;
            SetDiagonal( ZB, F(alpha*alpha) );
            // NOTE: This QR factorization could exploit the upper-triangular
            //       structure of the diagonal matrix ZB
            QR( Z );
            Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else
        {
            Matrix<F> U, V;
            Matrix<Base<F>> s; 
            U = A;
            SVD( U, s, V );
            EntrywiseMap
            ( s, [=]( Base<F> sigma ) 
                 { return sigma / (sigma*sigma + alpha*alpha); } );
            Gemm( ADJOINT, NORMAL, F(1), U, B, X );
            DiagonalScale( LEFT, NORMAL, s, X );
            U = X;
            Gemm( NORMAL, NORMAL, F(1), V, U, X );
        }
    }
    else
    {
        LogicError("This case not yet supported");
    }
}

template<typename F> 
void Ridge
( const DistMatrix<F>& A, const DistMatrix<F>& B, 
  Base<F> alpha, DistMatrix<F>& X, RidgeAlg alg )
{
    DEBUG_ONLY(CallStackEntry cse("Ridge"))
    const Int m = A.Height();
    const Int n = A.Width();

    if( m >= n )
    {
        DistMatrix<F> Z(A.Grid());
        if( alg == RIDGE_CHOLESKY )
        {
            Herk( LOWER, ADJOINT, F(1), A, Z );
            UpdateDiagonal( Z, F(alpha*alpha) );
            Cholesky( LOWER, Z );
            Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else if( alg == RIDGE_QR )
        {
            Zeros( Z, m+n, n );
            auto ZT = View( Z, 0, 0, m, n );
            auto ZB = View( Z, m, 0, n, n );
            ZT = A;
            SetDiagonal( ZB, F(alpha*alpha) );
            // NOTE: This QR factorization could exploit the upper-triangular
            //       structure of the diagonal matrix ZB
            QR( Z );
            Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else
        {
            DistMatrix<F> U(A.Grid()), V(A.Grid());
            DistMatrix<Base<F>,VR,STAR> s(A.Grid());
            U = A;
            SVD( U, s, V );
            EntrywiseMap
            ( s, [=]( Base<F> sigma )
                 { return sigma / (sigma*sigma + alpha*alpha); } );
            Gemm( ADJOINT, NORMAL, F(1), U, B, X );
            DiagonalScale( LEFT, NORMAL, s, X );
            U = X;
            Gemm( NORMAL, NORMAL, F(1), V, U, X );
        }
    }
    else
    {
        LogicError("This case not yet supported");
    }
}

#define PROTO(F) \
  template void Ridge \
  ( const Matrix<F>& A, const Matrix<F>& B, \
    Base<F> alpha, Matrix<F>& X, RidgeAlg alg ); \
  template void Ridge \
  ( const DistMatrix<F>& A, const DistMatrix<F>& B, \
    Base<F> alpha, DistMatrix<F>& X, RidgeAlg alg );

PROTO(float)
PROTO(double)
PROTO(Complex<float>)
PROTO(Complex<double>)

} // namespace El
