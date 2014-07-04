/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

template<typename F> 
void Tikhonov
( const Matrix<F>& A, const Matrix<F>& B, 
  const Matrix<F>& Gamma, Matrix<F>& X, TikhonovAlg alg )
{
    DEBUG_ONLY(CallStackEntry cse("Tikhonov"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( Gamma.Width() != n )
        LogicError("Tikhonov matrix was the wrong width");

    if( m >= n )
    {
        Matrix<F> Z;
        if( alg == TIKHONOV_CHOLESKY )
        {
            Herk( LOWER, ADJOINT, F(1), A, Z );
            Herk( LOWER, ADJOINT, F(1), Gamma, F(1), Z );
            Cholesky( LOWER, Z );
        }
        else
        {
            const Int mGamma = Gamma.Height();
            Zeros( Z, m+mGamma, n );
            auto ZT = View( Z, 0, 0, m,      n );
            auto ZB = View( Z, m, 0, mGamma, n );
            ZT = A;
            ZB = Gamma;
            QR( Z );
        }
        Gemm( ADJOINT, NORMAL, F(1), A, B, X );
        cholesky::SolveAfter( LOWER, NORMAL, Z, X );
    }
    else
    {
        LogicError("This case not yet supported");
    }
}

template<typename F> 
void Tikhonov
( const DistMatrix<F>& A, const DistMatrix<F>& B, 
  const DistMatrix<F>& Gamma, DistMatrix<F>& X, TikhonovAlg alg )
{
    DEBUG_ONLY(CallStackEntry cse("Tikhonov"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( Gamma.Width() != n )
        LogicError("Tikhonov matrix was the wrong width");

    if( m >= n )
    {
        DistMatrix<F> Z(A.Grid());
        if( alg == TIKHONOV_CHOLESKY )
        {
            Herk( LOWER, ADJOINT, F(1), A, Z );
            Herk( LOWER, ADJOINT, F(1), Gamma, F(1), Z );
            Cholesky( LOWER, Z );
        }
        else
        {
            const Int mGamma = Gamma.Height();
            Zeros( Z, m+mGamma, n );
            auto ZT = View( Z, 0, 0, m,      n );
            auto ZB = View( Z, m, 0, mGamma, n );
            ZT = A;
            ZB = Gamma;
            QR( Z );
        }
        Gemm( ADJOINT, NORMAL, F(1), A, B, X );
        cholesky::SolveAfter( LOWER, NORMAL, Z, X );
    }
    else
    {
        LogicError("This case not yet supported");
    }
}

#define PROTO(F) \
  template void Tikhonov \
  ( const Matrix<F>& A, const Matrix<F>& B, \
    const Matrix<F>& Gamma, Matrix<F>& X, TikhonovAlg alg ); \
  template void Tikhonov \
  ( const DistMatrix<F>& A, const DistMatrix<F>& B, \
    const DistMatrix<F>& Gamma, DistMatrix<F>& X, TikhonovAlg alg );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
