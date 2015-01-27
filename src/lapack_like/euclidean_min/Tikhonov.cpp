/*
   Copyright (c) 2009-2015, Jack Poulson
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
            Herk( LOWER, ADJOINT, Base<F>(1), A, Z );
            Herk( LOWER, ADJOINT, Base<F>(1), Gamma, Base<F>(1), Z );
            Cholesky( LOWER, Z );
        }
        else
        {
            const Int mGamma = Gamma.Height();
            Zeros( Z, m+mGamma, n );
            auto ZT = Z( IR(0,m),        IR(0,n) );
            auto ZB = Z( IR(m,m+mGamma), IR(0,n) );
            ZT = A;
            ZB = Gamma;
            qr::ExplicitTriang( Z ); 
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
( const AbstractDistMatrix<F>& APre, const AbstractDistMatrix<F>& BPre, 
  const AbstractDistMatrix<F>& Gamma, AbstractDistMatrix<F>& XPre, 
  TikhonovAlg alg )
{
    DEBUG_ONLY(CallStackEntry cse("Tikhonov"))

    auto APtr = ReadProxy<F,MC,MR>( &APre );  auto& A = *APtr;
    auto BPtr = ReadProxy<F,MC,MR>( &BPre );  auto& B = *BPtr;
    auto XPtr = WriteProxy<F,MC,MR>( &XPre ); auto& X = *XPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    if( Gamma.Width() != n )
        LogicError("Tikhonov matrix was the wrong width");

    if( m >= n )
    {
        DistMatrix<F> Z(A.Grid());
        if( alg == TIKHONOV_CHOLESKY )
        {
            Herk( LOWER, ADJOINT, Base<F>(1), A, Z );
            Herk( LOWER, ADJOINT, Base<F>(1), Gamma, Base<F>(1), Z );
            Cholesky( LOWER, Z );
        }
        else
        {
            const Int mGamma = Gamma.Height();
            Zeros( Z, m+mGamma, n );
            auto ZT = Z( IR(0,m),        IR(0,n) );
            auto ZB = Z( IR(m,m+mGamma), IR(0,n) );
            ZT = A;
            ZB = Gamma;
            qr::ExplicitTriang( Z ); 
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
( const DistSparseMatrix<F>& A, const DistMultiVec<F>& B,
  const DistSparseMatrix<F>& Gamma, DistMultiVec<F>& X, const BisectCtrl& ctrl )
{
    DEBUG_ONLY(
        CallStackEntry cse("Tikhonov");
        if( A.Height() != B.Height() )
            LogicError("Heights of A and B must match");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    DistSparseMatrix<F> C(A.Comm());

    X.SetComm( B.Comm() );
    Zeros( X, n, B.Width() );
    if( m >= n )
    {
        Herk( LOWER, ADJOINT, Base<F>(1), A, C );
        Herk( LOWER, ADJOINT, Base<F>(1), Gamma, Base<F>(1), C );
        MakeHermitian( LOWER, C );

        Multiply( ADJOINT, F(1), A, B, F(0), X );
        HermitianSolve( C, X, ctrl );
    }
    else
    {
        Herk( LOWER, NORMAL, Base<F>(1), A, C );
        Herk( LOWER, NORMAL, Base<F>(1), Gamma, Base<F>(1), C );
        MakeHermitian( LOWER, C );

        DistMultiVec<F> BCopy(B.Comm());
        BCopy = B;
        HermitianSolve( C, BCopy, ctrl );
        Multiply( ADJOINT, F(1), A, BCopy, F(0), X );
    }
}

#define PROTO(F) \
  template void Tikhonov \
  ( const Matrix<F>& A, const Matrix<F>& B, \
    const Matrix<F>& Gamma, Matrix<F>& X, TikhonovAlg alg ); \
  template void Tikhonov \
  ( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, \
    const AbstractDistMatrix<F>& Gamma, AbstractDistMatrix<F>& X, \
    TikhonovAlg alg ); \
  template void Tikhonov \
  ( const DistSparseMatrix<F>& A, const DistMultiVec<F>& B, \
    const DistSparseMatrix<F>& Gamma, DistMultiVec<F>& X, \
    const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
