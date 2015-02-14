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
            Herk( LOWER, ADJOINT, Base<F>(1), A, Z );
            ShiftDiagonal( Z, F(alpha*alpha) );
            Cholesky( LOWER, Z );
            Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else if( alg == RIDGE_QR )
        {
            Zeros( Z, m+n, n );
            auto ZT = Z( IR(0,m),   IR(0,n) );
            auto ZB = Z( IR(m,m+n), IR(0,n) );
            ZT = A;
            FillDiagonal( ZB, F(alpha*alpha) );
            // NOTE: This QR factorization could exploit the upper-triangular
            //       structure of the diagonal matrix ZB
            qr::ExplicitTriang( Z );
            Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else
        {
            Matrix<F> U, V;
            Matrix<Base<F>> s; 
            U = A;
            SVD( U, s, V );
            auto sigmaMap = 
              [=]( Base<F> sigma ) 
              { return sigma / (sigma*sigma + alpha*alpha); };
            EntrywiseMap( s, function<Base<F>(Base<F>)>(sigmaMap) );
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
( const AbstractDistMatrix<F>& APre, const AbstractDistMatrix<F>& BPre, 
  Base<F> alpha, AbstractDistMatrix<F>& XPre, RidgeAlg alg )
{
    DEBUG_ONLY(CallStackEntry cse("Ridge"))

    auto APtr = ReadProxy<F,MC,MR>( &APre );  auto& A = *APtr;
    auto BPtr = ReadProxy<F,MC,MR>( &BPre );  auto& B = *BPtr;
    auto XPtr = WriteProxy<F,MC,MR>(& XPre ); auto& X = *XPtr;

    const Int m = A.Height();
    const Int n = A.Width();

    if( m >= n )
    {
        DistMatrix<F> Z(A.Grid());
        if( alg == RIDGE_CHOLESKY )
        {
            Herk( LOWER, ADJOINT, Base<F>(1), A, Z );
            ShiftDiagonal( Z, F(alpha*alpha) );
            Cholesky( LOWER, Z );
            Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else if( alg == RIDGE_QR )
        {
            Zeros( Z, m+n, n );
            auto ZT = Z( IR(0,m),   IR(0,n) ); 
            auto ZB = Z( IR(m,m+n), IR(0,n) );
            ZT = A;
            FillDiagonal( ZB, F(alpha*alpha) );
            // NOTE: This QR factorization could exploit the upper-triangular
            //       structure of the diagonal matrix ZB
            qr::ExplicitTriang( Z );
            Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else
        {
            DistMatrix<F> U(A.Grid()), V(A.Grid());
            DistMatrix<Base<F>,VR,STAR> s(A.Grid());
            U = A;
            SVD( U, s, V );
            auto sigmaMap = 
              [=]( Base<F> sigma ) 
              { return sigma / (sigma*sigma + alpha*alpha); };
            EntrywiseMap( s, function<Base<F>(Base<F>)>(sigmaMap) );
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
( const DistSparseMatrix<F>& A, const DistMultiVec<F>& B, Base<F> alpha,
        DistMultiVec<F>& X, const BisectCtrl& ctrl )
{
    DEBUG_ONLY(
        CallStackEntry cse("Ridge");
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
        ShiftDiagonal( C, F(alpha*alpha) );
        MakeHermitian( LOWER, C );

        Multiply( ADJOINT, F(1), A, B, F(0), X );
        HermitianSolve( C, X, ctrl );
    }
    else
    {
        Herk( LOWER, NORMAL, Base<F>(1), A, C );
        ShiftDiagonal( C, F(alpha*alpha) );
        MakeHermitian( LOWER, C );

        DistMultiVec<F> BCopy(B.Comm());
        BCopy = B;
        HermitianSolve( C, BCopy, ctrl );
        Multiply( ADJOINT, F(1), A, BCopy, F(0), X );
    }
}

#define PROTO(F) \
  template void Ridge \
  ( const Matrix<F>& A, const Matrix<F>& B, \
    Base<F> alpha, Matrix<F>& X, RidgeAlg alg ); \
  template void Ridge \
  ( const AbstractDistMatrix<F>& A, const AbstractDistMatrix<F>& B, \
    Base<F> alpha, AbstractDistMatrix<F>& X, RidgeAlg alg ); \
  template void Ridge \
  ( const DistSparseMatrix<F>& A, const DistMultiVec<F>& B, Base<F> alpha, \
    DistMultiVec<F>& X, const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
