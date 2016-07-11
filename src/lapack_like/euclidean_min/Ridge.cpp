/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename F> 
void Ridge
( Orientation orientation,
  const Matrix<F>& A,
  const Matrix<F>& B, 
        Base<F> gamma,
        Matrix<F>& X, 
  RidgeAlg alg )
{
    DEBUG_CSE

    const bool normal = ( orientation==NORMAL );
    const Int m = ( normal ? A.Height() : A.Width()  );
    const Int n = ( normal ? A.Width()  : A.Height() );
    if( orientation == TRANSPOSE && IsComplex<F>::value )
        LogicError("Transpose version of complex Ridge not yet supported");

    if( m >= n )
    {
        Matrix<F> Z;
        if( alg == RIDGE_CHOLESKY )
        {
            if( orientation == NORMAL )
                Herk( LOWER, ADJOINT, Base<F>(1), A, Z );
            else
                Herk( LOWER, NORMAL, Base<F>(1), A, Z );
            ShiftDiagonal( Z, F(gamma*gamma) );
            Cholesky( LOWER, Z );
            if( orientation == NORMAL )
                Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            else
                Gemm( NORMAL, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else if( alg == RIDGE_QR )
        {
            Zeros( Z, m+n, n );
            auto ZT = Z( IR(0,m),   IR(0,n) );
            auto ZB = Z( IR(m,m+n), IR(0,n) );
            if( orientation == NORMAL )
                ZT = A;
            else
                Adjoint( A, ZT );
            FillDiagonal( ZB, F(gamma*gamma) );
            // NOTE: This QR factorization could exploit the upper-triangular
            //       structure of the diagonal matrix ZB
            qr::ExplicitTriang( Z );
            if( orientation == NORMAL )
                Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            else
                Gemm( NORMAL, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else
        {
            Matrix<F> U, V;
            Matrix<Base<F>> s; 
            if( orientation == NORMAL )
            {
                SVDCtrl<Base<F>> ctrl;
                ctrl.overwrite = false;
                SVD( A, U, s, V, ctrl );
            }
            else
            {
                Matrix<F> AAdj;
                Adjoint( A, AAdj );

                SVDCtrl<Base<F>> ctrl;
                ctrl.overwrite = true;
                SVD( AAdj, U, s, V, ctrl );
            }
            auto sigmaMap = 
              [=]( Base<F> sigma ) 
              { return sigma / (sigma*sigma + gamma*gamma); };
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
( Orientation orientation,
  const ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& BPre, 
        Base<F> gamma,
        ElementalMatrix<F>& XPre, 
        RidgeAlg alg )
{
    DEBUG_CSE

    DistMatrixReadProxy<F,F,MC,MR>
      AProx( APre ),
      BProx( BPre );
    DistMatrixWriteProxy<F,F,MC,MR>
      XProx( XPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& X = XProx.Get();

    const bool normal = ( orientation==NORMAL );
    const Int m = ( normal ? A.Height() : A.Width()  );
    const Int n = ( normal ? A.Width()  : A.Height() );
    if( orientation == TRANSPOSE && IsComplex<F>::value )
        LogicError("Transpose version of complex Ridge not yet supported");

    if( m >= n )
    {
        DistMatrix<F> Z(A.Grid());
        if( alg == RIDGE_CHOLESKY )
        {
            if( orientation == NORMAL )
                Herk( LOWER, ADJOINT, Base<F>(1), A, Z );
            else
                Herk( LOWER, NORMAL, Base<F>(1), A, Z );
            ShiftDiagonal( Z, F(gamma*gamma) );
            Cholesky( LOWER, Z );
            if( orientation == NORMAL )
                Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            else
                Gemm( NORMAL, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else if( alg == RIDGE_QR )
        {
            Zeros( Z, m+n, n );
            auto ZT = Z( IR(0,m),   IR(0,n) ); 
            auto ZB = Z( IR(m,m+n), IR(0,n) );
            if( orientation == NORMAL )
                ZT = A;
            else
                Adjoint( A, ZT );
            FillDiagonal( ZB, F(gamma*gamma) );
            // NOTE: This QR factorization could exploit the upper-triangular
            //       structure of the diagonal matrix ZB
            qr::ExplicitTriang( Z );
            if( orientation == NORMAL )
                Gemm( ADJOINT, NORMAL, F(1), A, B, X );
            else
                Gemm( NORMAL, NORMAL, F(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else
        {
            DistMatrix<F> U(A.Grid()), V(A.Grid());
            DistMatrix<Base<F>,VR,STAR> s(A.Grid());
            if( orientation == NORMAL )
            {
                SVDCtrl<Base<F>> ctrl;
                ctrl.overwrite = false;
                SVD( A, U, s, V, ctrl );
            }
            else
            {
                DistMatrix<F> AAdj(A.Grid());
                Adjoint( A, AAdj );

                SVDCtrl<Base<F>> ctrl;
                ctrl.overwrite = true;
                SVD( AAdj, U, s, V );
            }

            auto sigmaMap = 
              [=]( Base<F> sigma ) 
              { return sigma / (sigma*sigma + gamma*gamma); };
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
( Orientation orientation,
  const SparseMatrix<F>& A,
  const Matrix<F>& B, 
        Base<F> gamma,
        Matrix<F>& X, 
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
    )

    const Int n = A.Width();
    SparseMatrix<F> G;
    Zeros( G, n, n );
    ShiftDiagonal( G, gamma );

    Tikhonov( orientation, A, B, G, X, ctrl );
}

template<typename F>
void Ridge
( Orientation orientation,
  const DistSparseMatrix<F>& A,
  const DistMultiVec<F>& B, 
        Base<F> gamma,
        DistMultiVec<F>& X, 
  const LeastSquaresCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
    )

    const Int n = A.Width();
    DistSparseMatrix<F> G(A.Comm());
    Zeros( G, n, n );
    ShiftDiagonal( G, gamma );

    Tikhonov( orientation, A, B, G, X, ctrl );
}

#define PROTO(F) \
  template void Ridge \
  ( Orientation orientation, \
    const Matrix<F>& A, \
    const Matrix<F>& B, \
          Base<F> gamma, \
          Matrix<F>& X, \
          RidgeAlg alg ); \
  template void Ridge \
  ( Orientation orientation, \
    const ElementalMatrix<F>& A, \
    const ElementalMatrix<F>& B, \
          Base<F> gamma, \
          ElementalMatrix<F>& X, \
          RidgeAlg alg ); \
  template void Ridge \
  ( Orientation orientation, \
    const SparseMatrix<F>& A, \
    const Matrix<F>& B, \
          Base<F> gamma, \
          Matrix<F>& X, \
    const LeastSquaresCtrl<Base<F>>& ctrl ); \
  template void Ridge \
  ( Orientation orientation, \
    const DistSparseMatrix<F>& A, \
    const DistMultiVec<F>& B, \
          Base<F> gamma, \
          DistMultiVec<F>& X, \
    const LeastSquaresCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
