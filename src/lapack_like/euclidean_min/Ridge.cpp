/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

template<typename Field>
void Ridge
( Orientation orientation,
  const Matrix<Field>& A,
  const Matrix<Field>& B,
        Base<Field> gamma,
        Matrix<Field>& X,
  RidgeAlg alg )
{
    EL_DEBUG_CSE

    const bool normal = ( orientation==NORMAL );
    const Int m = ( normal ? A.Height() : A.Width()  );
    const Int n = ( normal ? A.Width()  : A.Height() );
    if( orientation == TRANSPOSE && IsComplex<Field>::value )
        LogicError("Transpose version of complex Ridge not yet supported");

    if( m >= n )
    {
        Matrix<Field> Z;
        if( alg == RIDGE_CHOLESKY )
        {
            if( orientation == NORMAL )
                Herk( LOWER, ADJOINT, Base<Field>(1), A, Z );
            else
                Herk( LOWER, NORMAL, Base<Field>(1), A, Z );
            ShiftDiagonal( Z, Field(gamma*gamma) );
            Cholesky( LOWER, Z );
            if( orientation == NORMAL )
                Gemm( ADJOINT, NORMAL, Field(1), A, B, X );
            else
                Gemm( NORMAL, NORMAL, Field(1), A, B, X );
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
            FillDiagonal( ZB, Field(gamma) );
            // NOTE: This QR factorization could exploit the upper-triangular
            //       structure of the diagonal matrix ZB
            qr::ExplicitTriang( Z );
            if( orientation == NORMAL )
                Gemm( ADJOINT, NORMAL, Field(1), A, B, X );
            else
                Gemm( NORMAL, NORMAL, Field(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else
        {
            Matrix<Field> U, V;
            Matrix<Base<Field>> s;
            if( orientation == NORMAL )
            {
                SVDCtrl<Base<Field>> ctrl;
                ctrl.overwrite = false;
                SVD( A, U, s, V, ctrl );
            }
            else
            {
                Matrix<Field> AAdj;
                Adjoint( A, AAdj );

                SVDCtrl<Base<Field>> ctrl;
                ctrl.overwrite = true;
                SVD( AAdj, U, s, V, ctrl );
            }
            auto sigmaMap =
              [=]( const Base<Field>& sigma )
              { return sigma / (sigma*sigma + gamma*gamma); };
            EntrywiseMap( s, MakeFunction(sigmaMap) );
            Gemm( ADJOINT, NORMAL, Field(1), U, B, X );
            DiagonalScale( LEFT, NORMAL, s, X );
            U = X;
            Gemm( NORMAL, NORMAL, Field(1), V, U, X );
        }
    }
    else
    {
        LogicError("This case not yet supported");
    }
}

template<typename Field>
void Ridge
( Orientation orientation,
  const AbstractDistMatrix<Field>& APre,
  const AbstractDistMatrix<Field>& BPre,
        Base<Field> gamma,
        AbstractDistMatrix<Field>& XPre,
        RidgeAlg alg )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Field,Field,MC,MR>
      AProx( APre ),
      BProx( BPre );
    DistMatrixWriteProxy<Field,Field,MC,MR>
      XProx( XPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& X = XProx.Get();

    const bool normal = ( orientation==NORMAL );
    const Int m = ( normal ? A.Height() : A.Width()  );
    const Int n = ( normal ? A.Width()  : A.Height() );
    if( orientation == TRANSPOSE && IsComplex<Field>::value )
        LogicError("Transpose version of complex Ridge not yet supported");

    if( m >= n )
    {
        DistMatrix<Field> Z(A.Grid());
        if( alg == RIDGE_CHOLESKY )
        {
            if( orientation == NORMAL )
                Herk( LOWER, ADJOINT, Base<Field>(1), A, Z );
            else
                Herk( LOWER, NORMAL, Base<Field>(1), A, Z );
            ShiftDiagonal( Z, Field(gamma*gamma) );
            Cholesky( LOWER, Z );
            if( orientation == NORMAL )
                Gemm( ADJOINT, NORMAL, Field(1), A, B, X );
            else
                Gemm( NORMAL, NORMAL, Field(1), A, B, X );
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
            FillDiagonal( ZB, Field(gamma) );
            // NOTE: This QR factorization could exploit the upper-triangular
            //       structure of the diagonal matrix ZB
            qr::ExplicitTriang( Z );
            if( orientation == NORMAL )
                Gemm( ADJOINT, NORMAL, Field(1), A, B, X );
            else
                Gemm( NORMAL, NORMAL, Field(1), A, B, X );
            cholesky::SolveAfter( LOWER, NORMAL, Z, X );
        }
        else
        {
            DistMatrix<Field> U(A.Grid()), V(A.Grid());
            DistMatrix<Base<Field>,VR,STAR> s(A.Grid());
            if( orientation == NORMAL )
            {
                SVDCtrl<Base<Field>> ctrl;
                ctrl.overwrite = false;
                SVD( A, U, s, V, ctrl );
            }
            else
            {
                DistMatrix<Field> AAdj(A.Grid());
                Adjoint( A, AAdj );

                SVDCtrl<Base<Field>> ctrl;
                ctrl.overwrite = true;
                SVD( AAdj, U, s, V );
            }

            auto sigmaMap =
              [=]( const Base<Field>& sigma )
              { return sigma / (sigma*sigma + gamma*gamma); };
            EntrywiseMap( s, MakeFunction(sigmaMap) );
            Gemm( ADJOINT, NORMAL, Field(1), U, B, X );
            DiagonalScale( LEFT, NORMAL, s, X );
            U = X;
            Gemm( NORMAL, NORMAL, Field(1), V, U, X );
        }
    }
    else
    {
        LogicError("This case not yet supported");
    }
}

template<typename Field>
void Ridge
( Orientation orientation,
  const SparseMatrix<Field>& A,
  const Matrix<Field>& B,
        Base<Field> gamma,
        Matrix<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
    )

    const Int n = A.Width();
    SparseMatrix<Field> G;
    Zeros( G, n, n );
    ShiftDiagonal( G, gamma );

    Tikhonov( orientation, A, B, G, X, ctrl );
}

template<typename Field>
void Ridge
( Orientation orientation,
  const DistSparseMatrix<Field>& A,
  const DistMultiVec<Field>& B,
        Base<Field> gamma,
        DistMultiVec<Field>& X,
  const LeastSquaresCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
    )

    const Int n = A.Width();
    DistSparseMatrix<Field> G(A.Grid());
    Zeros( G, n, n );
    ShiftDiagonal( G, gamma );

    Tikhonov( orientation, A, B, G, X, ctrl );
}

#define PROTO(Field) \
  template void Ridge \
  ( Orientation orientation, \
    const Matrix<Field>& A, \
    const Matrix<Field>& B, \
          Base<Field> gamma, \
          Matrix<Field>& X, \
          RidgeAlg alg ); \
  template void Ridge \
  ( Orientation orientation, \
    const AbstractDistMatrix<Field>& A, \
    const AbstractDistMatrix<Field>& B, \
          Base<Field> gamma, \
          AbstractDistMatrix<Field>& X, \
          RidgeAlg alg ); \
  template void Ridge \
  ( Orientation orientation, \
    const SparseMatrix<Field>& A, \
    const Matrix<Field>& B, \
          Base<Field> gamma, \
          Matrix<Field>& X, \
    const LeastSquaresCtrl<Base<Field>>& ctrl ); \
  template void Ridge \
  ( Orientation orientation, \
    const DistSparseMatrix<Field>& A, \
    const DistMultiVec<Field>& B, \
          Base<Field> gamma, \
          DistMultiVec<Field>& X, \
    const LeastSquaresCtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
