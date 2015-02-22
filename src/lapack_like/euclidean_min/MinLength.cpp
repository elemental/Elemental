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
void MinLength
( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, 
  Matrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MinLength"))
    const Int m = A.Height();
    const Int n = A.Width();
    if( m > n )
        LogicError("Assumed that height(A) <= width(A)");
    Matrix<F> t;
    Matrix<Base<F>> d;
    LQ( A, t, d );
    lq::SolveAfter( orientation, A, t, d, B, X );
}

template<typename F> 
void MinLength
( Orientation orientation, AbstractDistMatrix<F>& APre,
  const AbstractDistMatrix<F>& B, AbstractDistMatrix<F>& X )
{
    DEBUG_ONLY(CallStackEntry cse("MinLength"))

    auto APtr = ReadProxy<F,MC,MR>( &APre );
    auto& A = *APtr;

    const Int m = A.Height();
    const Int n = A.Width();
    if( m > n )
        LogicError("Assumed that height(A) <= width(A)");
    DistMatrix<F,MD,STAR> t(A.Grid());
    DistMatrix<Base<F>,MD,STAR> d(A.Grid());
    LQ( A, t, d );
    lq::SolveAfter( orientation, A, t, d, B, X );
}

template<typename F>
void MinLength
( Orientation orientation,
  const SparseMatrix<F>& A, const Matrix<F>& B, Matrix<F>& X,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(
      CallStackEntry cse("MinLength");
      if( orientation == NORMAL && A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
      if( orientation != NORMAL && A.Width() != B.Height() )
          LogicError("Width of A and height of B must match");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    if( m > n )
        LogicError("Assumed height(A) <= width(A)");

    // TODO: 
    // Given that the computational complexity of forming A^T A versus
    // forming A A^T depends upon the sparsity structure of A, rather than
    // just the matrix dimensions, the following branches should be exposed
    // as an option for the user until the time where the decision can be
    // automatically made in both an efficient and reliable manner.
    SparseMatrix<F> C;
    if( orientation == NORMAL )
    {
        Zeros( X, n, k );
        Herk( LOWER, NORMAL, Base<F>(1), A, C );
        MakeHermitian( LOWER, C );

        auto BCopy( B );
        HermitianSolve( C, BCopy, ctrl );
        Multiply( ADJOINT, F(1), A, BCopy, F(0), X );
    }
    else if( orientation == ADJOINT || !IsComplex<F>::val )
    {
        Zeros( X, m, k );
        Herk( LOWER, ADJOINT, Base<F>(1), A, C );
        MakeHermitian( LOWER, C );

        auto BCopy( B );
        HermitianSolve( C, BCopy, ctrl );
        Multiply( NORMAL, F(1), A, BCopy, F(0), X );
    }
    else
        LogicError("Complex transposed option not yet supported");
}

template<typename F>
void MinLength
( Orientation orientation,
  const DistSparseMatrix<F>& A, const DistMultiVec<F>& B, DistMultiVec<F>& X,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(
      CallStackEntry cse("MinLength");
      if( orientation == NORMAL && A.Height() != B.Height() )
          LogicError("Heights of A and B must match");
      if( orientation != NORMAL && A.Width() != B.Height() )
          LogicError("Width of A and height of B must match");
    )
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = B.Width();
    if( m > n )
        LogicError("Assumed height(A) <= width(A)");

    // TODO: 
    // Given that the computational complexity of forming A^T A versus
    // forming A A^T depends upon the sparsity structure of A, rather than
    // just the matrix dimensions, the following branches should be exposed
    // as an option for the user until the time where the decision can be
    // automatically made in both an efficient and reliable manner.
    DistSparseMatrix<F> C(A.Comm());
    X.SetComm( B.Comm() );
    if( orientation == NORMAL )
    {
        Zeros( X, n, k );
        Herk( LOWER, NORMAL, Base<F>(1), A, C );
        MakeHermitian( LOWER, C );

        DistMultiVec<F> BCopy(B.Comm());
        BCopy = B;
        HermitianSolve( C, BCopy, ctrl );
        Multiply( ADJOINT, F(1), A, BCopy, F(0), X );
    }
    else if( orientation == ADJOINT || !IsComplex<F>::val )
    {
        Zeros( X, m, k );
        Herk( LOWER, ADJOINT, Base<F>(1), A, C );
        MakeHermitian( LOWER, C );

        DistMultiVec<F> BCopy(B.Comm());
        BCopy = B;
        HermitianSolve( C, BCopy, ctrl );
        Multiply( NORMAL, F(1), A, BCopy, F(0), X );
    }
    else
        LogicError("Complex transposed option not yet supported");
}

#define PROTO(F) \
  template void MinLength \
  ( Orientation orientation, Matrix<F>& A, const Matrix<F>& B, \
    Matrix<F>& X ); \
  template void MinLength \
  ( Orientation orientation, AbstractDistMatrix<F>& A, \
    const AbstractDistMatrix<F>& B, AbstractDistMatrix<F>& X ); \
  template void MinLength \
  ( Orientation orientation, \
    const SparseMatrix<F>& A, const Matrix<F>& B, \
    Matrix<F>& X, const BisectCtrl& ctrl ); \
  template void MinLength \
  ( Orientation orientation, \
    const DistSparseMatrix<F>& A, const DistMultiVec<F>& B, \
    DistMultiVec<F>& X, const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace El
