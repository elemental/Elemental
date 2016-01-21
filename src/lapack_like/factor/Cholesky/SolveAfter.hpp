/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_CHOLESKY_SOLVEAFTER_HPP
#define EL_CHOLESKY_SOLVEAFTER_HPP

namespace El {
namespace cholesky {

template<typename F> 
void SolveAfter
( UpperOrLower uplo,
  Orientation orientation,
  const Matrix<F>& A,
        Matrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("cholesky::SolveAfter");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    if( orientation == TRANSPOSE )
        Conjugate( B );
    if( uplo == LOWER )
    {
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    if( orientation == TRANSPOSE )
        Conjugate( B );
}

template<typename F> 
void SolveAfter
( UpperOrLower uplo,
  Orientation orientation, 
  const Matrix<F>& A,
  const Permutation& P,
        Matrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("cholesky::SolveAfter");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    P.PermuteRows( B );
    if( orientation == TRANSPOSE )
        Conjugate( B );
    if( uplo == LOWER )
    {
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    if( orientation == TRANSPOSE )
        Conjugate( B );
    P.InversePermuteRows( B );
}

template<typename F> 
void SolveAfter
( UpperOrLower uplo,
  Orientation orientation, 
  const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("cholesky::SolveAfter");
      AssertSameGrids( A, B );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    if( orientation == TRANSPOSE )
        Conjugate( B );
    if( uplo == LOWER )
    {
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    if( orientation == TRANSPOSE )
        Conjugate( B );
}

template<typename F> 
void SolveAfter
( UpperOrLower uplo,
  Orientation orientation, 
  const AbstractDistMatrix<F>& A,
  const DistPermutation& P, 
        AbstractDistMatrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("cholesky::SolveAfter");
      AssertSameGrids( A, B );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    P.PermuteRows( B );
    if( orientation == TRANSPOSE )
        Conjugate( B );
    if( uplo == LOWER )
    {
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    if( orientation == TRANSPOSE )
        Conjugate( B );
    P.InversePermuteRows( B );
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_SOLVEAFTER_HPP
