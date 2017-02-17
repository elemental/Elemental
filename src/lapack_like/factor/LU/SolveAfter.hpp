/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LU_SOLVEAFTER_HPP
#define EL_LU_SOLVEAFTER_HPP

namespace El {
namespace lu {

template<typename F>
void SolveAfter
( Orientation orientation,
  const Matrix<F>& A,
        Matrix<F>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    if( orientation == NORMAL )
    {
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    }
}

template<typename F>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, B );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    if( orientation == NORMAL )
    {
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    }
}

template<typename F>
void SolveAfter
( Orientation orientation,
  const Matrix<F>& A,
  const Permutation& P,
        Matrix<F>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    if( orientation == NORMAL )
    {
        P.PermuteRows( B );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        P.InversePermuteRows( B );
    }
}

template<typename F>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<F>& A,
  const DistPermutation& P,
        AbstractDistMatrix<F>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, B );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    if( orientation == NORMAL )
    {
        P.PermuteRows( B );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        P.InversePermuteRows( B );
    }
}

template<typename F>
void SolveAfter
( Orientation orientation,
  const Matrix<F>& A,
  const Permutation& P,
  const Permutation& Q,
        Matrix<F>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    if( orientation == NORMAL )
    {
        P.PermuteRows( B );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        Q.InversePermuteRows( B );
    }
    else
    {
        Q.PermuteRows( B );
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        P.InversePermuteRows( B );
    }
}

template<typename F>
void SolveAfter
( Orientation orientation,
  const AbstractDistMatrix<F>& A,
  const DistPermutation& P,
  const DistPermutation& Q,
        AbstractDistMatrix<F>& B )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, B );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    if( orientation == NORMAL )
    {
        P.PermuteRows( B );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        Q.InversePermuteRows( B );
    }
    else
    {
        Q.PermuteRows( B );
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        P.InversePermuteRows( B );
    }
}

} // namespace lu
} // namespace El

#endif // ifndef EL_LU_SOLVEAFTER_HPP
