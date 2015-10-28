/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
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
    DEBUG_ONLY(
      CSE cse("lu::SolveAfter");
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
  const ElementalMatrix<F>& A,
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("lu::SolveAfter");
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
  const Matrix<Int>& rowPiv,
        Matrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("lu::SolveAfter");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
      if( rowPiv.Height() != A.Height() )
          LogicError("A and p must be the same height");
    )
    if( orientation == NORMAL )
    {
        ApplyRowPivots( B, rowPiv );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        ApplyInverseRowPivots( B, rowPiv );
    }
}

template<typename F> 
void SolveAfter
( Orientation orientation,
  const ElementalMatrix<F>& A, 
  const ElementalMatrix<Int>& rowPiv,
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("lu::SolveAfter");
      AssertSameGrids( A, B, rowPiv );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
      if( A.Height() != rowPiv.Height() )
          LogicError("A and rowPiv must be the same height");
    )
    if( orientation == NORMAL )
    {
        ApplyRowPivots( B, rowPiv );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        ApplyInverseRowPivots( B, rowPiv );
    }
}

template<typename F> 
void SolveAfter
( Orientation orientation,
  const Matrix<F>& A, 
  const Matrix<Int>& rowPiv, 
  const Matrix<Int>& colPiv, 
        Matrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("lu::SolveAfter");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
      if( rowPiv.Height() != A.Height() )
          LogicError("A and rowPiv must be the same height");
      if( colPiv.Height() != A.Height() )
          LogicError("A and colPiv must be the same height");
    )
    if( orientation == NORMAL )
    {
        ApplyRowPivots( B, rowPiv );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        ApplyInverseRowPivots( B, colPiv );
    }
    else
    {
        ApplyRowPivots( B, colPiv );
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        ApplyInverseRowPivots( B, rowPiv );
    }
}

template<typename F> 
void SolveAfter
( Orientation orientation,
  const ElementalMatrix<F>& A, 
  const ElementalMatrix<Int>& rowPiv,
  const ElementalMatrix<Int>& colPiv,
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(
      CSE cse("lu::SolveAfter");
      AssertSameGrids( A, B, rowPiv, colPiv );
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
      if( A.Height() != rowPiv.Height() )
          LogicError("A and rowPiv must be the same height");
      if( A.Height() != colPiv.Height() )
          LogicError("A and colPiv must be the same height");
    )
    if( orientation == NORMAL )
    {
        ApplyRowPivots( B, rowPiv );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        ApplyInverseRowPivots( B, colPiv );
    }
    else
    {
        ApplyRowPivots( B, colPiv );
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        ApplyInverseRowPivots( B, rowPiv );
    }
}

} // namespace lu
} // namespace El

#endif // ifndef EL_LU_SOLVEAFTER_HPP
