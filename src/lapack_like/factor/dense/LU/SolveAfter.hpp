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
void SolveAfter( Orientation orientation, const Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
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
  const AbstractDistMatrix<F>& A, AbstractDistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
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
( Orientation orientation, const Matrix<F>& A, 
  const Matrix<Int>& p, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( p.Height() != A.Height() )
            LogicError("A and p must be the same height");
    )
    if( orientation == NORMAL )
    {
        PermuteRows( B, p );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        InversePermuteRows( B, p );
    }
}

template<typename F> 
void SolveAfter
( Orientation orientation, const AbstractDistMatrix<F>& A, 
  const AbstractDistMatrix<Int>& p, AbstractDistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        AssertSameGrids( A, B, p );
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( A.Height() != p.Height() )
            LogicError("A and p must be the same height");
    )
    if( orientation == NORMAL )
    {
        PermuteRows( B, p );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        InversePermuteRows( B, p );
    }
}

template<typename F> 
void SolveAfter
( Orientation orientation, const Matrix<F>& A, 
  const Matrix<Int>& p, 
  const Matrix<Int>& q, 
        Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( p.Height() != A.Height() )
            LogicError("A and p must be the same height");
        if( q.Height() != A.Height() )
            LogicError("A and q must be the same height");
    )
    if( orientation == NORMAL )
    {
        PermuteRows( B, p );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        InversePermuteRows( B, q );
    }
    else
    {
        PermuteRows( B, q );
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        InversePermuteRows( B, p );
    }
}

template<typename F> 
void SolveAfter
( Orientation orientation, const AbstractDistMatrix<F>& A, 
  const AbstractDistMatrix<Int>& p, const AbstractDistMatrix<Int>& q,
        AbstractDistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        AssertSameGrids( A, B, p, q );
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( A.Height() != p.Height() )
            LogicError("A and p must be the same height");
        if( A.Height() != q.Height() )
            LogicError("A and q must be the same height");
    )
    if( orientation == NORMAL )
    {
        PermuteRows( B, p );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        InversePermuteRows( B, q );
    }
    else
    {
        PermuteRows( B, q );
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        InversePermuteRows( B, p );
    }
}

} // namespace lu
} // namespace El

#endif // ifndef EL_LU_SOLVEAFTER_HPP
