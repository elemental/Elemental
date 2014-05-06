/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LU_SOLVEAFTER_HPP
#define ELEM_LU_SOLVEAFTER_HPP

#include ELEM_APPLYROWPIVOTS_INC
#include ELEM_TRSM_INC

namespace elem {
namespace lu {

template<typename F> 
inline void
SolveAfter( Orientation orientation, const Matrix<F>& A, Matrix<F>& B )
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
inline void
SolveAfter
( Orientation orientation, const DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Grid() != B.Grid() )
            LogicError("{A,B} must be distributed over the same grid");
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
inline void
SolveAfter
( Orientation orientation, 
  const Matrix<F>& A, const Matrix<Int>& pPerm, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( pPerm.Height() != A.Height() )
            LogicError("A and pPerm must be the same height");
    )
    if( orientation == NORMAL )
    {
        PermuteRows( B, pPerm );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        InversePermuteRows( B, pPerm );
    }
}

template<typename F,Dist UPerm> 
inline void
SolveAfter
( Orientation orientation, 
  const DistMatrix<F>& A, 
  const DistMatrix<Int,UPerm,STAR>& pPerm, DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Grid() != B.Grid() || A.Grid() != pPerm.Grid() )
            LogicError("{A,B,pPerm} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( A.Height() != pPerm.Height() )
            LogicError("A and pPerm must be the same height");
    )
    if( orientation == NORMAL )
    {
        PermuteRows( B, pPerm );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
    }
    else
    {
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        InversePermuteRows( B, pPerm );
    }
}

template<typename F> 
inline void
SolveAfter
( Orientation orientation, 
  const Matrix<F>& A, 
  const Matrix<Int>& pPerm, 
  const Matrix<Int>& qPerm, 
        Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( pPerm.Height() != A.Height() )
            LogicError("A and p must be the same height");
        if( qPerm.Height() != A.Height() )
            LogicError("A and q must be the same height");
    )
    if( orientation == NORMAL )
    {
        PermuteRows( B, pPerm );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        InversePermuteRows( B, qPerm );
    }
    else
    {
        PermuteRows( B, qPerm );
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        InversePermuteRows( B, pPerm );
    }
}

template<typename F,Dist UPerm> 
inline void
SolveAfter
( Orientation orientation, 
  const DistMatrix<F>& A, 
  const DistMatrix<Int,UPerm,STAR>& pPerm, 
  const DistMatrix<Int,UPerm,STAR>& qPerm,
        DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Grid() != B.Grid() || 
            A.Grid() != pPerm.Grid() || 
            pPerm.Grid() != qPerm.Grid() )
            LogicError("{A,B,pPerm,qPerm} must be distributed over same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( A.Height() != pPerm.Height() )
            LogicError("A and pPerm must be the same height");
        if( A.Height() != qPerm.Height() )
            LogicError("A and qPerm must be the same height");
    )
    if( orientation == NORMAL )
    {
        PermuteRows( B, pPerm );
        Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), A, B );
        InversePermuteRows( B, qPerm );
    }
    else
    {
        PermuteRows( B, qPerm );
        Trsm( LEFT, UPPER, orientation, NON_UNIT, F(1), A, B );
        Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
        InversePermuteRows( B, pPerm );
    }
}

} // namespace lu
} // namespace elem

#endif // ifndef ELEM_LU_SOLVEAFTER_HPP
