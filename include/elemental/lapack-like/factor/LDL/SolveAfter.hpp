/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LDL_SOLVEAFTER_HPP
#define ELEM_LDL_SOLVEAFTER_HPP

#include ELEM_DIAGONALSOLVE_INC
#include ELEM_QUASIDIAGONALSOLVE_INC
#include ELEM_APPLYROWPIVOTS_INC
#include ELEM_TRSM_INC

namespace elem {
namespace ldl {

template<typename F> 
inline void
SolveAfter( const Matrix<F>& A, Matrix<F>& B, bool conjugated=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::SolveAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const bool checkIfSingular = false;
    const auto d = A.GetDiagonal();
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    DiagonalSolve( LEFT, NORMAL, d, B, checkIfSingular );
    Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
}

template<typename F> 
inline void
SolveAfter( const DistMatrix<F>& A, DistMatrix<F>& B, bool conjugated=false )
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
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const bool checkIfSingular = false;
    const auto d = A.GetDiagonal();
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    DiagonalSolve( LEFT, NORMAL, d, B, checkIfSingular );
    Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
}

template<typename F> 
inline void
SolveAfter
( const Matrix<F>& A, const Matrix<F>& dSub, const Matrix<Int>& p, 
  Matrix<F>& B, bool conjugated=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( p.Height() != A.Height() )
            LogicError("A and p must be the same height");
        // TODO: Check for dSub
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();
    ApplyRowPivots( B, p );
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    QuasiDiagonalSolve( LEFT, LOWER, NORMAL, d, dSub, B, conjugated );
    Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    ApplyInverseRowPivots( B, p );
}

template<typename F> 
inline void
SolveAfter
( const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dSub, 
  const DistMatrix<Int,VC,STAR>& p, DistMatrix<F>& B, bool conjugated=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::SolveAfter");
        if( A.Grid() != B.Grid() || A.Grid() != p.Grid() )
            LogicError("{A,B} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( A.Height() != p.Height() )
            LogicError("A and p must be the same height");
        // TODO: Check for dSub
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();
    ApplyRowPivots( B, p );
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    QuasiDiagonalSolve( LEFT, LOWER, NORMAL, d, dSub, B, conjugated );
    Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    ApplyInverseRowPivots( B, p );
}

} // namespace ldl
} // namespace elem

#endif // ifndef ELEM_LDL_SOLVEAFTER_HPP
