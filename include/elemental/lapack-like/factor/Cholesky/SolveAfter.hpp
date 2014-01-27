/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_CHOLESKY_SOLVEAFTER_HPP
#define ELEM_CHOLESKY_SOLVEAFTER_HPP

#include ELEM_CONJUGATE_INC
#include ELEM_APPLYROWPIVOTS_INC
#include ELEM_TRSM_INC

namespace elem {
namespace cholesky {

template<typename F> 
inline void
SolveAfter
( UpperOrLower uplo, Orientation orientation, const Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::SolveAfter");
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
inline void
SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, const Matrix<Int>& p, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::SolveAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( p.Height() != A.Height() )
            LogicError("Pivot vector is wrong size");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    ApplyRowPivots( B, p );
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
    ApplyInverseRowPivots( B, p );
}

template<typename F> 
inline void
SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::SolveAfter");
        if( A.Grid() != B.Grid() )
            LogicError("{A,B} must be distributed over the same grid");
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
inline void
SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<Int,VC,STAR>& p, DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::SolveAfter");
        if( A.Grid() != B.Grid() )
            LogicError("{A,B} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != p.Height() )
            LogicError("Pivot vector is wrong height");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    ApplyRowPivots( B, p );
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
    ApplyInverseRowPivots( B, p );
}

} // namespace cholesky
} // namespace elem

#endif // ifndef ELEM_CHOLESKY_SOLVEAFTER_HPP
