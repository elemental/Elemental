/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_CHOLESKY_SOLVEAFTER_HPP
#define EL_CHOLESKY_SOLVEAFTER_HPP

namespace El {
namespace cholesky {

template<typename F> 
void SolveAfter
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
void SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const Matrix<F>& A, const Matrix<Int>& pPerm, Matrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::SolveAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( pPerm.Height() != A.Height() )
            LogicError("Permutation vector is wrong size");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    Matrix<Int> pInvPerm;
    InvertPermutation( pPerm, pInvPerm );

    PermuteRows( B, pPerm, pInvPerm );
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
    PermuteRows( B, pInvPerm, pPerm );
}

template<typename F> 
void SolveAfter
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

template<typename F,Dist UPerm> 
void SolveAfter
( UpperOrLower uplo, Orientation orientation, 
  const DistMatrix<F>& A, const DistMatrix<Int,UPerm,STAR>& pPerm, 
        DistMatrix<F>& B )
{
    DEBUG_ONLY(
        CallStackEntry cse("cholesky::SolveAfter");
        if( A.Grid() != B.Grid() )
            LogicError("{A,B} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != pPerm.Height() )
            LogicError("Permutation vector is wrong height");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    DistMatrix<Int,UPerm,STAR> pInvPerm(pPerm.Grid());
    InvertPermutation( pPerm, pInvPerm );

    PermuteRows( B, pPerm, pInvPerm );
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
    PermuteRows( B, pInvPerm, pPerm );
}

} // namespace cholesky
} // namespace El

#endif // ifndef EL_CHOLESKY_SOLVEAFTER_HPP
