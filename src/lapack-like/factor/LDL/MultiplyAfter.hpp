/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_LDL_MULTIPLYAFTER_HPP
#define EL_LDL_MULTIPLYAFTER_HPP

namespace El {
namespace ldl {

template<typename F> 
void MultiplyAfter( const Matrix<F>& A, Matrix<F>& B, bool conjugated )
{
    DEBUG_ONLY(
        CallStackEntry cse("ldl::MultiplyAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    DiagonalScale( LEFT, NORMAL, d, B );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
}

template<typename F> 
void MultiplyAfter( const DistMatrix<F>& A, DistMatrix<F>& B, bool conjugated )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::MultiplyAfter");
        if( A.Grid() != B.Grid() )
            LogicError("{A,B} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    DiagonalScale( LEFT, NORMAL, d, B );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
}

template<typename F> 
void MultiplyAfter
( const Matrix<F>& A, const Matrix<F>& dSub, 
  const Matrix<Int>& pPerm, Matrix<F>& B, bool conjugated )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::MultiplyAfter");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( pPerm.Height() != A.Height() )
            LogicError("A and pPerm must be the same height");
        // TODO: Check for dSub
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();

    Matrix<Int> pInvPerm;
    InvertPermutation( pPerm, pInvPerm );
    
    PermuteRows( B, pPerm, pInvPerm );
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    QuasiDiagonalScale( LEFT, LOWER, d, dSub, B, conjugated );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    PermuteRows( B, pInvPerm, pPerm );
}

template<typename F,Dist UPerm> 
void MultiplyAfter
( const DistMatrix<F>& A, const DistMatrix<F,MD,STAR>& dSub, 
  const DistMatrix<Int,UPerm,STAR>& pPerm, DistMatrix<F>& B, bool conjugated )
{
    DEBUG_ONLY(
        CallStackEntry cse("lu::MultiplyAfter");
        if( A.Grid() != B.Grid() || A.Grid() != pPerm.Grid() )
            LogicError("{A,B,pPerm} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
        if( A.Height() != B.Height() )
            LogicError("A and B must be the same height");
        if( A.Height() != pPerm.Height() )
            LogicError("A and pPerm must be the same height");
        // TODO: Check for dSub
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = A.GetDiagonal();

    DistMatrix<Int,UPerm,STAR> pInvPerm(pPerm.Grid());
    InvertPermutation( pPerm, pInvPerm );

    PermuteRows( B, pPerm, pInvPerm );
    Trmm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    QuasiDiagonalScale( LEFT, LOWER, d, dSub, B, conjugated );
    Trmm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    PermuteRows( B, pInvPerm, pPerm );
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_MULTIPLYAFTER_HPP
