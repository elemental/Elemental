/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LDL_SOLVEAFTER_HPP
#define EL_LDL_SOLVEAFTER_HPP

namespace El {
namespace ldl {

template<typename F> 
void SolveAfter( const Matrix<F>& A, Matrix<F>& B, bool conjugated )
{
    DEBUG_ONLY(
      CSE cse("ldl::SolveAfter");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const bool checkIfSingular = false;
    const auto d = GetDiagonal(A);
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    DiagonalSolve( LEFT, NORMAL, d, B, checkIfSingular );
    Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
}

template<typename F> 
void SolveAfter
( const ElementalMatrix<F>& APre,
        ElementalMatrix<F>& B,
  bool conjugated )
{
    DEBUG_ONLY(
      CSE cse("ldl::SolveAfter");
      AssertSameGrids( APre, B );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
      if( APre.Height() != B.Height() )
          LogicError("A and B must be the same height");
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const bool checkIfSingular = false;

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.GetLocked();

    const auto d = GetDiagonal(A);

    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    DiagonalSolve( LEFT, NORMAL, d, B, checkIfSingular );
    Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
}

template<typename F> 
void SolveAfter
( const Matrix<F>& A,
  const Matrix<F>& dSub, 
  const Permutation& P,
        Matrix<F>& B,
  bool conjugated )
{
    DEBUG_ONLY(
      CSE cse("ldl::SolveAfter");
      if( A.Height() != A.Width() )
          LogicError("A must be square");
      if( A.Height() != B.Height() )
          LogicError("A and B must be the same height");
      // TODO: Check for dSub
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );
    const auto d = GetDiagonal(A);

    P.PermuteRows( B );
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    QuasiDiagonalSolve( LEFT, LOWER, d, dSub, B, conjugated );
    Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    P.InversePermuteRows( B );
}

template<typename F> 
void SolveAfter
( const ElementalMatrix<F>& APre,
  const ElementalMatrix<F>& dSub, 
  const DistPermutation& P,
        ElementalMatrix<F>& BPre, 
  bool conjugated )
{
    DEBUG_ONLY(
      CSE cse("ldl::SolveAfter");
      AssertSameGrids( APre, BPre );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
      if( APre.Height() != BPre.Height() )
          LogicError("A and B must be the same height");
      // TODO: Check for dSub
    )
    const Orientation orientation = ( conjugated ? ADJOINT : TRANSPOSE );

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixReadWriteProxy<F,F,MC,MR> BProx( BPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.Get();

    const auto d = GetDiagonal(A);

    P.PermuteRows( B );
    Trsm( LEFT, LOWER, NORMAL, UNIT, F(1), A, B );
    QuasiDiagonalSolve( LEFT, LOWER, d, dSub, B, conjugated );
    Trsm( LEFT, LOWER, orientation, UNIT, F(1), A, B );
    P.InversePermuteRows( B );
}

} // namespace ldl
} // namespace El

#endif // ifndef EL_LDL_SOLVEAFTER_HPP
