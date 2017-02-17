/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace herm_solve {

template<typename Field>
void Overwrite
( UpperOrLower uplo, Orientation orientation,
  Matrix<Field>& A, Matrix<Field>& B,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    symm_solve::Overwrite( uplo, orientation, A, B, true, ctrl );
}

template<typename Field>
void Overwrite
( UpperOrLower uplo, Orientation orientation,
  AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    symm_solve::Overwrite( uplo, orientation, A, B, true, ctrl );
}

} // namespace herm_solve

template<typename Field>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation,
  const Matrix<Field>& A, Matrix<Field>& B,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SymmetricSolve( uplo, orientation, A, B, true, ctrl );
}

template<typename Field>
void HermitianSolve
( UpperOrLower uplo, Orientation orientation,
  const AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    SymmetricSolve( uplo, orientation, A, B, true, ctrl );
}

// TODO(poulson): Add iterative refinement parameter
template<typename Field>
void HermitianSolve
( const SparseMatrix<Field>& A, Matrix<Field>& B,
  bool tryLDL, const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
    SymmetricSolve( A, B, true, tryLDL, ctrl );
}

// TODO(poulson): Add iterative refinement parameter
template<typename Field>
void HermitianSolve
( const DistSparseMatrix<Field>& A, DistMultiVec<Field>& B,
  bool tryLDL, const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
    SymmetricSolve( A, B, true, tryLDL, ctrl );
}

#define PROTO(Field) \
  template void herm_solve::Overwrite \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<Field>& A, Matrix<Field>& B, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void herm_solve::Overwrite \
  ( UpperOrLower uplo, Orientation orientation, \
    AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void HermitianSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<Field>& A, Matrix<Field>& B, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void HermitianSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void HermitianSolve \
  ( const SparseMatrix<Field>& A, Matrix<Field>& B, \
    bool tryLDL, const BisectCtrl& ctrl ); \
  template void HermitianSolve \
  ( const DistSparseMatrix<Field>& A, DistMultiVec<Field>& B, \
    bool tryLDL, const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
