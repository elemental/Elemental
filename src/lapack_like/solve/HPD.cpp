/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace hpd_solve {

template<typename Field>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation,
  Matrix<Field>& A,
  Matrix<Field>& B )
{
    EL_DEBUG_CSE
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

template<typename Field>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation,
  AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& BPre )
{
    EL_DEBUG_CSE

    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    DistMatrixWriteProxy<Field,Field,MC,MR> BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

} // namespace hpd_solve

template<typename Field>
void HPDSolve
( UpperOrLower uplo,
  Orientation orientation,
  const Matrix<Field>& A,
        Matrix<Field>& B )
{
    EL_DEBUG_CSE
    Matrix<Field> ACopy( A );
    hpd_solve::Overwrite( uplo, orientation, ACopy, B );
}

template<typename Field>
void HPDSolve
( UpperOrLower uplo,
  Orientation orientation,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B )
{
    EL_DEBUG_CSE
    DistMatrix<Field> ACopy( A );
    hpd_solve::Overwrite( uplo, orientation, ACopy, B );
}

// TODO(poulson): Add iterative refinement parameter
template<typename Field>
void HPDSolve
( const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
    SparseLDLFactorization<Field> sparseLDLFact;
    const bool hermitian = true;
    sparseLDLFact.Initialize( A, hermitian, ctrl );
    sparseLDLFact.Factor();
    /*
    sparseLDLFact.SolveWithIterativeRefinement
    ( A, B, relTolRefine, maxRefineIts );
    */
    sparseLDLFact.Solve( B );
}

// TODO(poulson): Add iterative refinement parameter
template<typename Field>
void HPDSolve
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
    DistSparseLDLFactorization<Field> sparseLDLFact;
    const bool hermitian = true;
    sparseLDLFact.Initialize( A, hermitian, ctrl );
    sparseLDLFact.Factor();
    /*
    sparseLDLFact.SolveWithIterativeRefinement
    ( A, B, relTolRefine, maxRefineIts );
    */
    sparseLDLFact.Solve( B );
}

#define PROTO(Field) \
  template void hpd_solve::Overwrite \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<Field>& A, Matrix<Field>& B ); \
  template void hpd_solve::Overwrite \
  ( UpperOrLower uplo, Orientation orientation, \
    AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B ); \
  template void HPDSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<Field>& A, Matrix<Field>& B ); \
  template void HPDSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    const AbstractDistMatrix<Field>& A, AbstractDistMatrix<Field>& B ); \
  template void HPDSolve \
  ( const SparseMatrix<Field>& A, Matrix<Field>& B, const BisectCtrl& ctrl ); \
  template void HPDSolve \
  ( const DistSparseMatrix<Field>& A, DistMultiVec<Field>& B, \
    const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
