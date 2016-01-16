/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {

namespace hpd_solve {

template<typename F>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation, 
  Matrix<F>& A,
  Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("hpd_solve::Overwrite"))
    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

template<typename F>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation, 
  ElementalMatrix<F>& APre,
  ElementalMatrix<F>& BPre )
{
    DEBUG_ONLY(CSE cse("hpd_solve::Overwrite"))

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,MC,MR> BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    Cholesky( uplo, A );
    cholesky::SolveAfter( uplo, orientation, A, B );
}

} // namespace hpd_solve

template<typename F>
void HPDSolve
( UpperOrLower uplo,
  Orientation orientation, 
  const Matrix<F>& A,
        Matrix<F>& B )
{
    DEBUG_ONLY(CSE cse("HPDSolve"))
    Matrix<F> ACopy( A );
    hpd_solve::Overwrite( uplo, orientation, ACopy, B );
}

template<typename F>
void HPDSolve
( UpperOrLower uplo,
  Orientation orientation, 
  const ElementalMatrix<F>& A,
        ElementalMatrix<F>& B )
{
    DEBUG_ONLY(CSE cse("HPDSolve"))
    DistMatrix<F> ACopy( A );
    hpd_solve::Overwrite( uplo, orientation, ACopy, B );
}

// TODO: Add iterative refinement parameter
template<typename F>
void HPDSolve
( const SparseMatrix<F>& A,
        Matrix<F>& B,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("HPDSolve"))
    ldl::NodeInfo info;
    ldl::Separator rootSep;
    vector<Int> map, invMap;
    ldl::NestedDissection( A.LockedGraph(), map, rootSep, info, ctrl );
    InvertMap( map, invMap );

    ldl::Front<F> front( A, map, info, true );
    LDL( info, front );

    // TODO: Extend ldl::SolveWithIterativeRefinement to support multiple
    //       right-hand sides
    /*
    ldl::SolveWithIterativeRefinement
    ( A, invMap, info, front, B, minReductionFactor, maxRefineIts );
    */
    ldl::SolveAfter( invMap, info, front, B );
}

// TODO: Add iterative refinement parameter
template<typename F>
void HPDSolve
( const DistSparseMatrix<F>& A,
        DistMultiVec<F>& B,
  const BisectCtrl& ctrl )
{
    DEBUG_ONLY(CSE cse("HPDSolve"))
    ldl::DistNodeInfo info;
    ldl::DistSeparator rootSep;
    DistMap map, invMap;
    ldl::NestedDissection( A.LockedDistGraph(), map, rootSep, info, ctrl );
    InvertMap( map, invMap );

    ldl::DistFront<F> front( A, map, rootSep, info, true );
    LDL( info, front );

    // TODO: Extend ldl::SolveWithIterativeRefinement to support multiple
    //       right-hand sides
    /*
    ldl::SolveWithIterativeRefinement
    ( A, invMap, info, front, B, minReductionFactor, maxRefineIts );
    */
    ldl::SolveAfter( invMap, info, front, B );
}

#define PROTO(F) \
  template void hpd_solve::Overwrite \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<F>& A, Matrix<F>& B ); \
  template void hpd_solve::Overwrite \
  ( UpperOrLower uplo, Orientation orientation, \
    ElementalMatrix<F>& A, ElementalMatrix<F>& B ); \
  template void HPDSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, Matrix<F>& B ); \
  template void HPDSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    const ElementalMatrix<F>& A, ElementalMatrix<F>& B ); \
  template void HPDSolve \
  ( const SparseMatrix<F>& A, Matrix<F>& B, const BisectCtrl& ctrl ); \
  template void HPDSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& B, const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include "El/macros/Instantiate.h"

} // namespace El
