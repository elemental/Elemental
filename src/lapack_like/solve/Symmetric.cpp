/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

namespace symm_solve {

template<typename Field>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation,
  Matrix<Field>& A,
  Matrix<Field>& B,
  bool conjugate,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    if( uplo == UPPER )
        LogicError("Upper Bunch-Kaufman is not yet supported");
    Permutation p;
    Matrix<Field> dSub;
    LDL( A, dSub, p, conjugate, ctrl );
    const bool conjFlip = ( (orientation == ADJOINT && conjugate == false) ||
                            (orientation == TRANSPOSE && conjugate == true) );
    if( conjFlip )
        Conjugate( B );
    ldl::SolveAfter( A, dSub, p, B, conjugate );
    if( conjFlip )
        Conjugate( B );
}

template<typename Field>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation,
  AbstractDistMatrix<Field>& APre,
  AbstractDistMatrix<Field>& BPre,
  bool conjugate,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    if( uplo == UPPER )
        LogicError("Upper Bunch-Kaufman is not yet supported");

    DistMatrixReadProxy<Field,Field,MC,MR> AProx( APre );
    DistMatrixReadWriteProxy<Field,Field,MC,MR> BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    DistPermutation p(A.Grid());
    DistMatrix<Field,MD,STAR> dSub(A.Grid());
    LDL( A, dSub, p, conjugate, ctrl );
    const bool conjFlip = ( (orientation == ADJOINT && conjugate == false) ||
                            (orientation == TRANSPOSE && conjugate == true) );
    if( conjFlip )
        Conjugate( B );
    ldl::SolveAfter( A, dSub, p, B, conjugate );
    if( conjFlip )
        Conjugate( B );
}

} // namespace symm_solve

template<typename Field>
void SymmetricSolve
( UpperOrLower uplo,
  Orientation orientation,
  const Matrix<Field>& A,
        Matrix<Field>& B,
  bool conjugate,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    Matrix<Field> ACopy( A );
    symm_solve::Overwrite( uplo, orientation, ACopy, B, conjugate, ctrl );
}

template<typename Field>
void SymmetricSolve
( UpperOrLower uplo,
  Orientation orientation,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B,
  bool conjugate,
  const LDLPivotCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    DistMatrix<Field> ACopy( A );
    symm_solve::Overwrite( uplo, orientation, ACopy, B, conjugate, ctrl );
}

// TODO(poulson): Add iterative refinement parameter
template<typename Field>
void SymmetricSolve
( const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  bool conjugate,
  bool tryLDL,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE

    if( tryLDL )
    {
        ldl::NodeInfo info;
        ldl::Separator rootSep;
        vector<Int> map, invMap;
        ldl::NestedDissection( A.LockedGraph(), map, rootSep, info, ctrl );
        InvertMap( map, invMap );

        ldl::Front<Field> front( A, map, info, conjugate );
        LDL( info, front );

        // TODO(poulson): Extend ldl::SolveWithIterativeRefinement to support
        // multiple right-hand sides
        /*
        ldl::SolveWithIterativeRefinement
        ( A, invMap, info, front, B, minReductionFactor, maxRefineIts );
        */
        ldl::SolveAfter( invMap, info, front, B );
    }
    else
    {
        LinearSolve( A, B );
    }
}

// TODO(poulson): Add iterative refinement parameter
template<typename Field>
void SymmetricSolve
( const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  bool conjugate,
  bool tryLDL,
  const BisectCtrl& ctrl )
{
    EL_DEBUG_CSE
    if( tryLDL )
    {
        ldl::DistNodeInfo info;
        ldl::DistSeparator rootSep;
        DistMap map, invMap;
        ldl::NestedDissection( A.LockedDistGraph(), map, rootSep, info, ctrl );
        InvertMap( map, invMap );

        ldl::DistFront<Field> front( A, map, rootSep, info, conjugate );
        LDL( info, front, LDL_INTRAPIV_1D );

        // TODO(poulson): Extend ldl::SolveWithIterativeRefinement to support
        // multiple right-hand sides
        /*
        ldl::SolveWithIterativeRefinement
        ( A, invMap, info, front, B, minReductionFactor, maxRefineIts );
        */
        ldl::SolveAfter( invMap, info, front, B );
    }
    else
    {
        LinearSolve( A, B );
    }
}

#define PROTO(Field) \
  template void symm_solve::Overwrite \
  ( UpperOrLower uplo, \
    Orientation orientation, \
    Matrix<Field>& A, \
    Matrix<Field>& B, \
    bool conjugate, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void symm_solve::Overwrite \
  ( UpperOrLower uplo, \
    Orientation orientation, \
    AbstractDistMatrix<Field>& A, \
    AbstractDistMatrix<Field>& B, \
    bool conjugate, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void SymmetricSolve \
  ( UpperOrLower uplo, \
    Orientation orientation, \
    const Matrix<Field>& A, \
          Matrix<Field>& B, \
    bool conjugate, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void SymmetricSolve \
  ( UpperOrLower uplo, \
    Orientation orientation, \
    const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& B, \
    bool conjugate, \
    const LDLPivotCtrl<Base<Field>>& ctrl ); \
  template void SymmetricSolve \
  ( const SparseMatrix<Field>& A, \
          Matrix<Field>& B, \
    bool conjugate, \
    bool tryLDL, \
    const BisectCtrl& ctrl ); \
  template void SymmetricSolve \
  ( const DistSparseMatrix<Field>& A, \
          DistMultiVec<Field>& B, \
    bool conjugate, \
    bool tryLDL, \
    const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
