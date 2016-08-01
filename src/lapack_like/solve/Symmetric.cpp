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

template<typename F>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation, 
  Matrix<F>& A,
  Matrix<F>& B, 
  bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    if( uplo == UPPER )
        LogicError("Upper Bunch-Kaufman is not yet supported");
    Permutation p; 
    Matrix<F> dSub;
    LDL( A, dSub, p, conjugate, ctrl );
    const bool conjFlip = ( (orientation == ADJOINT && conjugate == false) ||
                            (orientation == TRANSPOSE && conjugate == true) );
    if( conjFlip )
        Conjugate( B );
    ldl::SolveAfter( A, dSub, p, B, conjugate );
    if( conjFlip )
        Conjugate( B );
}

template<typename F>
void Overwrite
( UpperOrLower uplo,
  Orientation orientation, 
  ElementalMatrix<F>& APre,
  ElementalMatrix<F>& BPre, 
  bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    if( uplo == UPPER )
        LogicError("Upper Bunch-Kaufman is not yet supported");

    DistMatrixReadProxy<F,F,MC,MR> AProx( APre );
    DistMatrixReadWriteProxy<F,F,MC,MR> BProx( BPre );
    auto& A = AProx.Get();
    auto& B = BProx.Get();

    DistPermutation p(A.Grid()); 
    DistMatrix<F,MD,STAR> dSub(A.Grid());
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

template<typename F>
void SymmetricSolve
( UpperOrLower uplo,
  Orientation orientation, 
  const Matrix<F>& A,
        Matrix<F>& B, 
  bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    Matrix<F> ACopy( A );
    symm_solve::Overwrite( uplo, orientation, ACopy, B, conjugate, ctrl );
}

template<typename F>
void SymmetricSolve
( UpperOrLower uplo,
  Orientation orientation, 
  const ElementalMatrix<F>& A,
        ElementalMatrix<F>& B, 
  bool conjugate,
  const LDLPivotCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    DistMatrix<F> ACopy( A );
    symm_solve::Overwrite( uplo, orientation, ACopy, B, conjugate, ctrl );
}

// TODO: Add iterative refinement parameter
template<typename F>
void SymmetricSolve
( const SparseMatrix<F>& A,
        Matrix<F>& B,
  bool conjugate,
  bool tryLDL,
  const BisectCtrl& ctrl )
{
    DEBUG_CSE

    if( tryLDL )
    {
        ldl::NodeInfo info;
        ldl::Separator rootSep;
        vector<Int> map, invMap;
        ldl::NestedDissection( A.LockedGraph(), map, rootSep, info, ctrl );
        InvertMap( map, invMap );

        ldl::Front<F> front( A, map, info, conjugate );
        LDL( info, front );

        // TODO: Extend ldl::SolveWithIterativeRefinement to support multiple
        //       right-hand sides
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

// TODO: Add iterative refinement parameter
template<typename F>
void SymmetricSolve
( const DistSparseMatrix<F>& A,
        DistMultiVec<F>& B,
  bool conjugate,
  bool tryLDL,
  const BisectCtrl& ctrl )
{
    DEBUG_CSE
    if( tryLDL )
    {
        ldl::DistNodeInfo info;
        ldl::DistSeparator rootSep;
        DistMap map, invMap;
        ldl::NestedDissection( A.LockedDistGraph(), map, rootSep, info, ctrl );
        InvertMap( map, invMap );

        ldl::DistFront<F> front( A, map, rootSep, info, conjugate );
        LDL( info, front, LDL_INTRAPIV_1D );

        // TODO: Extend ldl::SolveWithIterativeRefinement to support multiple
        //       right-hand sides
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

#define PROTO(F) \
  template void symm_solve::Overwrite \
  ( UpperOrLower uplo, Orientation orientation, \
    Matrix<F>& A, Matrix<F>& B, bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void symm_solve::Overwrite \
  ( UpperOrLower uplo, Orientation orientation, \
    ElementalMatrix<F>& A, ElementalMatrix<F>& B, bool conjugate, \
    const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void SymmetricSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    const Matrix<F>& A, Matrix<F>& B, \
    bool conjugate, const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void SymmetricSolve \
  ( UpperOrLower uplo, Orientation orientation, \
    const ElementalMatrix<F>& A, ElementalMatrix<F>& B, \
    bool conjugate, const LDLPivotCtrl<Base<F>>& ctrl ); \
  template void SymmetricSolve \
  ( const SparseMatrix<F>& A, Matrix<F>& B, \
    bool conjugate, bool tryLDL, const BisectCtrl& ctrl ); \
  template void SymmetricSolve \
  ( const DistSparseMatrix<F>& A, DistMultiVec<F>& B, \
    bool conjugate, bool tryLDL, const BisectCtrl& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
