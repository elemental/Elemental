/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {

// Solvers for Symmetric Quasi Semi-Definite matrices,
//
//   J = | F,    A |,
//       | A^T, -G |
//
// where F and G are Symmetric Positive Semi-Definite (and F is n0 x n0).

template<typename Field>
void SQSDSolve
( Int n0, UpperOrLower uplo, const Matrix<Field>& A, Matrix<Field>& B )
{
    EL_DEBUG_CSE
    const Orientation orient = NORMAL;
    const bool conjugate = true;
    // TODO(poulson): LDLPivotCtrl control structure
    return SymmetricSolve( uplo, orient, A, B, conjugate );
}

template<typename Field>
void SQSDSolve
( Int n0,
  UpperOrLower uplo,
  const AbstractDistMatrix<Field>& A,
        AbstractDistMatrix<Field>& B )
{
    EL_DEBUG_CSE
    const Orientation orient = NORMAL;
    const bool conjugate = true;
    // TODO(poulson): LDLPivotCtrl control structure
    return SymmetricSolve( uplo, orient, A, B, conjugate );
}

template<typename Field>
void SQSDSolve
( Int n0,
        SparseMatrix<Field>& A,
        Matrix<Field>& B,
  const SQSDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    if( !ctrl.canOverwrite )
    {
        auto ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.canOverwrite = true;
        SQSDSolve( n0, ACopy, B, ctrlMod );
        return;
    }

    if( ctrl.scaleTwoNorm )
    {
        // Scale ABar down to roughly unit two-norm
        const Real normScale = TwoNormEstimate( A, ctrl.basisSize );
        if( ctrl.progress )
            Output("Estimated || A ||_2 ~= ",normScale);
        A *= Field(1)/normScale;
        B *= Field(1)/normScale;
    }
    // TODO(poulson): Add equilibration

    const Int n = A.Height();
    const Real scaledReg0Tmp =  ctrl.reg0Tmp*ctrl.reg0Tmp;
    const Real scaledReg1Tmp = -ctrl.reg1Tmp*ctrl.reg1Tmp;
    const Real scaledReg0Perm =  ctrl.reg0Perm*ctrl.reg0Perm;
    const Real scaledReg1Perm = -ctrl.reg1Perm*ctrl.reg1Perm;
    Timer timer;

    Matrix<Real> regTmp, regPerm;
    regTmp.Resize( n, 1 );
    regPerm.Resize( n, 1 );
    for( Int i=0; i<n0; ++i )
    {
        regTmp(i) = scaledReg0Tmp;
        regPerm(i) = scaledReg0Perm;
    }
    for( Int i=n0; i<n; ++i )
    {
        regTmp(i) = scaledReg1Tmp;
        regPerm(i) = scaledReg1Perm;
    }
    UpdateRealPartOfDiagonal( A, Real(1), regPerm );
    auto AMod( A );
    UpdateRealPartOfDiagonal( AMod, Real(1), regTmp );

    // Analyze the sparsity pattern and initialize the frontal tree
    // ============================================================
    if( ctrl.time )
        timer.Start();
    const bool hermitian = true;
    const BisectCtrl bisectCtrl;
    SparseLDLFactorization<Field> sparseLDLFact;
    sparseLDLFact.Initialize( AMod, hermitian, bisectCtrl );
    if( ctrl.time )
        Output("  Analysis: ",timer.Stop()," secs");

    // Factor the sparse system
    // ========================
    if( ctrl.time )
        timer.Start();
    sparseLDLFact.Factor();
    if( ctrl.time )
        Output("  LDL: ",timer.Stop()," secs");

    // Successively solve each of the linear systems
    // =============================================
    if( ctrl.time )
        timer.Start();
    reg_ldl::SolveAfter( A, regTmp, sparseLDLFact, B, ctrl.solveCtrl );
    if( ctrl.time )
        Output("  Solve: ",timer.Stop()," secs");
}

template<typename Field>
void SQSDSolve
( Int n0,
  const SparseMatrix<Field>& A,
        Matrix<Field>& B,
  const SQSDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    auto ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.canOverwrite = true;
    SQSDSolve( n0, ACopy, B, ctrlMod );
}

template<typename Field>
void SQSDSolve
( Int n0,
        DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  const SQSDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    if( !ctrl.canOverwrite )
    {
        auto ACopy( A );
        auto ctrlMod( ctrl );
        ctrlMod.canOverwrite = true;
        SQSDSolve( n0, ACopy, B, ctrlMod );
        return;
    }

    const Int n = A.Height();
    const Grid& grid = A.Grid();
    const int commRank = grid.Rank();
    Timer timer;

    if( ctrl.scaleTwoNorm )
    {
        // Scale ABar down to roughly unit two-norm
        const Real normScale = TwoNormEstimate( A, ctrl.basisSize );
        if( ctrl.progress && commRank == 0 )
            Output("Estimated || A ||_2 ~= ",normScale);
        A *= Field(1)/normScale;
        B *= Field(1)/normScale;
    }
    // TODO(poulson): Add equilibration

    const Real scaledReg0Tmp =  ctrl.reg0Tmp*ctrl.reg0Tmp;
    const Real scaledReg1Tmp = -ctrl.reg1Tmp*ctrl.reg1Tmp;
    const Real scaledReg0Perm =  ctrl.reg0Perm*ctrl.reg0Perm;
    const Real scaledReg1Perm = -ctrl.reg1Perm*ctrl.reg1Perm;

    DistMultiVec<Real> regTmp(grid), regPerm(grid);
    regTmp.Resize( n, 1 );
    regPerm.Resize( n, 1 );
    for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
    {
        const Int i = regTmp.GlobalRow(iLoc);
        regTmp.Set( i, 0, i < n0 ? scaledReg0Tmp : scaledReg1Tmp );
        regPerm.Set( i, 0, i < n0 ? scaledReg0Perm : scaledReg1Perm );
    }
    UpdateRealPartOfDiagonal( A, Real(1), regPerm );
    auto AMod( A );
    UpdateRealPartOfDiagonal( AMod, Real(1), regTmp );

    // Analyze the sparsity pattern and initialize the frontal tree
    // ============================================================
    if( commRank == 0 && ctrl.time )
        timer.Start();
    const bool hermitian = true;
    const BisectCtrl bisectCtrl;
    DistSparseLDLFactorization<Field> sparseLDLFact;
    sparseLDLFact.Initialize( AMod, hermitian, bisectCtrl );
    if( commRank == 0 && ctrl.time )
        Output("  Analysis: ",timer.Stop()," secs");

    // Factor the sparse system
    // ========================
    if( commRank == 0 && ctrl.time )
        timer.Start();
    sparseLDLFact.Factor( LDL_2D );
    if( commRank == 0 && ctrl.time )
        Output("  LDL: ",timer.Stop()," secs");

    // Solve the k linear systems
    // ==========================
    if( commRank == 0 && ctrl.time )
        timer.Start();
    reg_ldl::SolveAfter( A, regTmp, sparseLDLFact, B, ctrl.solveCtrl );
    if( commRank == 0 && ctrl.time )
        Output("  Solve: ",timer.Stop()," secs");
}

template<typename Field>
void SQSDSolve
( Int n0,
  const DistSparseMatrix<Field>& A,
        DistMultiVec<Field>& B,
  const SQSDCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    auto ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.canOverwrite = true;
    SQSDSolve( n0, ACopy, B, ctrlMod );
}

#define PROTO(Field) \
  template void SQSDSolve \
  ( Int n0, \
    UpperOrLower uplo, \
    const Matrix<Field>& A, \
          Matrix<Field>& B ); \
  template void SQSDSolve \
  ( Int n0, \
    UpperOrLower uplo, \
    const AbstractDistMatrix<Field>& A, \
          AbstractDistMatrix<Field>& B ); \
  template void SQSDSolve \
  ( Int n0, \
    const SparseMatrix<Field>& A, \
          Matrix<Field>& B, \
    const SQSDCtrl<Base<Field>>& ctrl ); \
  template void SQSDSolve \
  ( Int n0, \
          SparseMatrix<Field>& A, \
          Matrix<Field>& B, \
    const SQSDCtrl<Base<Field>>& ctrl ); \
  template void SQSDSolve \
  ( Int n0, \
    const DistSparseMatrix<Field>& A, \
          DistMultiVec<Field>& B, \
    const SQSDCtrl<Base<Field>>& ctrl ); \
  template void SQSDSolve \
  ( Int n0, \
          DistSparseMatrix<Field>& A, \
          DistMultiVec<Field>& B, \
    const SQSDCtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
