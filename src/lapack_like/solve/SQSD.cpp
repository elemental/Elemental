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

template<typename F>
void SQSDSolve( Int n0, UpperOrLower uplo, const Matrix<F>& A, Matrix<F>& B )
{
    DEBUG_CSE
    const Orientation orient = NORMAL;
    const bool conjugate = true;
    // TODO(poulson): LDLPivotCtrl control structure
    return SymmetricSolve( uplo, orient, A, B, conjugate );
}

template<typename F>
void SQSDSolve
( Int n0,
  UpperOrLower uplo,
  const AbstractDistMatrix<F>& A,
        AbstractDistMatrix<F>& B )
{
    DEBUG_CSE
    const Orientation orient = NORMAL;
    const bool conjugate = true;
    // TODO(poulson): LDLPivotCtrl control structure
    return SymmetricSolve( uplo, NORMAL, A, B, conjugate );
}

template<typename F>
void SQSDSolve
( Int n0,
        SparseMatrix<F>& A,
        Matrix<F>& B,
  const SQSDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
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
        A *= F(1)/normScale;
        B *= F(1)/normScale;
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

    if( ctrl.time )
        timer.Start();
    vector<Int> map, invMap;
    ldl::NodeInfo info;
    ldl::Separator rootSep;
    ldl::NestedDissection( AMod.LockedGraph(), map, rootSep, info );
    if( ctrl.time )
        Output("  ND: ",timer.Stop()," secs");
    InvertMap( map, invMap );
    ldl::Front<F> AModFront( AMod, map, info );
    if( ctrl.time )
        timer.Start();
    LDL( info, AModFront );
    if( ctrl.time )
        Output("  LDL: ",timer.Stop()," secs");

    // Successively solve each of the linear systems
    // =============================================
    if( ctrl.time )
        timer.Start();
    reg_ldl::SolveAfter
    ( A, regTmp, invMap, info, AModFront, B, ctrl.solveCtrl );
    if( ctrl.time )
        Output("  Solve: ",timer.Stop()," secs");
}

template<typename F>
void SQSDSolve
( Int n0,
  const SparseMatrix<F>& A,
        Matrix<F>& B,
  const SQSDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    auto ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.canOverwrite = true;
    SQSDSolve( n0, ACopy, B, ctrlMod );
}

template<typename F>
void SQSDSolve
( Int n0,
        DistSparseMatrix<F>& A,
        DistMultiVec<F>& B,
  const SQSDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    typedef Base<F> Real;
    if( !ctrl.canOverwrite )
    {
        auto ACopy( A ); 
        auto ctrlMod( ctrl );
        ctrlMod.canOverwrite = true;
        SQSDSolve( n0, ACopy, B, ctrlMod );
        return;
    }

    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const int commRank = mpi::Rank( comm );
    Timer timer;

    if( ctrl.scaleTwoNorm )
    {
        // Scale ABar down to roughly unit two-norm
        const Real normScale = TwoNormEstimate( A, ctrl.basisSize );
        if( ctrl.progress && commRank == 0 )
            Output("Estimated || A ||_2 ~= ",normScale);
        A *= F(1)/normScale;
        B *= F(1)/normScale;
    }
    // TODO(poulson): Add equilibration

    const Real scaledReg0Tmp =  ctrl.reg0Tmp*ctrl.reg0Tmp;
    const Real scaledReg1Tmp = -ctrl.reg1Tmp*ctrl.reg1Tmp;
    const Real scaledReg0Perm =  ctrl.reg0Perm*ctrl.reg0Perm;
    const Real scaledReg1Perm = -ctrl.reg1Perm*ctrl.reg1Perm;

    DistMultiVec<Real> regTmp(comm), regPerm(comm);
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

    if( commRank == 0 && ctrl.time )
        timer.Start();
    DistMap map, invMap;
    ldl::DistNodeInfo info;
    ldl::DistSeparator rootSep;
    ldl::NestedDissection( AMod.LockedDistGraph(), map, rootSep, info );
    if( commRank == 0 && ctrl.time )
        Output("  ND: ",timer.Stop()," secs");
    InvertMap( map, invMap );
    ldl::DistFront<F> AModFront( AMod, map, rootSep, info );

    if( commRank == 0 && ctrl.time )
        timer.Start();
    LDL( info, AModFront, LDL_2D );
    if( commRank == 0 && ctrl.time )
        Output("  LDL: ",timer.Stop()," secs");

    // Solve the k linear systems
    // ==========================
    if( commRank == 0 && ctrl.time )
        timer.Start();
    reg_ldl::SolveAfter
    ( A, regTmp, invMap, info, AModFront, B, ctrl.solveCtrl );
    if( commRank == 0 && ctrl.time )
        Output("  Solve: ",timer.Stop()," secs");
}

template<typename F>
void SQSDSolve
( Int n0,
  const DistSparseMatrix<F>& A,
        DistMultiVec<F>& B,
  const SQSDCtrl<Base<F>>& ctrl )
{
    DEBUG_CSE
    auto ACopy( A );
    auto ctrlMod( ctrl );
    ctrlMod.canOverwrite = true;
    SQSDSolve( n0, ACopy, B, ctrlMod );
}

#define PROTO(F) \
  template void SQSDSolve \
  ( Int n0, \
    UpperOrLower uplo, \
    const Matrix<F>& A, \
          Matrix<F>& B ); \
  template void SQSDSolve \
  ( Int n0, \
    UpperOrLower uplo, \
    const AbstractDistMatrix<F>& A, \
          AbstractDistMatrix<F>& B ); \
  template void SQSDSolve \
  ( Int n0, \
    const SparseMatrix<F>& A, \
          Matrix<F>& B, \
    const SQSDCtrl<Base<F>>& ctrl ); \
  template void SQSDSolve \
  ( Int n0, \
          SparseMatrix<F>& A, \
          Matrix<F>& B, \
    const SQSDCtrl<Base<F>>& ctrl ); \
  template void SQSDSolve \
  ( Int n0, \
    const DistSparseMatrix<F>& A, \
          DistMultiVec<F>& B, \
    const SQSDCtrl<Base<F>>& ctrl ); \
  template void SQSDSolve \
  ( Int n0, \
          DistSparseMatrix<F>& A, \
          DistMultiVec<F>& B, \
    const SQSDCtrl<Base<F>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace El
