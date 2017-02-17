/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace El {
namespace reg_ldl {

template<typename Field>
Int RegularizedSolveAfterNoPromote
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    EL_DEBUG_CSE
    auto applyA =
      [&]( const Matrix<Field>& X, Matrix<Field>& Y )
      {
        Y = X;
        DiagonalScale( LEFT, NORMAL, reg, Y ); 
        Multiply( NORMAL, Field(1), A, X, Field(1), Y );
      };
    auto applyAInv =
      [&]( Matrix<Field>& Y )
      {
        sparseLDLFact.Solve( Y );
      };
    return
      RefinedSolve 
      ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename Field>
Int RegularizedSolveAfterNoPromote
( const SparseMatrix<Field>& A,
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d, 
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    EL_DEBUG_CSE

    // TODO(poulson): Use time in these lambdas
    auto applyA =
      [&]( const Matrix<Field>& X, Matrix<Field>& Y )
      {
        Y = X;
        DiagonalScale( LEFT, NORMAL, reg, Y ); 
        Multiply( NORMAL, Field(1), A, X, Field(1), Y );
      };
    auto applyAInv = 
      [&]( Matrix<Field>& Y )
      {
        DiagonalSolve( LEFT, NORMAL, d, Y );
        sparseLDLFact.Solve( Y );
        DiagonalSolve( LEFT, NORMAL, d, Y );
      };
    return RefinedSolve( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename Field>
DisableIf<IsSame<Field,Promote<Field>>,Int>
RegularizedSolveAfterPromote
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    
    typedef Base<Field> Real; 
    typedef Promote<Real> PReal;
    typedef Promote<Field> PField;

    // TODO(poulson): Avoid reforming these each call
    SparseMatrix<PField> AProm;
    Copy( A, AProm );
    Matrix<PReal> regProm;
    Copy( reg, regProm );

    // TODO: Use time in these lambdas
    auto applyA =
      [&]( const Matrix<PField>& XProm, Matrix<PField>& YProm )
      {
        YProm = XProm; 
        DiagonalScale( LEFT, NORMAL, regProm, YProm ); 
        Multiply( NORMAL, PField(1), AProm, XProm, PField(1), YProm );
      }; 
    auto applyAInv =  
      [&]( Matrix<Field>& Y )
      {
        sparseLDLFact.Solve( Y );
      };

    return PromotedRefinedSolve
      ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename Field>
EnableIf<IsSame<Field,Promote<Field>>,Int>
RegularizedSolveAfterPromote
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    return RegularizedSolveAfterNoPromote
      ( A, reg, sparseLDLFact, B, relTol, maxRefineIts, progress, time );
}

template<typename Field>
DisableIf<IsSame<Field,Promote<Field>>,Int>
RegularizedSolveAfterPromote
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d, 
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Promote<Real> PReal;
    typedef Promote<Field> PField;

    // TODO(poulson): Avoid reforming these each call
    SparseMatrix<PField> AProm;
    Copy( A, AProm );
    Matrix<PReal> regProm;
    Copy( reg, regProm );

    // TODO(poulson): Use time in these lambdas
    auto applyA =
      [&]( const Matrix<PField>& XProm, Matrix<PField>& YProm )
      {
        YProm = XProm; 
        DiagonalScale( LEFT, NORMAL, regProm, YProm ); 
        Multiply( NORMAL, PField(1), AProm, XProm, PField(1), YProm );
      }; 
    auto applyAInv =  
      [&]( Matrix<Field>& Y )
      {
        DiagonalSolve( LEFT, NORMAL, d, Y );
        sparseLDLFact.Solve( Y );
        DiagonalSolve( LEFT, NORMAL, d, Y );
      };

    return PromotedRefinedSolve
      ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename Field>
EnableIf<IsSame<Field,Promote<Field>>,Int>
RegularizedSolveAfterPromote
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d, 
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    return RegularizedSolveAfter
      ( A, reg, d, sparseLDLFact, B, relTol, maxRefineIts, progress, time );
}

template<typename Field>
Int RegularizedSolveAfter
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts, 
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    return RegularizedSolveAfterPromote
           ( A, reg, sparseLDLFact, B, relTol, maxRefineIts, progress, time );
}

template<typename Field>
Int RegularizedSolveAfter
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d, 
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    return RegularizedSolveAfterPromote
           ( A, reg, d, sparseLDLFact, B, relTol, maxRefineIts,
             progress, time );
}

template<typename Field>
Int RegularizedSolveAfterNoPromote
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    auto applyA =
      [&]( const DistMultiVec<Field>& X, DistMultiVec<Field>& Y )
      {
        Y = X;
        DiagonalScale( LEFT, NORMAL, reg, Y ); 
        Multiply( NORMAL, Field(1), A, X, Field(1), Y );
      };
    auto applyAInv =
      [&]( DistMultiVec<Field>& Y )
      {
        sparseLDLFact.Solve( Y );
      };
    return RefinedSolve
      ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename Field>
Int RegularizedSolveAfterNoPromote
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    auto applyA =
      [&]( const DistMultiVec<Field>& X, DistMultiVec<Field>& Y )
      {
        Y = X;
        DiagonalScale( LEFT, NORMAL, reg, Y ); 
        Multiply( NORMAL, Field(1), A, X, Field(1), Y );
      };
    auto applyAInv = 
      [&]( DistMultiVec<Field>& Y )
      {
        DiagonalSolve( LEFT, NORMAL, d, Y );
        sparseLDLFact.Solve( Y );
        DiagonalSolve( LEFT, NORMAL, d, Y );
      };
    return RefinedSolve( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename Field>
DisableIf<IsSame<Field,Promote<Field>>,Int>
RegularizedSolveAfterPromote
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Promote<Real> PReal;
    typedef Promote<Field> PField;

    // TODO(poulson): Perform these conversions less frequently at a higher
    // level
    DistSparseMatrix<PField> AProm(A.Grid());
    Copy( A, AProm );
    DistMultiVec<PReal> regProm(reg.Grid());
    Copy( reg, regProm );

    // TODO(poulson): Use time in these lambdas
    auto applyA =
      [&]( const DistMultiVec<PField>& XProm, DistMultiVec<PField>& YProm )
      {
        YProm = XProm;
        DiagonalScale( LEFT, NORMAL, regProm, YProm ); 
        Multiply( NORMAL, PField(1), AProm, XProm, PField(1), YProm );
      };
    auto applyAInv =
      [&]( DistMultiVec<Field>& Y )
      {
        sparseLDLFact.Solve( Y );
      };
    return PromotedRefinedSolve
      ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename Field>
EnableIf<IsSame<Field,Promote<Field>>,Int>
RegularizedSolveAfterPromote
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    return RegularizedSolveAfterNoPromote
      ( A, reg, sparseLDLFact, B, relTol, maxRefineIts, progress, time );
}

template<typename Field>
DisableIf<IsSame<Field,Promote<Field>>,Int>
RegularizedSolveAfterPromote
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    typedef Base<Field> Real;
    typedef Promote<Real> PReal;
    typedef Promote<Field> PField;

    // TODO(poulson): Perform these conversions less frequently at a higher
    // level
    DistSparseMatrix<PField> AProm(A.Grid());
    Copy( A, AProm );
    DistMultiVec<PReal> regProm(reg.Grid());
    Copy( reg, regProm );

    // TODO(poulson): Use time in these lambdas
    auto applyA =
      [&]( const DistMultiVec<PField>& XProm, DistMultiVec<PField>& YProm )
      {
        YProm = XProm;
        DiagonalScale( LEFT, NORMAL, regProm, YProm ); 
        Multiply( NORMAL, PField(1), AProm, XProm, PField(1), YProm );
      };
    auto applyAInv = 
      [&]( DistMultiVec<Field>& Y )
      {
        DiagonalSolve( LEFT, NORMAL, d, Y );
        sparseLDLFact.Solve( Y );
        DiagonalSolve( LEFT, NORMAL, d, Y );
      };

    return PromotedRefinedSolve
           ( applyA, applyAInv, B, relTol, maxRefineIts, progress );
}

template<typename Field>
EnableIf<IsSame<Field,Promote<Field>>,Int>
RegularizedSolveAfterPromote
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    return RegularizedSolveAfterNoPromote
      ( A, reg, d, sparseLDLFact, B, relTol, maxRefineIts, progress, time );
}

template<typename Field>
Int RegularizedSolveAfter
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    return RegularizedSolveAfterPromote
      ( A, reg, sparseLDLFact, B, relTol, maxRefineIts, progress, time );
}

template<typename Field>
Int RegularizedSolveAfter
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d, 
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int maxRefineIts,
  bool progress,
  bool time )
{
    EL_DEBUG_CSE
    return RegularizedSolveAfterPromote
      ( A, reg, d, sparseLDLFact, B, relTol, maxRefineIts, progress, time );
}

template<typename Field>
Int LGMRESSolveAfter
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int restart,
  Int maxIts,
  Base<Field> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    EL_DEBUG_CSE

    auto applyA =
      [&]( Field alpha, const Matrix<Field>& X, Field beta, Matrix<Field>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( Matrix<Field>& W )
      {
        RegularizedSolveAfter
        ( A, reg, sparseLDLFact, W, relTolRefine, maxRefineIts, progress );
      };

    return LGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename Field>
Int LGMRESSolveAfter
( const SparseMatrix<Field>& A,
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  Base<Field> relTol,
  Int restart,
  Int maxIts,
  Base<Field> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    EL_DEBUG_CSE

    auto applyA =
      [&]( Field alpha, const Matrix<Field>& X, Field beta, Matrix<Field>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( Matrix<Field>& W )
      {
        RegularizedSolveAfter
        ( A, reg, d, sparseLDLFact, W, relTolRefine, maxRefineIts, progress );
      };

    return LGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename Field>
Int LGMRESSolveAfter
( const DistSparseMatrix<Field>& A,
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  Base<Field> relTol,
  Int restart,
  Int maxIts,
  Base<Field> relTolRefine,
  Int maxRefineIts,
  bool progress )
{
    EL_DEBUG_CSE

    auto applyA =
      [&]( Field alpha, const DistMultiVec<Field>& X,
           Field beta, DistMultiVec<Field>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( DistMultiVec<Field>& W )
      {
        RegularizedSolveAfter
        ( A, reg, sparseLDLFact, W, relTolRefine, maxRefineIts, progress );
      };

    return LGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename Field>
Int LGMRESSolveAfter
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        Base<Field> relTolRefine,
        Int maxRefineIts,
        bool progress )
{
    EL_DEBUG_CSE
    auto applyA =
      [&]( Field alpha, const DistMultiVec<Field>& X,
           Field beta, DistMultiVec<Field>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( DistMultiVec<Field>& W )
      {
        RegularizedSolveAfter
        ( A, reg, d, sparseLDLFact, W, relTolRefine, maxRefineIts, progress );
      };

    return LGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename Field>
Int FGMRESSolveAfter
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        Base<Field> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    EL_DEBUG_CSE
    auto applyA =
      [&]( Field alpha, const Matrix<Field>& X, Field beta, Matrix<Field>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( Matrix<Field>& W )
      {
        RegularizedSolveAfter
        ( A, reg, sparseLDLFact, W, relTolRefine, maxRefineIts, progress );
      };

    return FGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename Field>
Int FGMRESSolveAfter
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        Base<Field> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    EL_DEBUG_CSE

    auto applyA =
      [&]( Field alpha, const Matrix<Field>& X, Field beta, Matrix<Field>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( Matrix<Field>& W )
      {
        RegularizedSolveAfter
        ( A, reg, d, sparseLDLFact, W, relTolRefine, maxRefineIts, progress );
      };

    return FGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename Field>
Int FGMRESSolveAfter
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        Base<Field> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    EL_DEBUG_CSE

    auto applyA =
      [&]( Field alpha, const DistMultiVec<Field>& X,
           Field beta, DistMultiVec<Field>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( DistMultiVec<Field>& W )
      {
        RegularizedSolveAfter
        ( A, reg, sparseLDLFact, W, relTolRefine, maxRefineIts, progress );
      };

    return FGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

template<typename Field>
Int FGMRESSolveAfter
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        Base<Field> relTolRefine,
        Int maxRefineIts, 
        bool progress,
        bool time )
{
    EL_DEBUG_CSE

    auto applyA =
      [&]( Field alpha, const DistMultiVec<Field>& X,
           Field beta, DistMultiVec<Field>& Y )
      {
          Multiply( NORMAL, alpha, A, X, beta, Y );
      };
    auto precond =
      [&]( DistMultiVec<Field>& W )
      {
        RegularizedSolveAfter
        ( A, reg, d, sparseLDLFact, W, relTolRefine, maxRefineIts, progress );
      };

    return FGMRES( applyA, precond, B, relTol, restart, maxIts, progress );
}

// TODO(poulson): Add RGMRES

template<typename Field>
Int SolveAfter
( const SparseMatrix<Field>& A,
  const Matrix<Base<Field>>& reg,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  const RegSolveCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    switch( ctrl.alg )
    {
    case REG_SOLVE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, sparseLDLFact, B, 
          ctrl.relTol,
          ctrl.restart,
          ctrl.maxIts,
          ctrl.relTolRefine,
          ctrl.maxRefineIts, 
          ctrl.progress,
          ctrl.time );
    case REG_SOLVE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, sparseLDLFact, B, 
          ctrl.relTol,
          ctrl.restart,
          ctrl.maxIts,
          ctrl.relTolRefine,
          ctrl.maxRefineIts, 
          ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
        return -1;
    }
}

template<typename Field>
Int SolveAfter
( const SparseMatrix<Field>& A, 
  const Matrix<Base<Field>>& reg,
  const Matrix<Base<Field>>& d,
  const SparseLDLFactorization<Field>& sparseLDLFact,
        Matrix<Field>& B,
  const RegSolveCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    switch( ctrl.alg )
    {
    case REG_SOLVE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, d, sparseLDLFact, B, 
          ctrl.relTol,
          ctrl.restart,
          ctrl.maxIts,
          ctrl.relTolRefine,
          ctrl.maxRefineIts, 
          ctrl.progress,
          ctrl.time );
    case REG_SOLVE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, d, sparseLDLFact, B, 
          ctrl.relTol,
          ctrl.restart,
          ctrl.maxIts,
          ctrl.relTolRefine,
          ctrl.maxRefineIts, 
          ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
        return -1;
    }
}

template<typename Field>
Int SolveAfter
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  const RegSolveCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    switch( ctrl.alg )
    {
    case REG_SOLVE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, sparseLDLFact, B,
          ctrl.relTol,
          ctrl.restart,
          ctrl.maxIts,
          ctrl.relTolRefine,
          ctrl.maxRefineIts, 
          ctrl.progress,
          ctrl.time );
    case REG_SOLVE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, sparseLDLFact, B,
          ctrl.relTol,
          ctrl.restart,
          ctrl.maxIts,
          ctrl.relTolRefine,
          ctrl.maxRefineIts, 
          ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
        return -1;
    }
}

template<typename Field>
Int SolveAfter
( const DistSparseMatrix<Field>& A, 
  const DistMultiVec<Base<Field>>& reg,
  const DistMultiVec<Base<Field>>& d,
  const DistSparseLDLFactorization<Field>& sparseLDLFact,
        DistMultiVec<Field>& B,
  const RegSolveCtrl<Base<Field>>& ctrl )
{
    EL_DEBUG_CSE
    switch( ctrl.alg )
    {
    case REG_SOLVE_FGMRES:
        return FGMRESSolveAfter
        ( A, reg, d, sparseLDLFact, B,
          ctrl.relTol,
          ctrl.restart,
          ctrl.maxIts,
          ctrl.relTolRefine,
          ctrl.maxRefineIts, 
          ctrl.progress,
          ctrl.time );
    case REG_SOLVE_LGMRES:
        return LGMRESSolveAfter
        ( A, reg, d, sparseLDLFact, B,
          ctrl.relTol,
          ctrl.restart,
          ctrl.maxIts,
          ctrl.relTolRefine,
          ctrl.maxRefineIts, 
          ctrl.progress );
    default:
        LogicError("Invalid refinement algorithm");
        return -1;
    }
}

#define PROTO(Field) \
  template Int RegularizedSolveAfter \
  ( const SparseMatrix<Field>& A, \
    const Matrix<Base<Field>>& reg, \
    const SparseLDLFactorization<Field>& sparseLDLFact, \
          Matrix<Field>& B, \
    Base<Field> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int RegularizedSolveAfter \
  ( const SparseMatrix<Field>& A, \
    const Matrix<Base<Field>>& reg, \
    const Matrix<Base<Field>>& d, \
    const SparseLDLFactorization<Field>& sparseLDLFact, \
          Matrix<Field>& B, \
    Base<Field> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int RegularizedSolveAfter \
  ( const DistSparseMatrix<Field>& A, \
    const DistMultiVec<Base<Field>>& reg, \
    const DistSparseLDLFactorization<Field>& sparseLDLFact, \
          DistMultiVec<Field>& B, \
    Base<Field> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int RegularizedSolveAfter \
  ( const DistSparseMatrix<Field>& A, \
    const DistMultiVec<Base<Field>>& reg, \
    const DistMultiVec<Base<Field>>& d, \
    const DistSparseLDLFactorization<Field>& sparseLDLFact, \
          DistMultiVec<Field>& B, \
    Base<Field> relTol, Int maxRefineIts, bool progress, bool time ); \
  template Int SolveAfter \
  ( const SparseMatrix<Field>& A, \
    const Matrix<Base<Field>>& reg, \
    const SparseLDLFactorization<Field>& sparseLDLFact, \
          Matrix<Field>& B, \
    const RegSolveCtrl<Base<Field>>& ctrl ); \
  template Int SolveAfter \
  ( const SparseMatrix<Field>& A, \
    const Matrix<Base<Field>>& reg, \
    const Matrix<Base<Field>>& d, \
    const SparseLDLFactorization<Field>& sparseLDLFact, \
          Matrix<Field>& B, \
    const RegSolveCtrl<Base<Field>>& ctrl ); \
  template Int SolveAfter \
  ( const DistSparseMatrix<Field>& A, \
    const DistMultiVec<Base<Field>>& reg, \
    const DistSparseLDLFactorization<Field>& sparseLDLFact, \
          DistMultiVec<Field>& B, \
    const RegSolveCtrl<Base<Field>>& ctrl ); \
  template Int SolveAfter \
  ( const DistSparseMatrix<Field>& A, \
    const DistMultiVec<Base<Field>>& reg, \
    const DistMultiVec<Base<Field>>& d, \
    const DistSparseLDLFactorization<Field>& sparseLDLFact, \
          DistMultiVec<Field>& B, \
    const RegSolveCtrl<Base<Field>>& ctrl );

#define EL_NO_INT_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace reg_ldl
} // namespace El
