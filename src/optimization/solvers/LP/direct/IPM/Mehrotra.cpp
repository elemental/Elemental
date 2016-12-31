/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
#include "./util.hpp"

namespace El {

template<typename Real>
void CopyOrViewHelper( const DistMatrix<Real>& A, DistMatrix<Real>& B )
{
    EL_DEBUG_CSE
    if( A.ColAlign() == B.ColAlign() && A.RowAlign() == B.RowAlign() )
        LockedView( B, A );
    else
        B = A;
}

// TODO(poulson): Move this into a higher-level namespace

template<typename Ring>
void Gemv
( Orientation orientation,
  const Ring& alpha,
  const SparseMatrix<Ring>& A,
  const Matrix<Ring>& x,
  const Ring& beta,
        Matrix<Ring>& y )
{
    EL_DEBUG_CSE
    Multiply( orientation, alpha, A, x, beta, y );
}

template<typename Ring>
void Gemv
( Orientation orientation,
  const Ring& alpha,
  const DistSparseMatrix<Ring>& A,
  const DistMultiVec<Ring>& x,
  const Ring& beta,
        DistMultiVec<Ring>& y )
{
    EL_DEBUG_CSE
    Multiply( orientation, alpha, A, x, beta, y );
}

namespace lp {
namespace direct {

// The following solves the pair of linear programs in "direct" conic form:
//
//   min c^T x
//   s.t. A x = b, x >= 0,
//
//   max -b^T y
//   s.t. A^T y - z + c = 0, z >= 0,
//
// as opposed to the more general "affine" conic form:
//
//   min c^T x
//   s.t. A x = b, G x + s = h, s >= 0,
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z >= 0
//
// using a Mehrotra Predictor-Corrector scheme.
//

// TODO(poulson): Experiment with using the lagged factorization as a
// preconditioner.

template<typename Real>
struct DenseDirectLPEquilibration
{
    Real xScale;
    Real zScale;
    Matrix<Real> rowScale;
    Matrix<Real> colScale;
};
template<typename Real>
struct DistDenseDirectLPEquilibration
{
    Real xScale;
    Real zScale;
    DistMatrix<Real,MC,STAR> rowScale;
    DistMatrix<Real,MR,STAR> colScale;
};
template<typename Real>
struct SparseDirectLPEquilibration
{
    Real xScale;
    Real zScale;
    Matrix<Real> rowScale;
    Matrix<Real> colScale;
};
template<typename Real>
struct DistSparseDirectLPEquilibration
{
    Real xScale;
    Real zScale;
    DistMultiVec<Real> rowScale;
    DistMultiVec<Real> colScale;
};

template<typename Real>
void Equilibrate
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const DirectLPSolution<Matrix<Real>>& solution,
        DirectLPProblem<Matrix<Real>,Matrix<Real>>& equilibratedProblem,
        DirectLPSolution<Matrix<Real>>& equilibratedSolution,
        DenseDirectLPEquilibration<Real>& equilibration,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    equilibratedProblem = problem;
    equilibratedSolution = solution;

    RuizEquil
    ( equilibratedProblem.A, equilibration.rowScale, equilibration.colScale,
      ctrl.print );

    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScale, equilibratedProblem.b );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.colScale, equilibratedProblem.c );
    if( ctrl.primalInit )
        DiagonalScale
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.x );
    if( ctrl.dualInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScale, equilibratedSolution.y );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.z );
    }

    // Rescale || b ||_max to roughly one.
    equilibration.xScale = Max(MaxNorm(equilibratedProblem.b),Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.xScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.xScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.zScale = Max(MaxNorm(equilibratedProblem.c),Real(1));
    equilibratedProblem.c *= Real(1)/equilibration.zScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.zScale;
        equilibratedSolution.z *= Real(1)/equilibration.zScale;
    }
}

template<typename Real>
void Equilibrate
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const DirectLPSolution<DistMatrix<Real>>& solution,
        DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& equilibratedProblem,
        DirectLPSolution<DistMatrix<Real>>& equilibratedSolution,
        DistDenseDirectLPEquilibration<Real>& equilibration,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = problem.A.Grid();
    ForceSimpleAlignments( equilibratedProblem, grid );
    ForceSimpleAlignments( equilibratedSolution, grid );
    equilibratedProblem = problem;
    equilibratedSolution = solution;
    RuizEquil
    ( equilibratedProblem.A,
      equilibration.rowScale, equilibration.colScale, ctrl.print );

    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScale, equilibratedProblem.b );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.colScale, equilibratedProblem.c );
    if( ctrl.primalInit )
        DiagonalScale
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.x );
    if( ctrl.dualInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScale, equilibratedSolution.y );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.z );
    }

    // Rescale || b ||_max to roughly one.
    equilibration.xScale = Max(MaxNorm(equilibratedProblem.b),Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.xScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.xScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.zScale = Max(MaxNorm(equilibratedProblem.c),Real(1));
    equilibratedProblem.c *= Real(1)/equilibration.zScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.zScale;
        equilibratedSolution.z *= Real(1)/equilibration.zScale;
    }
}

template<typename Real>
void Equilibrate
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const DirectLPSolution<Matrix<Real>>& solution,
        DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& equilibratedProblem,
        DirectLPSolution<Matrix<Real>>& equilibratedSolution,
        SparseDirectLPEquilibration<Real>& equilibration,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    equilibratedProblem = problem;
    equilibratedSolution = solution;

    RuizEquil
    ( equilibratedProblem.A,
      equilibration.rowScale, equilibration.colScale, ctrl.print );

    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScale, equilibratedProblem.b );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.colScale, equilibratedProblem.c );
    if( ctrl.primalInit )
        DiagonalScale
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.x );
    if( ctrl.dualInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScale, equilibratedSolution.y );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.z );
    }

    // Rescale || b ||_max to roughly one.
    equilibration.xScale = Max(MaxNorm(equilibratedProblem.b),Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.xScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.xScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.zScale = Max(MaxNorm(equilibratedProblem.c),Real(1));
    equilibratedProblem.c *= Real(1)/equilibration.zScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.zScale;
        equilibratedSolution.z *= Real(1)/equilibration.zScale;
    }
}

template<typename Real>
void Equilibrate
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const DirectLPSolution<DistMultiVec<Real>>& solution,
        DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>&
          equilibratedProblem,
        DirectLPSolution<DistMultiVec<Real>>& equilibratedSolution,
        DistSparseDirectLPEquilibration<Real>& equilibration,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = problem.A.Grid();
    ForceSimpleAlignments( equilibratedProblem, grid );
    ForceSimpleAlignments( equilibratedSolution, grid );

    equilibratedProblem = problem;
    equilibratedSolution = solution;
    equilibration.rowScale.SetGrid( grid );
    equilibration.colScale.SetGrid( grid );
    RuizEquil
    ( equilibratedProblem.A,
      equilibration.rowScale, equilibration.colScale, ctrl.print );

    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScale, equilibratedProblem.b );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.colScale, equilibratedProblem.c );
    if( ctrl.primalInit )
        DiagonalScale
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.x );
    if( ctrl.dualInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScale, equilibratedSolution.y );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.z );
    }

    // Rescale || b ||_max to roughly one.
    equilibration.xScale = Max(MaxNorm(equilibratedProblem.b),Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.xScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.xScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.zScale = Max(MaxNorm(equilibratedProblem.c),Real(1));
    equilibratedProblem.c *= Real(1)/equilibration.zScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.zScale;
        equilibratedSolution.z *= Real(1)/equilibration.zScale;
    }
}

template<typename Real>
void UndoEquilibration
( const DirectLPSolution<Matrix<Real>>& equilibratedSolution,
  const DenseDirectLPEquilibration<Real>& equilibration,
        DirectLPSolution<Matrix<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;
    solution.x *= equilibration.xScale;
    solution.y *= equilibration.zScale;
    solution.z *= equilibration.zScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.colScale, solution.x );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScale, solution.y );
    DiagonalScale( LEFT, NORMAL, equilibration.colScale, solution.z );
}

template<typename Real>
void UndoEquilibration
( const DirectLPSolution<DistMatrix<Real>>& equilibratedSolution,
  const DistDenseDirectLPEquilibration<Real>& equilibration,
        DirectLPSolution<DistMatrix<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;
    solution.x *= equilibration.xScale;
    solution.y *= equilibration.zScale;
    solution.z *= equilibration.zScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.colScale, solution.x );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScale, solution.y );
    DiagonalScale( LEFT, NORMAL, equilibration.colScale, solution.z );
}

template<typename Real>
void UndoEquilibration
( const DirectLPSolution<Matrix<Real>>& equilibratedSolution,
  const SparseDirectLPEquilibration<Real>& equilibration,
        DirectLPSolution<Matrix<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;
    solution.x *= equilibration.xScale;
    solution.y *= equilibration.zScale;
    solution.z *= equilibration.zScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.colScale, solution.x );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScale, solution.y );
    DiagonalScale( LEFT, NORMAL, equilibration.colScale, solution.z );
}

template<typename Real>
void UndoEquilibration
( const DirectLPSolution<DistMultiVec<Real>>& equilibratedSolution,
  const DistSparseDirectLPEquilibration<Real>& equilibration,
        DirectLPSolution<DistMultiVec<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;
    solution.x *= equilibration.xScale;
    solution.y *= equilibration.zScale;
    solution.z *= equilibration.zScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.colScale, solution.x );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScale, solution.y );
    DiagonalScale( LEFT, NORMAL, equilibration.colScale, solution.z );
}

template<typename Real>
struct DirectRegularization
{
    Real primalEquality;
    Real dualEquality;
};

template<typename Real,typename MatrixType,typename VectorType>
struct DirectState;

template<typename Real>
struct DirectState<Real,Matrix<Real>,Matrix<Real>>
{
    Real cNorm;
    Real bNorm;

    Real barrier;
    Real barrierOld;
    Real barrierAffine;
    Real sigma;

    Real primalObjective;
    Real dualObjective;
    Real relativeGap;

    DirectLPResidual<Matrix<Real>> residual;
    Real primalEqualityNorm;
    Real dualEqualityNorm;
    Real dualConicNorm;
    Real relativePrimalEqualityNorm;
    Real relativeDualEqualityNorm;

    Int numIts;
    Real dimacsError;

    void Initialize
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const MehrotraCtrl<Real>& ctrl );

    void Update
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const DirectLPSolution<Matrix<Real>>& solution,
      const DirectRegularization<Real>& permReg,
      const MehrotraCtrl<Real>& ctrl );

    void PrintResiduals
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const DirectLPSolution<Matrix<Real>>& solution,
      const DirectLPSolution<Matrix<Real>>& correction,
      const DirectRegularization<Real>& permReg ) const;
};

template<typename Real>
void DirectState<Real,Matrix<Real>,Matrix<Real>>::Initialize
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    bNorm = FrobeniusNorm( problem.b );
    cNorm = FrobeniusNorm( problem.c );
    barrierOld = 0.1;
    dimacsError = 1;
    if( ctrl.print )
    {
        const Real ANrm1 = OneNorm( problem.A );
        Output("|| A ||_1 = ",ANrm1);
        Output("|| b ||_2 = ",bNorm);
        Output("|| c ||_2 = ",cNorm);
    }
}

template<typename Real>
void DirectState<Real,Matrix<Real>,Matrix<Real>>::Update
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const DirectLPSolution<Matrix<Real>>& solution,
  const DirectRegularization<Real>& permReg,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int degree = problem.A.Width();

    // Compute the new barrier parameter
    // ---------------------------------
    barrier = Dot(solution.x,solution.z) / degree;
    const Real compRatio = pos_orth::ComplementRatio( solution.x, solution.z );
    barrier = compRatio > ctrl.balanceTol ? barrierOld
      : Min(barrier,barrierOld);
    barrierOld = barrier;

    // Compute the objectives and relative duality gap
    primalObjective = Dot(problem.c,solution.x);
    dualObjective = -Dot(problem.b,solution.y);
    relativeGap = Abs(primalObjective-dualObjective) / (1+Abs(primalObjective));

    // Compute the primal equality residual,
    //
    //   r_b = A x - b,
    //
    // and its (relative) norm.
    //
    residual.primalEquality = problem.b;
    Gemv
    ( NORMAL, Real(1), problem.A, solution.x,
      Real(-1), residual.primalEquality );
    primalEqualityNorm = FrobeniusNorm(residual.primalEquality);
    relativePrimalEqualityNorm = primalEqualityNorm / (1 + bNorm);
    Axpy( -permReg.primalEquality, solution.y, residual.primalEquality );

    // Compute the dual equality residual,
    //
    //   r_c = A^T y - z + c,
    //
    // and its (relative) norm.
    //
    residual.dualEquality = problem.c;
    Gemv
    ( TRANSPOSE, Real(1), problem.A, solution.y,
      Real(1), residual.dualEquality );
    residual.dualEquality -= solution.z;
    dualEqualityNorm = FrobeniusNorm(residual.dualEquality);
    relativeDualEqualityNorm = dualEqualityNorm / (1 + cNorm);
    Axpy( permReg.dualEquality, solution.x, residual.dualEquality );

    // Compute the complimentarity vector,
    //
    //   r_mu := x o z,
    //
    // and its norm.
    //
    residual.dualConic = solution.z;
    DiagonalScale( LEFT, NORMAL, solution.x, residual.dualConic );
    dualConicNorm = FrobeniusNorm(residual.dualConic);

    // Now check the pieces
    // --------------------
    dimacsError =
      Max(Max(relativePrimalEqualityNorm,relativeDualEqualityNorm),relativeGap);
    if( ctrl.print )
    {
        const Real xNrm2 = FrobeniusNorm( solution.x );
        const Real yNrm2 = FrobeniusNorm( solution.y );
        const Real zNrm2 = FrobeniusNorm( solution.z );
        Output
        ("iter ",numIts,":\n",Indent(),
         "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
         "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
         "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
         "  || r_b ||_2 = ",primalEqualityNorm,"\n",Indent(),
         "  || r_c ||_2 = ",dualEqualityNorm,"\n",Indent(),
         "  || r_b ||_2 / (1 + || b ||_2) = ",relativePrimalEqualityNorm,
         "\n",Indent(),
         "  || r_c ||_2 / (1 + || c ||_2) = ",relativeDualEqualityNorm,
         "\n",Indent(),
         "  primal = ",primalObjective,"\n",Indent(),
         "  dual   = ",dualObjective,"\n",Indent(),
         "  |primal - dual| / (1 + |primal|) = ",relativeGap,
         "\n",Indent(),
         "  DIMACS: ",dimacsError);
    }
}

template<typename Real>
void DirectState<Real,Matrix<Real>,Matrix<Real>>::PrintResiduals
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const DirectLPSolution<Matrix<Real>>& solution,
  const DirectLPSolution<Matrix<Real>>& correction,
  const DirectRegularization<Real>& permReg ) const
{
    EL_DEBUG_CSE
    DirectLPResidual<Matrix<Real>> error;
    Matrix<Real> prod;

    error.primalEquality = residual.primalEquality;
    Gemv
    ( NORMAL, Real(1), problem.A, correction.x,
      Real(1), error.primalEquality );
    Axpy
    ( -permReg.primalEquality, correction.y, error.primalEquality );
    Real dxErrorNrm2 = FrobeniusNorm( error.primalEquality );

    error.dualEquality = residual.dualEquality;
    Gemv
    ( TRANSPOSE, Real(1), problem.A, correction.y,
      Real(1), error.dualEquality );
    Axpy( permReg.dualEquality, correction.x, error.dualEquality );
    error.dualEquality -= correction.z;
    Real dyErrorNrm2 = FrobeniusNorm( error.dualEquality );

    error.dualConic = residual.dualConic;
    prod = correction.z;
    DiagonalScale( LEFT, NORMAL, solution.x, prod );
    error.dualConic += prod;
    prod = correction.x;
    DiagonalScale( LEFT, NORMAL, solution.z, prod );
    error.dualConic += prod;
    Real dzErrorNrm2 = FrobeniusNorm( error.dualConic );

    Output
    ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
     dxErrorNrm2/(1+primalEqualityNorm),"\n",Indent(),
     "|| dyError ||_2 / (1 + || r_c ||_2) = ",
     dyErrorNrm2/(1+dualEqualityNorm),"\n",Indent(),
     "|| dzError ||_2 / (1 + || r_h ||_2) = ",
     dzErrorNrm2/(1+dualConicNorm));
}

template<typename Real,class MatrixType,class VectorType>
struct DirectKKTSolver;

template<typename Real>
struct DirectKKTSolver<Real,Matrix<Real>,Matrix<Real>>
{
    Matrix<Real> kktSystem;
    Matrix<Real> dSub;
    Permutation perm;

    mutable Matrix<Real> temp;

    void Factor()
    { LDL( kktSystem, dSub, perm, false ); }

    void Solve( Matrix<Real>& rhs ) const
    { ldl::SolveAfter( kktSystem, dSub, perm, rhs, false ); }

    void InitializeSystem
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const DirectRegularization<Real>& permReg,
      const DirectLPSolution<Matrix<Real>>& solution,
      KKTSystem system )
    {
        if( system == FULL_KKT )
        {
            KKT( problem.A, solution.x, solution.z, kktSystem );
        }
        else if( system == AUGMENTED_KKT )
        {
            AugmentedKKT( problem.A, solution.x, solution.z, kktSystem );
        }
        else if( system == NORMAL_KKT )
        {
            NormalKKT
            ( problem.A, Sqrt(permReg.dualEquality),
              Sqrt(permReg.primalEquality), solution.x, solution.z, kktSystem );
        }
        Factor();
    }

    void SolveSystem
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const DirectRegularization<Real>& permReg,
      const DirectLPResidual<Matrix<Real>>& residual,
      const DirectLPSolution<Matrix<Real>>& solution,
            DirectLPSolution<Matrix<Real>>& correction,
            KKTSystem system ) const
    {
        const Int m = solution.y.Height();
        const Int n = solution.x.Height();
        if( system == FULL_KKT )
        {
            KKTRHS
            ( residual.dualEquality,
              residual.primalEquality,
              residual.dualConic,
              solution.z, temp );
            Solve(temp);
            ExpandSolution
            ( m, n, temp, correction.x, correction.y, correction.z );
        }
        else if( system == AUGMENTED_KKT )
        {
            AugmentedKKTRHS
            ( solution.x, residual.dualEquality,
              residual.primalEquality, residual.dualConic, temp );
            Solve(temp);
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, temp,
              correction.x, correction.y, correction.z );
        }
        else if( system == NORMAL_KKT )
        {
            NormalKKTRHS
            ( problem.A, Sqrt(permReg.dualEquality), solution.x, solution.z,
              residual.dualEquality,
              residual.primalEquality,
              residual.dualConic,
              correction.y );
            Solve(correction.y);
            ExpandNormalSolution
            ( problem, Sqrt(permReg.dualEquality), solution, residual,
              correction );
        }
    }
};

template<typename Real>
void EquilibratedMehrotra
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl,
  bool outputRoot )
{
    EL_DEBUG_CSE
    const Int n = problem.A.Width();
    const Int degree = n;

    // TODO(poulson): Implement nonzero regularization
    DirectRegularization<Real> permReg;
    permReg.primalEquality = 0;
    permReg.dualEquality = 0;

    DirectState<Real,Matrix<Real>,Matrix<Real>> state;
    state.Initialize( problem, ctrl );

    const Int indent = PushIndent();
    try {
    Initialize
    ( problem, solution,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift );
    DirectKKTSolver<Real,Matrix<Real>,Matrix<Real>> solver;
    DirectLPSolution<Matrix<Real>> affineCorrection, correction;
    for( state.numIts=0; state.numIts<ctrl.maxIts; ++state.numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = pos_orth::NumOutside( solution.x );
        const Int zNumNonPos = pos_orth::NumOutside( solution.z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        state.Update( problem, solution, permReg, ctrl );

        // Check for convergence
        // =====================
        if( state.dimacsError <= ctrl.targetTol )
            break;

        // Set up and factor the KKT matrix
        // ================================
        solver.InitializeSystem( problem, permReg, solution, ctrl.system );

        // Compute the affine search direction
        // ===================================
        solver.SolveSystem
        ( problem, permReg, state.residual, solution, affineCorrection,
          ctrl.system );
        if( ctrl.checkResiduals && ctrl.print )
        {
            state.PrintResiduals
            ( problem, solution, affineCorrection, permReg );
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri =
          pos_orth::MaxStep( solution.x, affineCorrection.x, Real(1) );
        Real alphaAffDual =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && outputRoot )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: correction.z and correction.x are used as temporaries
        correction.x = solution.x;
        correction.z = solution.z;
        Axpy( alphaAffPri,  affineCorrection.x, correction.x );
        Axpy( alphaAffDual, affineCorrection.z, correction.z );
        state.barrierAffine = Dot(correction.x,correction.z) / degree;
        if( ctrl.print && outputRoot )
            Output
            ("barrierAffine = ",state.barrierAffine,", barrier=",state.barrier);
        state.sigma =
          ctrl.centralityRule
          (state.barrier,state.barrierAffine,alphaAffPri,alphaAffDual);
        if( ctrl.print && outputRoot )
            Output("sigma=",state.sigma);

        // Solve for the combined direction
        // ================================
        state.residual.primalEquality *= 1-state.sigma;
        state.residual.dualEquality *= 1-state.sigma;
        Shift( state.residual.dualConic, -state.sigma*state.barrier );
        if( ctrl.mehrotra )
        {
            // r_mu += dxAff o dzAff
            // ---------------------
            // NOTE: We are using correction.z as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.x, correction.z );
            state.residual.dualConic += correction.z;
        }
        solver.SolveSystem
        ( problem, permReg, state.residual, solution, correction,
          ctrl.system );
        if( ctrl.checkResiduals && ctrl.print )
        {
            state.PrintResiduals( problem, solution, correction, permReg );
        }

        // Update the current estimates
        // ============================
        Real alphaPri =
          pos_orth::MaxStep( solution.x, correction.x, 1/ctrl.maxStepRatio );
        Real alphaDual =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && outputRoot )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaDual, correction.y, solution.y );
        Axpy( alphaDual, correction.z, solution.z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
            RuntimeError("Zero step size");
    }
    } catch(...) { }
    if( state.dimacsError > ctrl.minTol )
        RuntimeError
        ("Unable to achieve minimum tolerance ",ctrl.minTol);

    SetIndent( indent );
}

template<typename Real>
void Mehrotra
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const bool outputRoot = true;
    if( ctrl.outerEquil )
    {
        DirectLPProblem<Matrix<Real>,Matrix<Real>> equilibratedProblem;
        DirectLPSolution<Matrix<Real>> equilibratedSolution;
        DenseDirectLPEquilibration<Real> equilibration;
        Equilibrate
        ( problem, solution,
          equilibratedProblem, equilibratedSolution, equilibration, ctrl );
        EquilibratedMehrotra
        ( equilibratedProblem, equilibratedSolution, ctrl, outputRoot );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedMehrotra( problem, solution, ctrl, outputRoot );
    }
    if( ctrl.print && outputRoot )
    {
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        const Real xNrm2 = FrobeniusNorm( solution.x );
        const Real yNrm2 = FrobeniusNorm( solution.y );
        const Real zNrm2 = FrobeniusNorm( solution.z );
        Output
        ("Exiting with:\n",Indent(),
         "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
         "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
         "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
         "  primal = ",primObj,"\n",Indent(),
         "  dual   = ",dualObj,"\n",Indent(),
         "  |primal - dual| / (1 + |primal|) = ",objConv);
    }
}

// This interface is now deprecated.
template<typename Real>
void Mehrotra
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    DirectLPProblem<Matrix<Real>,Matrix<Real>> problem;
    DirectLPSolution<Matrix<Real>> solution;
    LockedView( problem.c, c );
    LockedView( problem.A, A );
    LockedView( problem.b, b );
    solution.x = x;
    solution.y = y;
    solution.z = z;
    Mehrotra( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

template<typename Real>
void EquilibratedMehrotra
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int degree = n;
    const Grid& grid = problem.A.Grid();
    const int commRank = grid.Rank();

    // TODO(poulson): Implement nonzero regularization
    const Real gammaPerm = 0;
    const Real deltaPerm = 0;

    const Real bNrm2 = FrobeniusNorm( problem.b );
    const Real cNrm2 = FrobeniusNorm( problem.c );
    if( ctrl.print )
    {
        const Real ANrm1 = OneNorm( problem.A );
        if( commRank == 0 )
        {
            Output("|| A ||_1 = ",ANrm1);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
        }
    }

    Initialize
    ( problem, solution,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift );

    Real muOld = 0.1;
    Real relError = 1;
    DistMatrix<Real> J(grid), d(grid);
    DistMatrix<Real> dSub(grid);
    DistPermutation p(grid);
    auto attemptToFactor = [&]()
      {
        try { LDL( J, dSub, p, false ); }
        catch(...)
        {
            if( relError > ctrl.minTol )
                RuntimeError
                ("Unable to achieve minimum tolerance ",ctrl.minTol);
            return false;
        }
        return true;
      };
    auto attemptToSolve = [&]( DistMatrix<Real>& rhs )
      {
        try { ldl::SolveAfter( J, dSub, p, rhs, false ); }
        catch(...)
        {
            if( relError > ctrl.minTol )
                RuntimeError
                ("Unable to achieve minimum tolerance ",ctrl.minTol);
            return false;
        }
        return true;
      };

    DirectLPSolution<DistMatrix<Real>> affineCorrection, correction;
    ForceSimpleAlignments( affineCorrection, grid );
    ForceSimpleAlignments( correction, grid );

    DirectLPResidual<DistMatrix<Real>> residual, error;
    ForceSimpleAlignments( residual, grid );
    ForceSimpleAlignments( error, grid );

    DistMatrix<Real> prod(grid);
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = pos_orth::NumOutside( solution.x );
        const Int zNumNonPos = pos_orth::NumOutside( solution.z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the barrier parameter
        // =============================
        Real mu = Dot(solution.x,solution.z) / degree;
        const Real compRatio =
          pos_orth::ComplementRatio( solution.x, solution.z );
        mu = compRatio > ctrl.balanceTol ? muOld : Min(mu,muOld);
        muOld = mu;

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        residual.primalEquality = problem.b;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(-1), residual.primalEquality );
        const Real rbNrm2 = FrobeniusNorm( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        Axpy( -deltaPerm*deltaPerm, solution.y, residual.primalEquality );
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        residual.dualEquality -= solution.z;
        const Real rcNrm2 = FrobeniusNorm( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        Axpy( gammaPerm*gammaPerm, solution.x, residual.dualEquality );
        // Now check the pieces
        // --------------------
        relError = Max(Max(objConv,rbConv),rcConv);
        if( ctrl.print )
        {
            const Real xNrm2 = FrobeniusNorm( solution.x );
            const Real yNrm2 = FrobeniusNorm( solution.y );
            const Real zNrm2 = FrobeniusNorm( solution.z );
            if( commRank == 0 )
                Output
                ("iter ",numIts,":\n",Indent(),
                 "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
                 "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
                 "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
                 "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
                 "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
                 "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
                 "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
                 "  primal = ",primObj,"\n",Indent(),
                 "  dual   = ",dualObj,"\n",Indent(),
                 "  |primal - dual| / (1 + |primal|) = ",objConv);
        }
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := x o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.x, residual.dualConic );

        if( ctrl.system == FULL_KKT )
        {
            // Construct the KKT system
            // ------------------------
            KKT( problem.A, solution.x, solution.z, J );
            KKTRHS
            ( residual.dualEquality, residual.primalEquality,
              residual.dualConic, solution.z, d );

            // Solve for the direction
            // -----------------------
            if( !attemptToFactor() )
                break;
            if( !attemptToSolve(d) )
                break;
            ExpandSolution
            ( m, n, d,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the KKT system
            // ------------------------
            AugmentedKKT( problem.A, solution.x, solution.z, J );
            AugmentedKKTRHS
            ( solution.x, residual.dualEquality, residual.primalEquality,
              residual.dualConic, d );

            // Solve for the direction
            // -----------------------
            if( !attemptToFactor() )
                break;
            if( !attemptToSolve(d) )
                break;
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, d,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }
        else if( ctrl.system == NORMAL_KKT )
        {
            // Construct the KKT system
            // ------------------------
            NormalKKT
            ( problem.A, gammaPerm, deltaPerm, solution.x, solution.z, J );
            NormalKKTRHS
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, affineCorrection.y );

            // Solve for the direction
            // -----------------------
            if( !attemptToFactor() )
                break;
            if( !attemptToSolve(affineCorrection.y) )
                break;
            ExpandNormalSolution
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }

        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Gemv
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            Axpy
            ( -deltaPerm*deltaPerm, affineCorrection.y, error.primalEquality );
            Real dxErrorNrm2 = FrobeniusNorm( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Gemv
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Axpy( gammaPerm*gammaPerm, affineCorrection.x, error.dualEquality );
            error.dualEquality -= affineCorrection.z;
            Real dyErrorNrm2 = FrobeniusNorm( error.dualEquality );

            Real rmuNrm2 = FrobeniusNorm( residual.dualConic );
            error.dualConic = residual.dualConic;
            prod = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, solution.x, prod );
            error.dualConic += prod;
            prod = affineCorrection.x;
            DiagonalScale( LEFT, NORMAL, solution.z, prod );
            error.dualConic += prod;
            Real dzErrorNrm2 = FrobeniusNorm( error.dualConic );

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri =
          pos_orth::MaxStep( solution.x, affineCorrection.x, Real(1) );
        Real alphaAffDual =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: correction.z and correction.x are used as temporaries
        correction.x = solution.x;
        correction.z = solution.z;
        Axpy( alphaAffPri,  affineCorrection.x, correction.x );
        Axpy( alphaAffDual, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.x,correction.z) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma =
          ctrl.centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        residual.primalEquality *= 1-sigma;
        residual.dualEquality *= 1-sigma;
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dxAff o dzAff
            // ---------------------
            // NOTE: correction.z is used as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.x, correction.z );
            residual.dualConic += correction.z;
        }

        if( ctrl.system == FULL_KKT )
        {
            // Construct the new KKT RHS
            // -------------------------
            KKTRHS
            ( residual.dualEquality, residual.primalEquality,
              residual.dualConic, solution.z, d );

            // Solve for the direction
            // -----------------------
            if( !attemptToSolve(d) )
                break;
            ExpandSolution( m, n, d, correction.x, correction.y, correction.z );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the new KKT RHS
            // -------------------------
            AugmentedKKTRHS
            ( solution.x, residual.dualEquality, residual.primalEquality,
              residual.dualConic, d );

            // Solve for the direction
            // -----------------------
            if( !attemptToSolve(d) )
                break;
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, d,
              correction.x, correction.y, correction.z );
        }
        else if( ctrl.system == NORMAL_KKT )
        {
            // Construct the new KKT RHS
            // -------------------------
            NormalKKTRHS
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, correction.y );

            // Solve for the direction
            // -----------------------
            if( !attemptToSolve(correction.y) )
                break;
            ExpandNormalSolution
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              correction.x, correction.y, correction.z );
        }
        // TODO(poulson): Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri =
          pos_orth::MaxStep( solution.x, correction.x, 1/ctrl.maxStepRatio );
        Real alphaDual =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaDual, correction.y, solution.y );
        Axpy( alphaDual, correction.z, solution.z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
    }
    SetIndent( indent );
}

template<typename Real>
void Mehrotra
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = problem.A.Grid();
    if( ctrl.outerEquil )
    {
        DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>> equilibratedProblem;
        DirectLPSolution<DistMatrix<Real>> equilibratedSolution;
        DistDenseDirectLPEquilibration<Real> equilibration;
        ForceSimpleAlignments( equilibratedProblem, grid );
        ForceSimpleAlignments( equilibratedSolution, grid );
        Equilibrate
        ( problem, solution,
          equilibratedProblem, equilibratedSolution,
          equilibration, ctrl );
        EquilibratedMehrotra( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        // Avoid creating unnecessary copies where we can.
        if( SimpleAlignments(problem) && SimpleAlignments(solution) )
        {
            EquilibratedMehrotra( problem, solution, ctrl );
        }
        else if( SimpleAlignments(problem) )
        {
            DirectLPSolution<DistMatrix<Real>> alignedSolution;
            ForceSimpleAlignments( alignedSolution, grid );
            alignedSolution = solution;
            EquilibratedMehrotra( problem, alignedSolution, ctrl );
            solution = alignedSolution;
        }
        else if( SimpleAlignments(solution) )
        {
            DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>> alignedProblem;
            ForceSimpleAlignments( alignedProblem, grid );
            CopyOrViewHelper( problem.c, alignedProblem.c );
            CopyOrViewHelper( problem.A, alignedProblem.A );
            CopyOrViewHelper( problem.b, alignedProblem.b );
            EquilibratedMehrotra( alignedProblem, solution, ctrl );
        }
        else
        {
            DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>> alignedProblem;
            ForceSimpleAlignments( alignedProblem, grid );
            CopyOrViewHelper( problem.c, alignedProblem.c );
            CopyOrViewHelper( problem.A, alignedProblem.A );
            CopyOrViewHelper( problem.b, alignedProblem.b );
            DirectLPSolution<DistMatrix<Real>> alignedSolution;
            ForceSimpleAlignments( alignedSolution, grid );
            alignedSolution = solution;
            EquilibratedMehrotra( alignedProblem, alignedSolution, ctrl );
            solution = alignedSolution;
        }
    }
    if( ctrl.print )
    {
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        const Real xNrm2 = FrobeniusNorm( solution.x );
        const Real yNrm2 = FrobeniusNorm( solution.y );
        const Real zNrm2 = FrobeniusNorm( solution.z );
        OutputFromRoot
        (grid.Comm(),
         "Exiting with:\n",Indent(),
         "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
         "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
         "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
         "  primal = ",primObj,"\n",Indent(),
         "  dual   = ",dualObj,"\n",Indent(),
         "  |primal - dual| / (1 + |primal|) = ",objConv);
    }
}

// This interface is now deprecated.
template<typename Real>
void Mehrotra
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = A.Grid();
    DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>> problem;
    DirectLPSolution<DistMatrix<Real>> solution;
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    Copy( c, problem.c );
    Copy( A, problem.A );
    Copy( b, problem.b );
    Copy( x, solution.x );
    Copy( y, solution.y );
    Copy( z, solution.z );
    Mehrotra( problem, solution, ctrl );
    Copy( solution.x, x );
    Copy( solution.y, y );
    Copy( solution.z, z );
}

template<typename Real>
void EquilibratedMehrotra
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int degree = n;

    // TODO(poulson): Move these into the control structure
    Real gammaPerm, deltaPerm, betaPerm, gammaTmp, deltaTmp, betaTmp;
    if( ctrl.system == NORMAL_KKT )
    {
        gammaPerm = deltaPerm = betaPerm = gammaTmp = deltaTmp = betaTmp = 0;
    }
    else
    {
        gammaPerm = ctrl.reg0Perm;
        deltaPerm = ctrl.reg1Perm;
        betaPerm = ctrl.reg2Perm;
        gammaTmp = ctrl.reg0Tmp;
        deltaTmp = ctrl.reg1Tmp;
        betaTmp = ctrl.reg2Tmp;
    }

    const Real bNrm2 = FrobeniusNorm( problem.b );
    const Real cNrm2 = FrobeniusNorm( problem.c );
    const Real twoNormEstA = TwoNormEstimate( problem.A, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstA + 1;
    if( ctrl.print )
    {
        Output("|| A ||_2 estimate: ",twoNormEstA);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
    }

    SparseLDLFactorization<Real> sparseLDLFact;
    // The initialization involves an augmented KKT system, and so we can
    // only reuse the factorization metadata if the this IPM is using the
    // augmented formulation
    if( ctrl.system == AUGMENTED_KKT )
    {
        Initialize
        ( problem, solution, sparseLDLFact,
          ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift,
          ctrl.solveCtrl );
    }
    else
    {
        SparseLDLFactorization<Real> augmentedSparseLDLFact;
        Initialize
        ( problem, solution, augmentedSparseLDLFact,
          ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift,
          ctrl.solveCtrl );
    }

    Matrix<Real> regTmp;
    if( ctrl.system == FULL_KKT )
    {
        regTmp.Resize( m+2*n, 1 );
        for( Int i=0; i<m+2*n; ++i )
        {
            if( i < n )        regTmp(i) = gammaTmp*gammaTmp;
            else if( i < n+m ) regTmp(i) = -deltaTmp*deltaTmp;
            else               regTmp(i) = -betaTmp*betaTmp;
        }
    }
    else if( ctrl.system == AUGMENTED_KKT )
    {
        regTmp.Resize( n+m, 1 );
        for( Int i=0; i<n+m; ++i )
        {
            if( i < n ) regTmp(i) = gammaTmp*gammaTmp;
            else        regTmp(i) = -deltaTmp*deltaTmp;
        }
    }
    else if( ctrl.system == NORMAL_KKT )
    {
        regTmp.Resize( m, 1 );
        Fill( regTmp, deltaTmp*deltaTmp );
    }
    regTmp *= origTwoNormEst;

    Real muOld = 0.1;
    Real relError = 1;
    SparseMatrix<Real> J, JOrig;
    Matrix<Real> d, w;
    Matrix<Real> dInner;

    DirectLPSolution<Matrix<Real>> affineCorrection, correction;
    DirectLPResidual<Matrix<Real>> residual, error;

    Matrix<Real> prod;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = pos_orth::NumOutside( solution.x );
        const Int zNumNonPos = pos_orth::NumOutside( solution.z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        residual.primalEquality = problem.b;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(-1), residual.primalEquality );
        const Real rbNrm2 = FrobeniusNorm( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        Axpy( -deltaPerm*deltaPerm, solution.y, residual.primalEquality );
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        residual.dualEquality -= solution.z;
        const Real rcNrm2 = FrobeniusNorm( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        Axpy( gammaPerm*gammaPerm, solution.x, residual.dualEquality );
        // Now check the pieces
        // --------------------
        relError = Max(Max(objConv,rbConv),rcConv);

        // Compute the scaling point
        // =========================
        pos_orth::NesterovTodd( solution.x, solution.z, w );
        const Real wMaxNorm = MaxNorm( w );

        // Compute the barrier parameter
        // =============================
        Real mu = Dot(solution.x,solution.z) / degree;
        const Real compRatio =
          pos_orth::ComplementRatio( solution.x, solution.z );
        mu = compRatio > ctrl.balanceTol ? muOld : Min(mu,muOld);
        muOld = mu;

        if( ctrl.print )
        {
            const Real xNrm2 = FrobeniusNorm( solution.x );
            const Real yNrm2 = FrobeniusNorm( solution.y );
            const Real zNrm2 = FrobeniusNorm( solution.z );
            Output
            ("iter ",numIts,":\n",Indent(),
             "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
             "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
             "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
             "  ||  w  ||_max = ",wMaxNorm,"\n",Indent(),
             "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
             "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
             "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
             "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
             "  mu        = ",mu,"\n",Indent(),
             "  primal    = ",primObj,"\n",Indent(),
             "  dual      = ",dualObj,"\n",Indent(),
             "  |primal - dual| / (1 + |primal|) = ",objConv);
        }
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := x o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.x, residual.dualConic );

        if( ctrl.system == FULL_KKT || ctrl.system == AUGMENTED_KKT )
        {
            // Construct the KKT system
            // ------------------------
            if( ctrl.system == FULL_KKT )
            {
                KKT
                ( problem.A, gammaPerm, deltaPerm, betaPerm,
                  solution.x, solution.z, JOrig, false );
                KKTRHS
                ( residual.dualEquality, residual.primalEquality,
                  residual.dualConic, solution.z, d );
            }
            else
            {
                AugmentedKKT
                ( problem.A, gammaPerm, deltaPerm, solution.x, solution.z,
                  JOrig, false );
                AugmentedKKTRHS
                ( solution.x, residual.dualEquality, residual.primalEquality,
                  residual.dualConic, d );
            }
            J = JOrig;
            UpdateDiagonal( J, Real(1), regTmp );

            // Solve for the direction
            // -----------------------
            try
            {
                if( wMaxNorm >= ctrl.ruizEquilTol )
                {
                    if( ctrl.print )
                        Output("Running SymmetricRuizEquil");
                    SymmetricRuizEquil
                    ( J, dInner, ctrl.ruizMaxIter, ctrl.print );
                }
                else if( wMaxNorm >= ctrl.diagEquilTol )
                {
                    if( ctrl.print )
                        Output("Running SymmetricDiagonalEquil");
                    SymmetricDiagonalEquil( J, dInner, ctrl.print );
                }
                else
                    Ones( dInner, J.Height(), 1 );

                if( numIts == 0 &&
                    (ctrl.system != AUGMENTED_KKT ||
                     (ctrl.primalInit && ctrl.dualInit)) )
                {
                    const bool hermitian = true;
                    const BisectCtrl bisectCtrl;
                    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );
                }
                else
                {
                    sparseLDLFact.ChangeNonzeroValues( J );
                }

                sparseLDLFact.Factor( LDL_2D );
                if( ctrl.resolveReg )
                    reg_ldl::SolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d, ctrl.solveCtrl );
                else
                    reg_ldl::RegularizedSolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d,
                      ctrl.solveCtrl.relTol,
                      ctrl.solveCtrl.maxRefineIts,
                      ctrl.solveCtrl.progress );
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            if( ctrl.system == FULL_KKT )
                ExpandSolution
                ( m, n, d,
                  affineCorrection.x, affineCorrection.y, affineCorrection.z );
            else
                ExpandAugmentedSolution
                ( solution.x, solution.z, residual.dualConic, d,
                  affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }
        else // ctrl.system == NORMAL_KKT
        {
            // Construct the KKT system
            // ------------------------
            // TODO(poulson): Apply updates to a matrix of explicit zeros
            // (with the correct sparsity pattern)
            NormalKKT
            ( problem.A, gammaPerm, deltaPerm,
              solution.x, solution.z, J, false );
            NormalKKTRHS
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, affineCorrection.y );

            // Solve for the direction
            // -----------------------
            try
            {
                if( numIts == 0 )
                {
                    const bool hermitian = true;
                    const BisectCtrl bisectCtrl;
                    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );
                }
                else
                {
                    sparseLDLFact.ChangeNonzeroValues( J );
                }

                sparseLDLFact.Factor( LDL_2D );

                // NOTE: regTmp should be all zeros; replace with unregularized
                reg_ldl::RegularizedSolveAfter
                ( J, regTmp, sparseLDLFact, affineCorrection.y,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress,
                  ctrl.solveCtrl.time );
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandNormalSolution
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }

        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Gemv
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            Axpy
            ( -deltaPerm*deltaPerm, affineCorrection.y, error.primalEquality );
            Real dxErrorNrm2 = FrobeniusNorm( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Gemv
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Axpy( gammaPerm*gammaPerm, affineCorrection.x, error.dualEquality );
            error.dualEquality -= affineCorrection.z;
            Real dyErrorNrm2 = FrobeniusNorm( error.dualEquality );

            Real rmuNrm2 = FrobeniusNorm( residual.dualConic );
            error.dualConic = residual.dualConic;
            prod = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, solution.x, prod );
            error.dualConic += prod;
            prod = affineCorrection.x;
            DiagonalScale( LEFT, NORMAL, solution.z, prod );
            error.dualConic += prod;
            Real dzErrorNrm2 = FrobeniusNorm( error.dualConic );

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri =
          pos_orth::MaxStep( solution.x, affineCorrection.x, Real(1) );
        Real alphaAffDual =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: correction.z and correction.x are used as temporaries
        correction.x = solution.x;
        correction.z = solution.z;
        Axpy( alphaAffPri,  affineCorrection.x, correction.x );
        Axpy( alphaAffDual, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.x,correction.z) / degree;
        if( ctrl.print )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma =
          ctrl.centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        residual.primalEquality *= 1-sigma;
        residual.dualEquality *= 1-sigma;
        Shift( residual.dualConic, -sigma*mu );
        // TODO(poulson): Gondzio's corrections
        if( ctrl.mehrotra )
        {
            // r_mu += dxAff o dzAff
            // ---------------------
            // NOTE: correction.z is used as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.x, correction.z );
            residual.dualConic += correction.z;
        }

        if( ctrl.system == FULL_KKT )
        {
            KKTRHS
            ( residual.dualEquality, residual.primalEquality,
              residual.dualConic, solution.z, d );
            try
            {
                if( ctrl.resolveReg )
                    reg_ldl::SolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d, ctrl.solveCtrl );
                else
                    reg_ldl::RegularizedSolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d,
                      ctrl.solveCtrl.relTol,
                      ctrl.solveCtrl.maxRefineIts,
                      ctrl.solveCtrl.progress );
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandSolution( m, n, d, correction.x, correction.y, correction.z );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            AugmentedKKTRHS
            ( solution.x, residual.dualEquality, residual.primalEquality,
              residual.dualConic, d );
            try
            {
                if( ctrl.resolveReg )
                    reg_ldl::SolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d, ctrl.solveCtrl );
                else
                    reg_ldl::RegularizedSolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d,
                      ctrl.solveCtrl.relTol,
                      ctrl.solveCtrl.maxRefineIts,
                      ctrl.solveCtrl.progress );
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, d,
              correction.x, correction.y, correction.z );
        }
        else
        {
            NormalKKTRHS
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, correction.y );
            try
            {
                // NOTE: regTmp should be all zeros; replace with unregularized
                reg_ldl::RegularizedSolveAfter
                ( J, regTmp, sparseLDLFact, correction.y,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress,
                  ctrl.solveCtrl.time );
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandNormalSolution
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              correction.x, correction.y, correction.z );
        }
        // TODO(poulson): Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri =
          pos_orth::MaxStep( solution.x, correction.x, 1/ctrl.maxStepRatio );
        Real alphaDual =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaDual, correction.y, solution.y );
        Axpy( alphaDual, correction.z, solution.z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
    }
    SetIndent( indent );
}

template<typename Real>
void Mehrotra
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.outerEquil )
    {
        DirectLPProblem<SparseMatrix<Real>,Matrix<Real>> equilibratedProblem;
        DirectLPSolution<Matrix<Real>> equilibratedSolution;
        SparseDirectLPEquilibration<Real> equilibration;
        Equilibrate
        ( problem, solution,
          equilibratedProblem, equilibratedSolution, equilibration, ctrl );
        EquilibratedMehrotra( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedMehrotra( problem, solution, ctrl );
    }
    if( ctrl.print )
    {
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        const Real xNrm2 = FrobeniusNorm( solution.x );
        const Real yNrm2 = FrobeniusNorm( solution.y );
        const Real zNrm2 = FrobeniusNorm( solution.z );
        Output
        ("Exiting with:\n",Indent(),
         "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
         "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
         "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
         "  primal = ",primObj,"\n",Indent(),
         "  dual   = ",dualObj,"\n",Indent(),
         "  |primal - dual| / (1 + |primal|) = ",objConv);
    }
}

// This interface is now deprecated.
template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    DirectLPProblem<SparseMatrix<Real>,Matrix<Real>> problem;
    DirectLPSolution<Matrix<Real>> solution;
    LockedView( problem.c, c );
    problem.A = A;
    LockedView( problem.b, b );
    solution.x = x;
    solution.y = y;
    solution.z = z;
    Mehrotra( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

// TODO(poulson): Not use temporary regularization except in final iterations?
template<typename Real>
void EquilibratedMehrotra
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int degree = n;
    const Grid& grid = problem.A.Grid();
    const int commRank = grid.Rank();
    Timer timer;

    // TODO(poulson): Move these into the control structure
    Real gammaPerm, deltaPerm, betaPerm, gammaTmp, deltaTmp, betaTmp;
    if( ctrl.system == NORMAL_KKT )
    {
        gammaPerm = deltaPerm = betaPerm = gammaTmp = deltaTmp = betaTmp = 0;
    }
    else
    {
        gammaPerm = ctrl.reg0Perm;
        deltaPerm = ctrl.reg1Perm;
        betaPerm = ctrl.reg2Perm;
        gammaTmp = ctrl.reg0Tmp;
        deltaTmp = ctrl.reg1Tmp;
        betaTmp = ctrl.reg2Tmp;
    }

    const Real bNrm2 = FrobeniusNorm( problem.b );
    const Real cNrm2 = FrobeniusNorm( problem.c );
    const Real twoNormEstA = TwoNormEstimate( problem.A, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstA + 1;
    if( ctrl.print )
    {
        const double imbalanceA = problem.A.Imbalance();
        if( commRank == 0 )
        {
            Output("|| A ||_2 estimate: ",twoNormEstA);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
            Output("Imbalance factor of A: ",imbalanceA);
        }
    }

    DistSparseLDLFactorization<Real> sparseLDLFact;
    // The initialization involves an augmented KKT system, and so we can
    // only reuse the factorization metadata if the this IPM is using the
    // augmented formulation
    if( commRank == 0 && ctrl.time )
        timer.Start();
    if( ctrl.system == AUGMENTED_KKT )
    {
        Initialize
        ( problem, solution, sparseLDLFact,
          ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift,
          ctrl.solveCtrl );
    }
    else
    {
        DistSparseLDLFactorization<Real> augmentedSparseLDLFact;
        Initialize
        ( problem, solution, augmentedSparseLDLFact,
          ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift,
          ctrl.solveCtrl );
    }
    if( commRank == 0 && ctrl.time )
        Output("Init: ",timer.Stop()," secs");

    DistMultiVec<Real> regTmp(grid);
    if( ctrl.system == FULL_KKT )
    {
        regTmp.Resize( m+2*n, 1 );
        for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
        {
            const Int i = regTmp.GlobalRow(iLoc);
            if( i < n )        regTmp.SetLocal( iLoc, 0,  gammaTmp*gammaTmp );
            else if( i < n+m ) regTmp.SetLocal( iLoc, 0, -deltaTmp*deltaTmp );
            else               regTmp.SetLocal( iLoc, 0, -betaTmp*betaTmp );
        }
    }
    else if( ctrl.system == AUGMENTED_KKT )
    {
        regTmp.Resize( n+m, 1 );
        for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
        {
            const Int i = regTmp.GlobalRow(iLoc);
            if( i < n ) regTmp.SetLocal( iLoc, 0,  gammaTmp*gammaTmp );
            else        regTmp.SetLocal( iLoc, 0, -deltaTmp*deltaTmp );
        }
    }
    else if( ctrl.system == NORMAL_KKT )
    {
        regTmp.Resize( m, 1 );
        Fill( regTmp, deltaTmp*deltaTmp );
    }
    regTmp *= origTwoNormEst;

    Real muOld = 0.1;
    Real relError = 1;

    DistGraphMultMeta metaOrig, meta;
    DistSparseMatrix<Real> J(grid), JOrig(grid);
    DistMultiVec<Real> d(grid), w(grid);
    DistMultiVec<Real> dInner(grid);

    DirectLPSolution<DistMultiVec<Real>> affineCorrection, correction;
    DirectLPResidual<DistMultiVec<Real>> residual, error;
    ForceSimpleAlignments( affineCorrection, grid );
    ForceSimpleAlignments( correction, grid );
    ForceSimpleAlignments( residual, grid );
    ForceSimpleAlignments( error, grid );

    DistMultiVec<Real> prod(grid);
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = pos_orth::NumOutside( solution.x );
        const Int zNumNonPos = pos_orth::NumOutside( solution.z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the barrier parameter
        // =============================
        Real mu = Dot(solution.x,solution.z) / degree;
        const Real compRatio =
          pos_orth::ComplementRatio( solution.x, solution.z );
        mu = compRatio > ctrl.balanceTol ? muOld : Min(mu,muOld);
        muOld = mu;

        pos_orth::NesterovTodd( solution.x, solution.z, w );
        const Real wMaxNorm = MaxNorm( w );

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        residual.primalEquality = problem.b;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(-1), residual.primalEquality );
        const Real rbNrm2 = FrobeniusNorm( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        Axpy( -deltaPerm*deltaPerm, solution.y, residual.primalEquality );
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        residual.dualEquality -= solution.z;
        const Real rcNrm2 = FrobeniusNorm( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        Axpy( gammaPerm*gammaPerm, solution.x, residual.dualEquality );
        // Now check the pieces
        // --------------------
        relError = Max(Max(objConv,rbConv),rcConv);
        if( ctrl.print )
        {
            const Real xNrm2 = FrobeniusNorm( solution.x );
            const Real yNrm2 = FrobeniusNorm( solution.y );
            const Real zNrm2 = FrobeniusNorm( solution.z );
            if( commRank == 0 )
                Output
                ("iter ",numIts,":\n",Indent(),
                 "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
                 "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
                 "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
                 "  ||  w  ||_max = ",wMaxNorm,"\n",Indent(),
                 "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
                 "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
                 "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
                 "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
                 "  primal = ",primObj,"\n",Indent(),
                 "  dual   = ",dualObj,"\n",Indent(),
                 "  |primal - dual| / (1 + |primal|) = ",objConv);
        }
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := x o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.x, residual.dualConic );

        if( ctrl.system == FULL_KKT || ctrl.system == AUGMENTED_KKT )
        {
            // Assemble the KKT system
            // -----------------------
            if( ctrl.system == FULL_KKT )
            {
                KKT
                ( problem.A, gammaPerm, deltaPerm, betaPerm,
                  solution.x, solution.z, JOrig, false );
                KKTRHS
                ( residual.dualEquality, residual.primalEquality,
                  residual.dualConic, solution.z, d );
            }
            else
            {
                AugmentedKKT
                ( problem.A, gammaPerm, deltaPerm, solution.x, solution.z,
                  JOrig, false );
                AugmentedKKTRHS
                ( solution.x, residual.dualEquality, residual.primalEquality,
                  residual.dualConic, d );
            }
            J = JOrig;
            UpdateDiagonal( J, Real(1), regTmp );
            if( numIts == 0 )
            {
                metaOrig = JOrig.InitializeMultMeta();
                meta = J.InitializeMultMeta();
                if( ctrl.print )
                {
                    const double imbalanceJ = J.Imbalance();
                    if( commRank == 0 )
                        Output("Imbalance factor of J: ",imbalanceJ);
                }
            }
            else
            {
                JOrig.LockedDistGraph().multMeta = metaOrig;
                J.LockedDistGraph().multMeta = meta;
            }

            // Solve for the direction
            // -----------------------
            try
            {
                if( commRank == 0 && ctrl.time )
                    timer.Start();
                if( wMaxNorm >= ctrl.ruizEquilTol )
                {
                    if( ctrl.print && commRank == 0 )
                        Output("Running SymmetricRuizEquil");
                    SymmetricRuizEquil
                    ( J, dInner, ctrl.ruizMaxIter, ctrl.print );
                }
                else if( wMaxNorm >= ctrl.diagEquilTol )
                {
                    if( ctrl.print && commRank == 0 )
                        Output("Running SymmetricDiagonalEquil");
                    SymmetricDiagonalEquil( J, dInner, ctrl.print );
                }
                else
                    Ones( dInner, J.Height(), 1 );
                if( commRank == 0 && ctrl.time )
                    Output("Equilibration: ",timer.Stop()," secs");

                if( numIts == 0 &&
                    (ctrl.system != AUGMENTED_KKT ||
                     (ctrl.primalInit && ctrl.dualInit)) )
                {
                    if( commRank == 0 && ctrl.time )
                        timer.Start();
                    const bool hermitian = true;
                    const BisectCtrl bisectCtrl;
                    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );
                    if( commRank == 0 && ctrl.time )
                        Output("Analysis: ",timer.Stop()," secs");
                }
                else
                    sparseLDLFact.ChangeNonzeroValues( J );

                if( commRank == 0 && ctrl.time )
                    timer.Start();
                sparseLDLFact.Factor( LDL_2D );
                if( commRank == 0 && ctrl.time )
                    Output("LDL: ",timer.Stop()," secs");

                if( commRank == 0 && ctrl.time )
                    timer.Start();
                if( ctrl.resolveReg )
                    reg_ldl::SolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d, ctrl.solveCtrl );
                else
                    reg_ldl::RegularizedSolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d,
                      ctrl.solveCtrl.relTol,
                      ctrl.solveCtrl.maxRefineIts,
                      ctrl.solveCtrl.progress );
                if( commRank == 0 && ctrl.time )
                    Output("Affine: ",timer.Stop()," secs");
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }

            if( ctrl.system == FULL_KKT )
                ExpandSolution
                ( m, n, d,
                  affineCorrection.x, affineCorrection.y, affineCorrection.z );
            else
                ExpandAugmentedSolution
                ( solution.x, solution.z, residual.dualConic, d,
                  affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }
        else // ctrl.system == NORMAL_KKT
        {
            // Assemble the KKT system
            // -----------------------
            // TODO(poulson): Apply updates on top of explicit zeros
            NormalKKT
            ( problem.A, gammaPerm, deltaPerm, solution.x, solution.z,
              J, false );
            NormalKKTRHS
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, affineCorrection.y );
            if( numIts == 0 )
            {
                if( ctrl.print )
                {
                    const double imbalanceJ = J.Imbalance();
                    if( commRank == 0 )
                        Output("Imbalance factor of J: ",imbalanceJ);
                }
                meta = J.InitializeMultMeta();
            }
            else
            {
                J.LockedDistGraph().multMeta = meta;
            }

            // Solve for the direction
            // -----------------------
            try
            {
                if( numIts == 0 )
                {
                    if( commRank == 0 && ctrl.time )
                        timer.Start();
                    const bool hermitian = true;
                    const BisectCtrl bisectCtrl;
                    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );
                    if( commRank == 0 && ctrl.time )
                        Output("Analysis: ",timer.Stop()," secs");
                }
                else
                {
                    sparseLDLFact.ChangeNonzeroValues( J );
                }

                if( commRank == 0 && ctrl.time )
                    timer.Start();
                sparseLDLFact.Factor( LDL_2D );
                if( commRank == 0 && ctrl.time )
                    Output("LDL: ",timer.Stop()," secs");

                if( commRank == 0 && ctrl.time )
                    timer.Start();
                reg_ldl::RegularizedSolveAfter
                ( J, regTmp, sparseLDLFact, affineCorrection.y,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress,
                  ctrl.solveCtrl.time );
                if( commRank == 0 && ctrl.time )
                    Output("Affine: ",timer.Stop()," secs");
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandNormalSolution
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }

        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Gemv
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            Axpy
            ( -deltaPerm*deltaPerm, affineCorrection.y, error.primalEquality );
            Real dxErrorNrm2 = FrobeniusNorm( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Gemv
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Axpy( gammaPerm*gammaPerm, affineCorrection.x, error.dualEquality );
            error.dualEquality -= affineCorrection.z;
            Real dyErrorNrm2 = FrobeniusNorm( error.dualEquality );

            Real rmuNrm2 = FrobeniusNorm( residual.dualConic );
            error.dualConic = residual.dualConic;
            prod = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, solution.x, prod );
            error.dualConic += prod;
            prod = affineCorrection.x;
            DiagonalScale( LEFT, NORMAL, solution.z, prod );
            error.dualConic += prod;
            Real dzErrorNrm2 = FrobeniusNorm( error.dualConic );

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri =
          pos_orth::MaxStep( solution.x, affineCorrection.x, Real(1) );
        Real alphaAffDual =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: correction.z and correction.x are used as temporaries
        correction.x = solution.x;
        correction.z = solution.z;
        Axpy( alphaAffPri,  affineCorrection.x, correction.x );
        Axpy( alphaAffDual, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.x,correction.z) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma =
          ctrl.centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        residual.primalEquality *= 1-sigma;
        residual.dualEquality *= 1-sigma;
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dxAff o dzAff
            // ---------------------
            // NOTE: correction.z is used as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.x, correction.z );
            residual.dualConic += correction.z;
        }

        if( ctrl.system == FULL_KKT )
        {
            KKTRHS
            ( residual.dualEquality, residual.primalEquality,
              residual.dualConic, solution.z, d );
            try
            {
                if( commRank == 0 && ctrl.time )
                    timer.Start();
                if( ctrl.resolveReg )
                    reg_ldl::SolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d, ctrl.solveCtrl );
                else
                    reg_ldl::RegularizedSolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d,
                      ctrl.solveCtrl.relTol,
                      ctrl.solveCtrl.maxRefineIts,
                      ctrl.solveCtrl.progress );
                if( commRank == 0 && ctrl.time )
                    Output("Corrector: ",timer.Stop()," secs");
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandSolution( m, n, d, correction.x, correction.y, correction.z );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            AugmentedKKTRHS
            ( solution.x, residual.dualEquality, residual.primalEquality,
              residual.dualConic, d );
            try
            {
                if( commRank == 0 && ctrl.time )
                    timer.Start();
                if( ctrl.resolveReg )
                    reg_ldl::SolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d, ctrl.solveCtrl );
                else
                    reg_ldl::RegularizedSolveAfter
                    ( JOrig, regTmp, dInner, sparseLDLFact, d,
                      ctrl.solveCtrl.relTol,
                      ctrl.solveCtrl.maxRefineIts,
                      ctrl.solveCtrl.progress );
                if( commRank == 0 && ctrl.time )
                    Output("Corrector: ",timer.Stop()," secs");
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, d,
              correction.x, correction.y, correction.z );
        }
        else
        {
            NormalKKTRHS
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, correction.y );
            try
            {
                if( commRank == 0 && ctrl.time )
                    timer.Start();
                reg_ldl::RegularizedSolveAfter
                ( J, regTmp, sparseLDLFact, correction.y,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress,
                  ctrl.solveCtrl.time );
                if( commRank == 0 && ctrl.time )
                    Output("Corrector: ",timer.Stop()," secs");
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandNormalSolution
            ( problem.A, gammaPerm, solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              correction.x, correction.y, correction.z );
        }
        // TODO(poulson): Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri =
          pos_orth::MaxStep( solution.x, correction.x, 1/ctrl.maxStepRatio );
        Real alphaDual =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaDual, correction.y, solution.y );
        Axpy( alphaDual, correction.z, solution.z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
    }
    SetIndent( indent );
}

template<typename Real>
void Mehrotra
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.outerEquil )
    {
        DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>
          equilibratedProblem;
        DirectLPSolution<DistMultiVec<Real>> equilibratedSolution;
        DistSparseDirectLPEquilibration<Real> equilibration;
        Equilibrate
        ( problem, solution,
          equilibratedProblem, equilibratedSolution, equilibration, ctrl );
        EquilibratedMehrotra( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedMehrotra( problem, solution, ctrl );
    }
    if( ctrl.print )
    {
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        const Real xNrm2 = FrobeniusNorm( solution.x );
        const Real yNrm2 = FrobeniusNorm( solution.y );
        const Real zNrm2 = FrobeniusNorm( solution.z );
        OutputFromRoot
        (problem.A.Grid().Comm(),
         "Exiting with:\n",Indent(),
         "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
         "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
         "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
         "  primal = ",primObj,"\n",Indent(),
         "  dual   = ",dualObj,"\n",Indent(),
         "  |primal - dual| / (1 + |primal|) = ",objConv);
    }
}

// This interface is now deprecated.
template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = A.Grid();

    DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>> problem;
    DirectLPSolution<DistMultiVec<Real>> solution;
    ForceSimpleAlignments( problem, grid );
    ForceSimpleAlignments( solution, grid );

    problem.c = c;
    problem.A = A;
    problem.b = b;
    solution.x = x;
    solution.y = y;
    solution.z = z;
    Mehrotra( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
          DirectLPSolution<Matrix<Real>>& solution, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
          DirectLPSolution<DistMatrix<Real>>& solution, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
          DirectLPSolution<Matrix<Real>>& solution, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
          DirectLPSolution<DistMultiVec<Real>>& solution, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace direct
} // namespace lp
} // namespace El
