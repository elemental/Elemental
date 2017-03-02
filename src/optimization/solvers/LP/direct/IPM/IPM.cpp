/*
   Copyright (c) 2009-2017, Jack Poulson
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
//   s.t. A^T y + G^T z + c = 0, z >= 0.
//
// We make use of the regularized Lagrangian
//
//   L(x;y,z) = c^T x + y^T (A x - b) - z^T x
//              + (1/2) gamma_x || x - x_0 ||_2^2
//              - (1/2) gamma_y || y - y_0 ||_2^2
//              - (1/2) gamma_z || z - z_0 ||_2^2
//              + mu Phi(z).
//
// where we note that the two-norm regularization is positive for the primal
// variable x and *negative* for the dual variables y and z. NOTE: While z is
// regularized in the affine solver, it is not yet implemented for direct
// solvers.
//
// The subsequent first-order optimality conditions for x, y, and z become
//
//   Nabla_x L = c + A^T y - z + gamma_x (x - x_0) = 0,
//   Nabla_y L = A x - b - gamma_y (y - y_0) = 0,
//
// These can be arranged into the symmetric quasi-definite form
//
//   | gamma_x I,      A^T,  | | x | = | -c + gamma_x x_0 |.
//   |      A,    -gamma_y I | | y |   |  b - gamma_y y_0 |
//
// The regularization on z implies the optimality condition
//
//   Nabla_z L = -x - gamma_z (z - z_0) - mu inv(z) = 0,
//
// which simplifies to
//
//   - z o (x + gamma_z (z - z_0)) = mu e.
//

// TODO(poulson): Experiment with using the lagged factorization as a
// preconditioner.

// TODO(poulson): Make these variables consistent with the affine
// implementation.
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
  const IPMCtrl<Real>& ctrl )
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

    if( ctrl.print )
    {
        Output("xScale=",equilibration.xScale);
        Output("zScale=",equilibration.zScale);
        Output("|| rowScale ||_2 = ",FrobeniusNorm(equilibration.rowScale));
        Output("|| colScale ||_2 = ",FrobeniusNorm(equilibration.colScale));
    }
}

template<typename Real>
void Equilibrate
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const DirectLPSolution<DistMatrix<Real>>& solution,
        DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& equilibratedProblem,
        DirectLPSolution<DistMatrix<Real>>& equilibratedSolution,
        DistDenseDirectLPEquilibration<Real>& equilibration,
  const IPMCtrl<Real>& ctrl )
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

    if( ctrl.print && problem.A.Grid().Rank() == 0 )
    {
        Output("xScale=",equilibration.xScale);
        Output("zScale=",equilibration.zScale);
        Output("|| rowScale ||_2 = ",FrobeniusNorm(equilibration.rowScale));
        Output("|| colScale ||_2 = ",FrobeniusNorm(equilibration.colScale));
    }
}

template<typename Real>
void Equilibrate
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const DirectLPSolution<Matrix<Real>>& solution,
        DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& equilibratedProblem,
        DirectLPSolution<Matrix<Real>>& equilibratedSolution,
        SparseDirectLPEquilibration<Real>& equilibration,
  const IPMCtrl<Real>& ctrl )
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

    if( ctrl.print )
    {
        Output("xScale=",equilibration.xScale);
        Output("zScale=",equilibration.zScale);
        Output("|| rowScale ||_2 = ",FrobeniusNorm(equilibration.rowScale));
        Output("|| colScale ||_2 = ",FrobeniusNorm(equilibration.colScale));
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
  const IPMCtrl<Real>& ctrl )
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

    if( ctrl.print && problem.A.Grid().Rank() == 0 )
    {
        Output("xScale=",equilibration.xScale);
        Output("zScale=",equilibration.zScale);
        Output("|| rowScale ||_2 = ",FrobeniusNorm(equilibration.rowScale));
        Output("|| colScale ||_2 = ",FrobeniusNorm(equilibration.colScale));
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
    Real relCompGap;
    Real relObjGap;

    DirectLPResidual<Matrix<Real>> residual;
    Real primalEqualityNorm;
    Real dualEqualityNorm;
    Real dualConicNorm;
    Real relativePrimalEqualityNorm;
    Real relativeDualEqualityNorm;
    Real infeasError;

    Int numIts;
    Real dimacsError;
    Real dimacsErrorOld;

    bool metTolerances;

    void Initialize
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const IPMCtrl<Real>& ctrl );

    void Update
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const DirectLPSolution<Matrix<Real>>& solution,
      const DirectRegularization<Real>& smallReg,
      const IPMCtrl<Real>& ctrl );

    void PrintResiduals
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const DirectLPSolution<Matrix<Real>>& solution,
      const DirectLPSolution<Matrix<Real>>& correction,
      const DirectRegularization<Real>& smallReg ) const;
};

template<typename Real>
void DirectState<Real,Matrix<Real>,Matrix<Real>>::Initialize
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    bNorm = FrobeniusNorm( problem.b );
    cNorm = FrobeniusNorm( problem.c );
    barrierOld = 0.1;
    infeasError = 1;
    relCompGap = 1;
    relObjGap = 1;
    dimacsError = 1;
    dimacsErrorOld = 1;
    metTolerances = false;
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
  const DirectRegularization<Real>& smallReg,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int degree = problem.A.Width();
    const Real epsilon = limits::Epsilon<Real>();

    // Compute the new barrier parameter
    // ---------------------------------
    const Real dualProd = Dot(solution.x,solution.z);
    barrier = dualProd / degree;
    const Real compRatio = pos_orth::ComplementRatio( solution.x, solution.z );
    barrier = compRatio > ctrl.maxComplementRatio ? barrierOld
      : Min(barrier,barrierOld);
    barrierOld = barrier;

    // Compute the objectives and relative duality gap
    primalObjective = Dot(problem.c,solution.x);
    dualObjective = -Dot(problem.b,solution.y);
    relObjGap = RelativeObjectiveGap( primalObjective, dualObjective );
    relCompGap =
      RelativeComplementarityGap( primalObjective, dualObjective, dualProd );
    const Real maxRelGap = Max( relCompGap, relObjGap );

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
    dimacsErrorOld = dimacsError;
    infeasError = Max(relativePrimalEqualityNorm,relativeDualEqualityNorm);
    dimacsError = Max(infeasError,maxRelGap);

    metTolerances =
      infeasError <= Pow(epsilon,ctrl.infeasibilityTolLogEps) &&
      relCompGap <= Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
      relObjGap <= Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);

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
         "  || primalInfeas ||_2 / (1 + || b ||_2) = ",
         relativePrimalEqualityNorm,"\n",Indent(),
         "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",
         relativeDualEqualityNorm,"\n",Indent(),
         "  primal = ",primalObjective,"\n",Indent(),
         "  dual   = ",dualObjective,"\n",Indent(),
         "  relative gap = ",maxRelGap,"\n",Indent(),
         "  DIMACS: ",dimacsError);
    }
}

template<typename Real>
void DirectState<Real,Matrix<Real>,Matrix<Real>>::PrintResiduals
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const DirectLPSolution<Matrix<Real>>& solution,
  const DirectLPSolution<Matrix<Real>>& correction,
  const DirectRegularization<Real>& smallReg ) const
{
    EL_DEBUG_CSE
    DirectLPResidual<Matrix<Real>> error;
    Matrix<Real> prod;

    error.primalEquality = residual.primalEquality;
    Gemv
    ( NORMAL, Real(1), problem.A, correction.x,
      Real(1), error.primalEquality );
    Axpy
    ( -smallReg.primalEquality, correction.y, error.primalEquality );
    Real dxErrorNrm2 = FrobeniusNorm( error.primalEquality );

    error.dualEquality = residual.dualEquality;
    Gemv
    ( TRANSPOSE, Real(1), problem.A, correction.y,
      Real(1), error.dualEquality );
    Axpy( smallReg.dualEquality, correction.x, error.dualEquality );
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
    ("|| dxError ||_2 / (1 + || primalInfeas ||_2) = ",
     dxErrorNrm2/(1+primalEqualityNorm),"\n",Indent(),
     "|| dyError ||_2 / (1 + || dualInfeas   ||_2) = ",
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
      const DirectRegularization<Real>& smallReg,
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
            ( problem.A,
              Sqrt(smallReg.dualEquality),
              Sqrt(smallReg.primalEquality),
              solution.x, solution.z, kktSystem );
        }
        Factor();
    }

    void SolveSystem
    ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
      const DirectRegularization<Real>& smallReg,
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
            ( problem.A, Sqrt(smallReg.dualEquality), solution.x, solution.z,
              residual.dualEquality,
              residual.primalEquality,
              residual.dualConic,
              correction.y );
            Solve(correction.y);
            ExpandNormalSolution
            ( problem, Sqrt(smallReg.dualEquality), solution, residual,
              correction );
        }
    }
};

template<typename Real>
void EquilibratedIPM
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl,
  bool outputRoot )
{
    EL_DEBUG_CSE
    const Int n = problem.A.Width();
    const Int degree = n;

    // TODO(poulson): Implement nonzero regularization
    DirectRegularization<Real> smallReg;
    smallReg.primalEquality = 0;
    smallReg.dualEquality = 0;

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

        state.Update( problem, solution, smallReg, ctrl );

        // Check for convergence
        // =====================
        if( state.metTolerances )
        {
            if( state.dimacsError >=
                ctrl.minDimacsDecreaseRatio*state.dimacsErrorOld )
            {
                // We have met the tolerances and progress in the last iteration
                // was not significant.
                break;
            }
            else if( state.numIts == ctrl.maxIts )
            {
                // We have hit the iteration limit but can declare success.
                break;
            }
        }
        else if( state.numIts == ctrl.maxIts )
        {
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving tolerances");
        }

        // Set up and factor the KKT matrix
        // ================================
        solver.InitializeSystem( problem, smallReg, solution, ctrl.system );

        // Compute the affine search direction
        // ===================================
        solver.SolveSystem
        ( problem, smallReg, state.residual, solution, affineCorrection,
          ctrl.system );
        /*
        if( ctrl.checkResiduals && ctrl.print )
        {
            state.PrintResiduals
            ( problem, solution, affineCorrection, smallReg );
        }
        */

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
        if( ctrl.compositeNewton )
        {
            // r_mu += dxAff o dzAff
            // ---------------------
            // NOTE: We are using correction.z as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.x, correction.z );
            state.residual.dualConic += correction.z;
        }
        solver.SolveSystem
        ( problem, smallReg, state.residual, solution, correction,
          ctrl.system );
        /*
        if( ctrl.checkResiduals && ctrl.print )
        {
            state.PrintResiduals( problem, solution, correction, smallReg );
        }
        */

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
    if( !state.metTolerances )
        RuntimeError("Unable to achieve tolerances");

    SetIndent( indent );
}

template<typename Real>
void IPM
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
        EquilibratedIPM
        ( equilibratedProblem, equilibratedSolution, ctrl, outputRoot );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedIPM( problem, solution, ctrl, outputRoot );
    }
    if( ctrl.print && outputRoot )
    {
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real dualProd = Dot(solution.x,solution.z);
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relObjGap, relCompGap );

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
         "  relativeDualityGap = ",maxRelGap);
    }
}

// This interface is now deprecated.
template<typename Real>
void IPM
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const IPMCtrl<Real>& ctrl )
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
    IPM( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

template<typename Real>
void EquilibratedIPM
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int degree = n;
    const Grid& grid = problem.A.Grid();
    const int commRank = grid.Rank();
    const Real epsilon = limits::Epsilon<Real>();

    const Real bNrm2 = FrobeniusNorm( problem.b );
    const Real cNrm2 = FrobeniusNorm( problem.c );
    const Real ANrm1 = OneNorm( problem.A );
    if( ctrl.print )
    {
        if( commRank == 0 )
        {
            Output("|| A ||_1 = ",ANrm1);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
        }
    }

    const Real xRegSmall0 = ANrm1*Pow(epsilon,ctrl.xRegSmallLogEps);
    const Real yRegSmall0 = ANrm1*Pow(epsilon,ctrl.yRegSmallLogEps);
    Real xRegSmall = xRegSmall0;
    Real yRegSmall = yRegSmall0;

    Initialize
    ( problem, solution,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift );

    Real muOld = 0.1;
    Real infeasError = 1;
    Real dimacsError = 1, dimacsErrorOld = 1;
    DistMatrix<Real> J(grid), d(grid);
    DistMatrix<Real> dSub(grid);
    DistPermutation p(grid);
    auto attemptToFactor = [&]()
      {
        try { LDL( J, dSub, p, false ); }
        catch(...) { return false; }
        return true;
      };
    auto attemptToSolve = [&]( DistMatrix<Real>& rhs )
      {
        try { ldl::SolveAfter( J, dSub, p, rhs, false ); }
        catch(...) { return false; }
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
    Int numIts=0;
    for( ; numIts<=ctrl.maxIts; ++numIts, dimacsErrorOld=dimacsError )
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
        const Real dualProd = Dot(solution.x,solution.z);
        Real mu = dualProd / degree;
        const Real compRatio =
          pos_orth::ComplementRatio( solution.x, solution.z );
        mu = compRatio > ctrl.maxComplementRatio ? muOld : Min(mu,muOld);
        muOld = mu;

        // Check for convergence
        // =====================

        // Check the duality gaps
        // ----------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relObjGap, relCompGap );

        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(-1), residual.primalEquality );
        const Real rbNrm2 = FrobeniusNorm( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);

        // || A^T y - z + c ||_2 / (1 + || c ||_2)
        // ---------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        residual.dualEquality -= solution.z;
        const Real rcNrm2 = FrobeniusNorm( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);

        // Now check the pieces
        // --------------------
        infeasError = Max(rbConv,rcConv);
        dimacsError = Max(maxRelGap,infeasError);
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
                 "  || primalInfeas ||_2 / (1 + || b ||_2) = ",
                 rbConv,"\n",Indent(),
                 "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",
                 rcConv,"\n",Indent(),
                 "  primal = ",primObj,"\n",Indent(),
                 "  dual   = ",dualObj,"\n",Indent(),
                 "  relative duality gap = ",maxRelGap);
        }

        const bool metTolerances =
          infeasError <= Pow(epsilon,ctrl.infeasibilityTolLogEps) &&
          relCompGap <= Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
          relObjGap <= Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);
        if( metTolerances )
        {
            if( dimacsError >= ctrl.minDimacsDecreaseRatio*dimacsErrorOld )
            {
                // We have met the tolerances and progress in the last iteration
                // was not significant.
                break;
            }
            else if( numIts == ctrl.maxIts )
            {
                // We have hit the iteration limit but can declare success.
                break;
            }
        }
        else if( numIts == ctrl.maxIts )
        {
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving tolerances");
        }

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
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
            if( !attemptToSolve(d) )
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
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
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
            if( !attemptToSolve(d) )
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, d,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }
        else if( ctrl.system == NORMAL_KKT )
        {
            // Construct the KKT system
            // ------------------------
            NormalKKT
            ( problem.A, Sqrt(xRegSmall), Sqrt(yRegSmall),
              solution.x, solution.z, J );
            NormalKKTRHS
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, affineCorrection.y );

            // Solve for the direction
            // -----------------------
            if( !attemptToFactor() )
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
            if( !attemptToSolve(affineCorrection.y) )
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
            ExpandNormalSolution
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }

        /*
        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Gemv
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            Axpy( -yRegSmall, affineCorrection.y, error.primalEquality );
            Real dxErrorNrm2 = FrobeniusNorm( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Gemv
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Axpy( xRegSmall, affineCorrection.x, error.dualEquality );
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
                ("|| dxError ||_2 / (1 + || primalInfeas ||_2) = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || dualInfeas   ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rmuNrm2));
        }
        */

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
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.compositeNewton )
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
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
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
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, d,
              correction.x, correction.y, correction.z );
        }
        else if( ctrl.system == NORMAL_KKT )
        {
            // Construct the new KKT RHS
            // -------------------------
            NormalKKTRHS
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, correction.y );

            // Solve for the direction
            // -----------------------
            if( !attemptToSolve(correction.y) )
            {
                // TODO(poulson): Increase regularization and continue.
                break;
            }
            ExpandNormalSolution
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
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
            if( metTolerances )
            {
                break;
            }
            else
            {
                RuntimeError("Could not achieve tolerances");
            }
        }
    }
    SetIndent( indent );
}

template<typename Real>
void IPM
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
        EquilibratedIPM( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        // Avoid creating unnecessary copies where we can.
        if( SimpleAlignments(problem) && SimpleAlignments(solution) )
        {
            EquilibratedIPM( problem, solution, ctrl );
        }
        else if( SimpleAlignments(problem) )
        {
            DirectLPSolution<DistMatrix<Real>> alignedSolution;
            ForceSimpleAlignments( alignedSolution, grid );
            alignedSolution = solution;
            EquilibratedIPM( problem, alignedSolution, ctrl );
            solution = alignedSolution;
        }
        else if( SimpleAlignments(solution) )
        {
            DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>> alignedProblem;
            ForceSimpleAlignments( alignedProblem, grid );
            CopyOrViewHelper( problem.c, alignedProblem.c );
            CopyOrViewHelper( problem.A, alignedProblem.A );
            CopyOrViewHelper( problem.b, alignedProblem.b );
            EquilibratedIPM( alignedProblem, solution, ctrl );
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
            EquilibratedIPM( alignedProblem, alignedSolution, ctrl );
            solution = alignedSolution;
        }
    }
    if( ctrl.print )
    {
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real dualProd = Dot(solution.x,solution.z);
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relObjGap, relCompGap );

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
         "  relative duality gap = ",maxRelGap);
    }
}

// This interface is now deprecated.
template<typename Real>
void IPM
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const IPMCtrl<Real>& ctrl )
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
    IPM( problem, solution, ctrl );
    Copy( solution.x, x );
    Copy( solution.y, y );
    Copy( solution.z, z );
}

template<typename Real>
void EquilibratedIPM
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int degree = n;
    const Real epsilon = limits::Epsilon<Real>();

    const Real bNrm2 = FrobeniusNorm( problem.b );
    const Real cNrm2 = FrobeniusNorm( problem.c );
    const Real twoNormEstA =
      TwoNormEstimate( problem.A, ctrl.twoNormKrylovBasisSize );
    const Real origTwoNormEst = twoNormEstA + 1;
    if( ctrl.print )
    {
        Output("|| A ||_2 estimate: ",twoNormEstA);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
    }

    const Real xRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.xRegSmallLogEps);
    const Real yRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.yRegSmallLogEps);
    const Real zRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.zRegSmallLogEps);
    const Real xRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.xRegLargeLogEps);
    const Real yRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.yRegLargeLogEps);
    const Real zRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.zRegLargeLogEps);
    Real xRegLarge = xRegLarge0;
    Real yRegLarge = yRegLarge0;
    Real zRegLarge = zRegLarge0;
    Real xRegSmall = xRegSmall0;
    Real yRegSmall = yRegSmall0;
    Real zRegSmall = zRegSmall0;

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

    Matrix<Real> regLarge;
    if( ctrl.system == FULL_KKT )
    {
        regLarge.Resize( m+2*n, 1 );
        for( Int i=0; i<m+2*n; ++i )
        {
            if( i < n )        regLarge(i) =  xRegLarge;
            else if( i < n+m ) regLarge(i) = -yRegLarge;
            else               regLarge(i) = -zRegLarge;
        }
    }
    else if( ctrl.system == AUGMENTED_KKT )
    {
        regLarge.Resize( n+m, 1 );
        for( Int i=0; i<n+m; ++i )
        {
            if( i < n ) regLarge(i) =  xRegLarge;
            else        regLarge(i) = -yRegLarge;
        }
    }
    else if( ctrl.system == NORMAL_KKT )
    {
        regLarge.Resize( m, 1 );
        Fill( regLarge, yRegLarge );
    }
    regLarge *= origTwoNormEst;

    Real muOld = 0.1;
    Real infeasError = 1;
    Real dimacsError = 1, dimacsErrorOld = 1;
    SparseMatrix<Real> J, JOrig;
    Matrix<Real> d, w;
    Matrix<Real> dInner;

    DirectLPSolution<Matrix<Real>> affineCorrection, correction;
    DirectLPResidual<Matrix<Real>> residual, error;

    Matrix<Real> prod;
    const Int indent = PushIndent();
    Int numIts=0;
    for( ; numIts<=ctrl.maxIts; ++numIts, dimacsErrorOld=dimacsError )
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

        // Check the duality gap
        // ---------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real dualProd = Dot(solution.x,solution.z);
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relObjGap, relCompGap );

        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(-1), residual.primalEquality );
        const Real rbNrm2 = FrobeniusNorm( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);

        // || A^T y - z + c ||_2 / (1 + || c ||_2)
        // ---------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        residual.dualEquality -= solution.z;
        const Real rcNrm2 = FrobeniusNorm( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);

        // Now check the pieces
        // --------------------
        infeasError = Max(rbConv,rcConv);
        dimacsError = Max(maxRelGap,infeasError);

        // Compute the scaling point
        // =========================
        pos_orth::NesterovTodd( solution.x, solution.z, w );
        const Real wMaxNorm = MaxNorm( w );

        // Compute the barrier parameter
        // =============================
        Real mu = Dot(solution.x,solution.z) / degree;
        const Real compRatio =
          pos_orth::ComplementRatio( solution.x, solution.z );
        mu = compRatio > ctrl.maxComplementRatio ? muOld : Min(mu,muOld);
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
             "  || primalInfeas ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
             "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
             "  mu        = ",mu,"\n",Indent(),
             "  primal    = ",primObj,"\n",Indent(),
             "  dual      = ",dualObj,"\n",Indent(),
             "  relative duality gap = ",maxRelGap);
        }

        const bool metTolerances =
          infeasError <= Pow(epsilon,ctrl.infeasibilityTolLogEps) &&
          relCompGap <= Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
          relObjGap <= Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);
        if( metTolerances )
        {
            if( dimacsError >= ctrl.minDimacsDecreaseRatio*dimacsErrorOld )
            {
                // We have met the tolerances and progress in the last iteration
                // was not significant.
                break;
            }
            else if( numIts == ctrl.maxIts )
            {
                // We have hit the iteration limit but can declare success.
                break;
            }
        }
        else if( numIts == ctrl.maxIts )
        {
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving tolerances");
        }

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
                ( problem.A,
                  Sqrt(xRegSmall),
                  Sqrt(yRegSmall),
                  Sqrt(zRegSmall),
                  solution.x, solution.z, JOrig, false );
                KKTRHS
                ( residual.dualEquality, residual.primalEquality,
                  residual.dualConic, solution.z, d );
            }
            else
            {
                AugmentedKKT
                ( problem.A, Sqrt(xRegSmall), Sqrt(yRegSmall),
                  solution.x, solution.z, JOrig, false );
                AugmentedKKTRHS
                ( solution.x, residual.dualEquality, residual.primalEquality,
                  residual.dualConic, d );
            }
            J = JOrig;
            UpdateDiagonal( J, Real(1), regLarge );

            // Solve for the direction
            // -----------------------
            /*
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
            */
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
            RegSolveInfo<Real> solveInfo;
            if( ctrl.twoStage )
            {
                solveInfo = reg_ldl::SolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d, ctrl.solveCtrl );
            }
            if( !solveInfo.metRequestedTol )
            {
                solveInfo = reg_ldl::RegularizedSolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress );
                if( !solveInfo.metRequestedTol )
                {
                    if( metTolerances )
                    {
                        break;
                    }
                    else
                    {
                        // TODO(poulson): Increase regularization and continue.
                        RuntimeError("Could not achieve tolerances");
                    }
                }
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
            ( problem.A, Sqrt(xRegSmall), Sqrt(yRegSmall),
              solution.x, solution.z, J, false );
            NormalKKTRHS
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, affineCorrection.y );

            // Solve for the direction
            // -----------------------
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

            // NOTE: regLarge should be all zeros; replace with unregularized
            auto solveInfo = reg_ldl::RegularizedSolveAfter
            ( J, regLarge, sparseLDLFact, affineCorrection.y,
              ctrl.solveCtrl.relTol,
              ctrl.solveCtrl.maxRefineIts,
              ctrl.solveCtrl.progress,
              ctrl.solveCtrl.time );
            if( !solveInfo.metRequestedTol )
            {
                if( metTolerances )
                {
                    break;
                }
                else
                {
                    // TODO(poulson): Increase regularization and continue.
                    RuntimeError("Could not achieve tolerances");
                }
            }
            ExpandNormalSolution
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }

        /*
        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Gemv
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            Axpy
            ( -yRegSmall, affineCorrection.y, error.primalEquality );
            Real dxErrorNrm2 = FrobeniusNorm( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Gemv
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Axpy( xRegSmall, affineCorrection.x, error.dualEquality );
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
            ("|| dxError ||_2 / (1 + || primalInfeas ||_2) = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || dualInfeas   ||_2) = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+rmuNrm2));
        }
        */

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
        Shift( residual.dualConic, -sigma*mu );
        // TODO(poulson): Gondzio's corrections
        if( ctrl.compositeNewton )
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
            RegSolveInfo<Real> solveInfo;
            if( ctrl.twoStage )
            {
                solveInfo = reg_ldl::SolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d, ctrl.solveCtrl );
            }
            if( !solveInfo.metRequestedTol )
            {
                solveInfo = reg_ldl::RegularizedSolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress );
                if( !solveInfo.metRequestedTol )
                {
                    if( metTolerances )
                    {
                        break;
                    }
                    else
                    {
                        // TODO(poulson): Increase regularization and continue.
                        RuntimeError("Could not achieve tolerances");
                    }
                }
            }
            ExpandSolution( m, n, d, correction.x, correction.y, correction.z );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            AugmentedKKTRHS
            ( solution.x, residual.dualEquality, residual.primalEquality,
              residual.dualConic, d );
            RegSolveInfo<Real> solveInfo;
            if( ctrl.twoStage )
            {
                solveInfo = reg_ldl::SolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d, ctrl.solveCtrl );
            }
            if( !solveInfo.metRequestedTol )
            {
                solveInfo = reg_ldl::RegularizedSolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress );
                if( !solveInfo.metRequestedTol )
                {
                    if( metTolerances )
                    {
                        break;
                    }
                    else
                    {
                        // TODO(poulson): Increase regularization and continue.
                        RuntimeError("Could not achieve tolerances");
                    }
                }
            }
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, d,
              correction.x, correction.y, correction.z );
        }
        else
        {
            NormalKKTRHS
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, correction.y );
            // NOTE: regLarge should be all zeros; replace with unregularized
            auto solveInfo = reg_ldl::RegularizedSolveAfter
            ( J, regLarge, sparseLDLFact, correction.y,
              ctrl.solveCtrl.relTol,
              ctrl.solveCtrl.maxRefineIts,
              ctrl.solveCtrl.progress,
              ctrl.solveCtrl.time );
            if( !solveInfo.metRequestedTol )
            {
                if( metTolerances )
                {
                    break;
                }
                else
                {
                    // TODO(poulson): Increase regularization and continue.
                    RuntimeError("Could not achieve tolerances");
                }
            }
            ExpandNormalSolution
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
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
            if( metTolerances )
            {
                break;
            }
            else
            {
                // TODO(poulson): Increase regularization and continue.
                RuntimeError("Could not achieve tolerances");
            }
        }
    }
    SetIndent( indent );
}

template<typename Real>
void IPM
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
        EquilibratedIPM( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedIPM( problem, solution, ctrl );
    }
    if( ctrl.print )
    {
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real dualProd = Dot(solution.x,solution.z);
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relObjGap, relCompGap );

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
         "  relative duality gap = ",maxRelGap);
    }
}

// This interface is now deprecated.
template<typename Real>
void IPM
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const IPMCtrl<Real>& ctrl )
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
    IPM( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

// TODO(poulson): Not use temporary regularization except in final iterations?
template<typename Real>
void EquilibratedIPM
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int degree = n;
    const Real epsilon = limits::Epsilon<Real>();
    const Grid& grid = problem.A.Grid();
    const int commRank = grid.Rank();
    Timer timer;

    const Real bNrm2 = FrobeniusNorm( problem.b );
    const Real cNrm2 = FrobeniusNorm( problem.c );
    const Real twoNormEstA =
      TwoNormEstimate( problem.A, ctrl.twoNormKrylovBasisSize );
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

    const Real xRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.xRegSmallLogEps);
    const Real yRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.yRegSmallLogEps);
    const Real zRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.zRegSmallLogEps);
    const Real xRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.xRegLargeLogEps);
    const Real yRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.yRegLargeLogEps);
    const Real zRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.zRegLargeLogEps);
    Real xRegLarge = xRegLarge0;
    Real yRegLarge = yRegLarge0;
    Real zRegLarge = zRegLarge0;
    Real xRegSmall = xRegSmall0;
    Real yRegSmall = yRegSmall0;
    Real zRegSmall = zRegSmall0;

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

    DistMultiVec<Real> regLarge(grid);
    if( ctrl.system == FULL_KKT )
    {
        regLarge.Resize( m+2*n, 1 );
        for( Int iLoc=0; iLoc<regLarge.LocalHeight(); ++iLoc )
        {
            const Int i = regLarge.GlobalRow(iLoc);
            if( i < n )        regLarge.SetLocal( iLoc, 0,  xRegLarge );
            else if( i < n+m ) regLarge.SetLocal( iLoc, 0, -yRegLarge );
            else               regLarge.SetLocal( iLoc, 0, -zRegLarge );
        }
    }
    else if( ctrl.system == AUGMENTED_KKT )
    {
        regLarge.Resize( n+m, 1 );
        for( Int iLoc=0; iLoc<regLarge.LocalHeight(); ++iLoc )
        {
            const Int i = regLarge.GlobalRow(iLoc);
            if( i < n ) regLarge.SetLocal( iLoc, 0,  xRegLarge );
            else        regLarge.SetLocal( iLoc, 0, -yRegLarge );
        }
    }
    else if( ctrl.system == NORMAL_KKT )
    {
        regLarge.Resize( m, 1 );
        Fill( regLarge, yRegLarge );
    }
    regLarge *= origTwoNormEst;

    Real muOld = 0.1;
    Real infeasError = 1;
    Real dimacsError = 1, dimacsErrorOld = 1;

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
    Int numIts=0;
    for( ; numIts<=ctrl.maxIts; ++numIts, dimacsErrorOld=dimacsError )
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
        mu = compRatio > ctrl.maxComplementRatio ? muOld : Min(mu,muOld);
        muOld = mu;

        pos_orth::NesterovTodd( solution.x, solution.z, w );
        const Real wMaxNorm = MaxNorm( w );

        // Check for convergence
        // =====================

        // Check the duality gap
        // ---------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real dualProd = Dot(solution.x,solution.z);
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relObjGap, relCompGap );

        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(-1), residual.primalEquality );
        const Real rbNrm2 = FrobeniusNorm( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);

        // || A^T y - z + c ||_2 / (1 + || c ||_2)
        // ---------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        residual.dualEquality -= solution.z;
        const Real rcNrm2 = FrobeniusNorm( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);

        // Now check the pieces
        // --------------------
        infeasError = Max(rbConv,rcConv);
        dimacsError = Max(maxRelGap,infeasError);
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
                 "  || primalInfeas ||_2 / (1 + || b ||_2) = ",
                 rbConv,"\n",Indent(),
                 "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",
                 rcConv,"\n",Indent(),
                 "  primal = ",primObj,"\n",Indent(),
                 "  dual   = ",dualObj,"\n",Indent(),
                 "  relative duality gap = ",maxRelGap);
        }

        const bool metTolerances =
          infeasError <= Pow(epsilon,ctrl.infeasibilityTolLogEps) &&
          relCompGap <= Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
          relObjGap <= Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);
        if( metTolerances )
        {
            if( dimacsError >= ctrl.minDimacsDecreaseRatio*dimacsErrorOld )
            {
                // We have met the tolerances and progress in the last iteration
                // was not significant.
                break;
            }
            else if( numIts == ctrl.maxIts )
            {
                // We have hit the iteration limit but can declare success.
                break;
            }
        }
        else if( numIts == ctrl.maxIts )
        {
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving tolerances");
        }

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
                ( problem.A,
                  Sqrt(xRegSmall),
                  Sqrt(yRegSmall),
                  Sqrt(zRegSmall),
                  solution.x, solution.z, JOrig, false );
                KKTRHS
                ( residual.dualEquality, residual.primalEquality,
                  residual.dualConic, solution.z, d );
            }
            else
            {
                AugmentedKKT
                ( problem.A, Sqrt(xRegSmall), Sqrt(yRegSmall),
                  solution.x, solution.z,
                  JOrig, false );
                AugmentedKKTRHS
                ( solution.x, residual.dualEquality, residual.primalEquality,
                  residual.dualConic, d );
            }
            J = JOrig;
            UpdateDiagonal( J, Real(1), regLarge );
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
            // TODO(poulson): Handle the equilibration consistently.
            /*
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
            */
            Ones( dInner, J.Height(), 1 );

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
            RegSolveInfo<Real> solveInfo;
            if( ctrl.twoStage )
            {
                solveInfo = reg_ldl::SolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d, ctrl.solveCtrl );
            }
            if( !solveInfo.metRequestedTol )
            {
                solveInfo = reg_ldl::RegularizedSolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress );
                if( !solveInfo.metRequestedTol )
                {
                    if( metTolerances )
                    {
                        break;
                    }
                    else
                    {
                        // TODO(poulson): Increase regularization and continue.
                        RuntimeError("Could not achieve tolerances");
                    }
                }
            }
            if( commRank == 0 && ctrl.time )
                Output("Affine: ",timer.Stop()," secs");
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
            ( problem.A, Sqrt(xRegSmall), Sqrt(yRegSmall),
              solution.x, solution.z, J, false );
            NormalKKTRHS
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
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
            auto solveInfo = reg_ldl::RegularizedSolveAfter
            ( J, regLarge, sparseLDLFact, affineCorrection.y,
              ctrl.solveCtrl.relTol,
              ctrl.solveCtrl.maxRefineIts,
              ctrl.solveCtrl.progress,
              ctrl.solveCtrl.time );
            if( !solveInfo.metRequestedTol )
            {
                if( metTolerances )
                {
                    break;
                }
                else
                {
                    // TODO(poulson): Increase regularization and continue.
                    RuntimeError("Could not achieve tolerances");
                }
            }
            if( commRank == 0 && ctrl.time )
                Output("Affine: ",timer.Stop()," secs");
            ExpandNormalSolution
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
              residual.dualEquality, residual.dualConic,
              affineCorrection.x, affineCorrection.y, affineCorrection.z );
        }

        /*
        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Gemv
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            Axpy( -yRegSmall, affineCorrection.y, error.primalEquality );
            Real dxErrorNrm2 = FrobeniusNorm( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Gemv
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Axpy( xRegSmall, affineCorrection.x, error.dualEquality );
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
                ("|| dxError ||_2 / (1 + || primalInfeas ||_2) = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || dualInfeas   ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rmuNrm2));
        }
        */

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
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.compositeNewton )
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
            if( commRank == 0 && ctrl.time )
                timer.Start();
            RegSolveInfo<Real> solveInfo;
            if( ctrl.twoStage )
            {
                solveInfo = reg_ldl::SolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d, ctrl.solveCtrl );
            }
            if( !solveInfo.metRequestedTol )
            {
                solveInfo = reg_ldl::RegularizedSolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress );
                if( !solveInfo.metRequestedTol )
                {
                    if( metTolerances )
                    {
                        break;
                    }
                    else
                    {
                        // TODO(poulson): Increase regularization and continue.
                        RuntimeError("Could not achieve tolerances");
                    }
                }
            }
            if( commRank == 0 && ctrl.time )
                Output("Corrector: ",timer.Stop()," secs");
            ExpandSolution( m, n, d, correction.x, correction.y, correction.z );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            AugmentedKKTRHS
            ( solution.x, residual.dualEquality, residual.primalEquality,
              residual.dualConic, d );
            if( commRank == 0 && ctrl.time )
                timer.Start();
            RegSolveInfo<Real> solveInfo;
            if( ctrl.twoStage )
            {
                solveInfo = reg_ldl::SolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d, ctrl.solveCtrl );
            }
            if( !solveInfo.metRequestedTol )
            {
                auto solveInfo = reg_ldl::RegularizedSolveAfter
                ( JOrig, regLarge, dInner, sparseLDLFact, d,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress );
                if( !solveInfo.metRequestedTol )
                {
                    if( metTolerances )
                    {
                        break;
                    }
                    else
                    {
                        // TODO(poulson): Increase regularization and continue.
                        RuntimeError("Could not achieve tolerances");
                    }
                }
            }
            if( commRank == 0 && ctrl.time )
                Output("Corrector: ",timer.Stop()," secs");
            ExpandAugmentedSolution
            ( solution.x, solution.z, residual.dualConic, d,
              correction.x, correction.y, correction.z );
        }
        else
        {
            NormalKKTRHS
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
              residual.dualEquality, residual.primalEquality,
              residual.dualConic, correction.y );
            if( commRank == 0 && ctrl.time )
                timer.Start();
            auto solveInfo = reg_ldl::RegularizedSolveAfter
            ( J, regLarge, sparseLDLFact, correction.y,
              ctrl.solveCtrl.relTol,
              ctrl.solveCtrl.maxRefineIts,
              ctrl.solveCtrl.progress,
              ctrl.solveCtrl.time );
            if( !solveInfo.metRequestedTol )
            {
                if( metTolerances )
                {
                    break;
                }
                else
                {
                    // TODO(poulson): Increase regularization and continue.
                    RuntimeError("Could not achieve tolerances");
                }
            }
            if( commRank == 0 && ctrl.time )
                Output("Corrector: ",timer.Stop()," secs");
            ExpandNormalSolution
            ( problem.A, Sqrt(xRegSmall), solution.x, solution.z,
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
            if( metTolerances )
            {
                break;
            }
            else
            {
                RuntimeError("Could not achieve tolerances");
            }
        }
    }
    SetIndent( indent );
}

template<typename Real>
void IPM
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
        EquilibratedIPM( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedIPM( problem, solution, ctrl );
    }
    if( ctrl.print )
    {
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y);
        const Real dualProd = Dot(solution.x,solution.z);
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relObjGap, relCompGap );

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
         "  relative duality gap = ",maxRelGap);
    }
}

// This interface is now deprecated.
template<typename Real>
void IPM
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const IPMCtrl<Real>& ctrl )
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
    IPM( problem, solution, ctrl );
    x = solution.x;
    y = solution.y;
    z = solution.z;
}

#define PROTO(Real) \
  template void IPM \
  ( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
          DirectLPSolution<Matrix<Real>>& solution, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
          DirectLPSolution<DistMatrix<Real>>& solution, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
          DirectLPSolution<Matrix<Real>>& solution, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
          DirectLPSolution<DistMultiVec<Real>>& solution, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const IPMCtrl<Real>& ctrl );

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
