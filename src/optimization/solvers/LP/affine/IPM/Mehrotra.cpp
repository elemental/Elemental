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

namespace lp {

template<typename Real>
Real RelativeDualityGap
( const Real& primalObj, const Real& dualObj, const Real& dualityProduct )
{
    EL_DEBUG_CSE
    Real relGap;
    if( primalObj < Real(0) )
        relGap = dualityProduct / -primalObj;
    else if( dualObj > Real(0) )
        relGap = dualityProduct / dualObj;
    else
        relGap = 2; // 200% error if the signs differ wrong.
    return relGap;
}

namespace affine {

// The following solves a pair of linear programs in "affine" conic form:
//
//   min c^T x
//   s.t. A x = b, G x + s = h, s >= 0,
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z >= 0,
//
// as opposed to the more specific "direct" conic form:
//
//   min c^T x
//   s.t. A x = b, x >= 0,
//
//   max -b^T y
//   s.t. A^T y - z + c = 0, z >= 0,
//
// which corresponds to G = -I and h = 0, using a Mehrotra Predictor-Corrector
// scheme.
//
// We make use of the regularized Lagrangian
//
//   L(x,s;y,z) = c^T x + y^T (A x - b) + z^T (G x + s - h) + mu Phi(s) +
//                (1/2) gamma_x^2 || x ||_2^2 -
//                (1/2) gamma_y^2 || y ||_2^2 -
//                (1/2) gamma_z^2 || z ||_2^2,
//
// where we note that the two-norm regularization is positive for the primal
// variable x and *negative* for the dual variables y and z. There is not yet
// any regularization on the primal slack variable s (though it may be
// investigated in the future).
//
// The subsequent first-order optimality conditions for x, y, and z become
//
//   Delta_x L = c + A^T y + G^T z + gamma_x^2 x = 0,
//   Delta_y L = A x - b - gamma_y^2 y = 0,
//   Delta_z L = G x + s - h - gamma_z^2 z = 0.
//
// These can be arranged into the symmetric quasi-definite form
//
//   | gamma_x^2 I,      A^T,          G^T     | | x | = | -c  |
//   |      A,      -gamma_y^2 I,       0      | | y |   |  b  |
//   |      G,            0,      -gamma_z^2 I | | z |   | h-s |.
//

template<typename Real,class MatrixType,class VectorType>
Real PrimalObjective
( const AffineLPProblem<MatrixType,VectorType>& problem,
  const AffineLPSolution<VectorType>& solution )
{
    EL_DEBUG_CSE
    const Real primalObjective = Dot(problem.c,solution.x);
    return primalObjective;
}

template<typename Real,class MatrixType,class VectorType>
Real DualObjective
( const AffineLPProblem<MatrixType,VectorType>& problem,
  const AffineLPSolution<VectorType>& solution )
{
    EL_DEBUG_CSE
    const Real dualObjective = -Dot(problem.b,solution.y) -
      Dot(problem.h,solution.z);
    return dualObjective;
}

/*
template<typename Real,class MatrixType,class VectorType>
Real RelativeDualityGap
( const AffineLPProblem<MatrixType,VectorType>& problem,
  const AffineLPSolution<VectorType>& solution )
{
    EL_DEBUG_CSE
    const Real primalObj = PrimalObjective<Real>( problem, solution );
    const Real dualObj = DualObjective<Real>( problem, solution );
    const Real dualProd = Dot(solution.s,solution.z);
    const Real relGap = RelativeDualityGap( primalObj, dualObj, dualProd );
    return relGap;
}
*/

template<typename Real>
struct DenseAffineLPEquilibration
{
    Real sScale;
    Real zScale;
    Matrix<Real> rowScaleA;
    Matrix<Real> rowScaleG;
    Matrix<Real> colScale;
};
template<typename Real>
struct DistDenseAffineLPEquilibration
{
    Real sScale;
    Real zScale;
    DistMatrix<Real,MC,STAR> rowScaleA;
    DistMatrix<Real,MC,STAR> rowScaleG;
    DistMatrix<Real,MR,STAR> colScale;
};
template<typename Real>
struct SparseAffineLPEquilibration
{
    Real sScale;
    Real zScale;
    Matrix<Real> rowScaleA;
    Matrix<Real> rowScaleG;
    Matrix<Real> colScale;
};
template<typename Real>
struct DistSparseAffineLPEquilibration
{
    Real sScale;
    Real zScale;
    DistMultiVec<Real> rowScaleA;
    DistMultiVec<Real> rowScaleG;
    DistMultiVec<Real> colScale;
};

template<typename Real>
void Equilibrate
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const AffineLPSolution<Matrix<Real>>& solution,
        AffineLPProblem<Matrix<Real>,Matrix<Real>>& equilibratedProblem,
        AffineLPSolution<Matrix<Real>>& equilibratedSolution,
        DenseAffineLPEquilibration<Real>& equilibration,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    equilibratedProblem = problem;
    equilibratedSolution = solution;

    // Equilibrate the LP by diagonally scaling [A;G]
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      equilibration.rowScaleA,
      equilibration.rowScaleG,
      equilibration.colScale,
      ctrl.print );

    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScaleA, equilibratedProblem.b );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedProblem.h );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.colScale, equilibratedProblem.c );
    if( ctrl.primalInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.x );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedSolution.s );
    }
    if( ctrl.dualInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScaleA, equilibratedSolution.y );
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedSolution.z );
    }

    // Rescale max(|| b ||_max,|| h ||_max) to roughly one.
    equilibration.sScale =
      Max(Max(MaxNorm(equilibratedProblem.b),MaxNorm(equilibratedProblem.h)),
        Real(1e-6));
    equilibratedProblem.b *= Real(1)/equilibration.sScale;
    equilibratedProblem.h *= Real(1)/equilibration.sScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.sScale;
        equilibratedSolution.s *= Real(1)/equilibration.sScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.zScale = Max(MaxNorm(equilibratedProblem.c),Real(1e-6));
    equilibratedProblem.c *= Real(1)/equilibration.zScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.zScale;
        equilibratedSolution.z *= Real(1)/equilibration.zScale;
    }

    if( ctrl.print )
    {
        Output("sScale=",equilibration.sScale);
        Output("zScale=",equilibration.zScale);
        Output("|| rowScaleA ||_2 = ",FrobeniusNorm(equilibration.rowScaleA));
        Output("|| rowScaleG ||_2 = ",FrobeniusNorm(equilibration.rowScaleG));
        Output("|| colScale  ||_2 = ",FrobeniusNorm(equilibration.colScale));
    }
}

template<typename Real>
void Equilibrate
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const AffineLPSolution<DistMatrix<Real>>& solution,
        AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& equilibratedProblem,
        AffineLPSolution<DistMatrix<Real>>& equilibratedSolution,
        DistDenseAffineLPEquilibration<Real>& equilibration,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    equilibratedProblem = problem;
    equilibratedSolution = solution;

    // Equilibrate the LP by diagonally scaling [A;G]
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      equilibration.rowScaleA,
      equilibration.rowScaleG,
      equilibration.colScale,
      ctrl.print );

    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScaleA, equilibratedProblem.b );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedProblem.h );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.colScale, equilibratedProblem.c );
    if( ctrl.primalInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.x );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedSolution.s );
    }
    if( ctrl.dualInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScaleA, equilibratedSolution.y );
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedSolution.z );
    }

    // Rescale max(|| b ||_max,|| h ||_max) to roughly one.
    equilibration.sScale =
      Max(Max(MaxNorm(equilibratedProblem.b),MaxNorm(equilibratedProblem.h)),
        Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.sScale;
    equilibratedProblem.h *= Real(1)/equilibration.sScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.sScale;
        equilibratedSolution.s *= Real(1)/equilibration.sScale;
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
        Output("sScale=",equilibration.sScale);
        Output("zScale=",equilibration.zScale);
        Output("|| rowScaleA ||_2 = ",FrobeniusNorm(equilibration.rowScaleA));
        Output("|| rowScaleG ||_2 = ",FrobeniusNorm(equilibration.rowScaleG));
        Output("|| colScale  ||_2 = ",FrobeniusNorm(equilibration.colScale));
    }
}

template<typename Real>
void Equilibrate
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const AffineLPSolution<Matrix<Real>>& solution,
        AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& equilibratedProblem,
        AffineLPSolution<Matrix<Real>>& equilibratedSolution,
        SparseAffineLPEquilibration<Real>& equilibration,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    equilibratedProblem = problem;
    equilibratedSolution = solution;

    // Equilibrate the LP by diagonally scaling [A;G]
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      equilibration.rowScaleA,
      equilibration.rowScaleG,
      equilibration.colScale,
      ctrl.print );

    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScaleA, equilibratedProblem.b );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedProblem.h );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.colScale, equilibratedProblem.c );
    if( ctrl.primalInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.x );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedSolution.s );
    }
    if( ctrl.dualInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScaleA, equilibratedSolution.y );
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedSolution.z );
    }

    // Rescale max(|| b ||_max,|| h ||_max) to roughly one.
    equilibration.sScale =
      Max(Max(MaxNorm(equilibratedProblem.b),MaxNorm(equilibratedProblem.h)),
        Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.sScale;
    equilibratedProblem.h *= Real(1)/equilibration.sScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.sScale;
        equilibratedSolution.s *= Real(1)/equilibration.sScale;
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
        Output("sScale=",equilibration.sScale);
        Output("zScale=",equilibration.zScale);
        Output("|| rowScaleA ||_2 = ",FrobeniusNorm(equilibration.rowScaleA));
        Output("|| rowScaleG ||_2 = ",FrobeniusNorm(equilibration.rowScaleG));
        Output("|| colScale  ||_2 = ",FrobeniusNorm(equilibration.colScale));
    }
}

template<typename Real>
void Equilibrate
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const AffineLPSolution<DistMultiVec<Real>>& solution,
        AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>&
          equilibratedProblem,
        AffineLPSolution<DistMultiVec<Real>>& equilibratedSolution,
        DistSparseAffineLPEquilibration<Real>& equilibration,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = problem.A.Grid();
    ForceSimpleAlignments( equilibratedProblem, grid );
    ForceSimpleAlignments( equilibratedSolution, grid );

    equilibratedProblem = problem;
    equilibratedSolution = solution;

    // Equilibrate the LP by diagonally scaling [A;G]
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      equilibration.rowScaleA,
      equilibration.rowScaleG,
      equilibration.colScale,
      ctrl.print );

    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScaleA, equilibratedProblem.b );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedProblem.h );
    DiagonalSolve
    ( LEFT, NORMAL, equilibration.colScale, equilibratedProblem.c );
    if( ctrl.primalInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.colScale, equilibratedSolution.x );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedSolution.s );
    }
    if( ctrl.dualInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScaleA, equilibratedSolution.y );
        DiagonalScale
        ( LEFT, NORMAL, equilibration.rowScaleG, equilibratedSolution.z );
    }

    // Rescale max(|| b ||_max,|| h ||_max) to roughly one.
    equilibration.sScale =
      Max(Max(MaxNorm(equilibratedProblem.b),MaxNorm(equilibratedProblem.h)),
        Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.sScale;
    equilibratedProblem.h *= Real(1)/equilibration.sScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.sScale;
        equilibratedSolution.s *= Real(1)/equilibration.sScale;
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
        Output("sScale=",equilibration.sScale);
        Output("zScale=",equilibration.zScale);
        Output("|| rowScaleA ||_2 = ",FrobeniusNorm(equilibration.rowScaleA));
        Output("|| rowScaleG ||_2 = ",FrobeniusNorm(equilibration.rowScaleG));
        Output("|| colScale  ||_2 = ",FrobeniusNorm(equilibration.colScale));
    }
}

template<typename Real>
void UndoEquilibration
( const AffineLPSolution<Matrix<Real>>& equilibratedSolution,
  const DenseAffineLPEquilibration<Real>& equilibration,
        AffineLPSolution<Matrix<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;
    solution.x *= equilibration.sScale;
    solution.s *= equilibration.sScale;
    solution.y *= equilibration.zScale;
    solution.z *= equilibration.zScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.colScale,  solution.x );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScaleA, solution.y );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScaleG, solution.z );
    DiagonalScale( LEFT, NORMAL, equilibration.rowScaleG, solution.s );
}

template<typename Real>
void UndoEquilibration
( const AffineLPSolution<DistMatrix<Real>>& equilibratedSolution,
  const DistDenseAffineLPEquilibration<Real>& equilibration,
        AffineLPSolution<DistMatrix<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;
    solution.x *= equilibration.sScale;
    solution.s *= equilibration.sScale;
    solution.y *= equilibration.zScale;
    solution.z *= equilibration.zScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.colScale,  solution.x );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScaleA, solution.y );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScaleG, solution.z );
    DiagonalScale( LEFT, NORMAL, equilibration.rowScaleG, solution.s );
}

template<typename Real>
void UndoEquilibration
( const AffineLPSolution<Matrix<Real>>& equilibratedSolution,
  const SparseAffineLPEquilibration<Real>& equilibration,
        AffineLPSolution<Matrix<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;
    solution.x *= equilibration.sScale;
    solution.s *= equilibration.sScale;
    solution.y *= equilibration.zScale;
    solution.z *= equilibration.zScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.colScale,  solution.x );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScaleA, solution.y );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScaleG, solution.z );
    DiagonalScale( LEFT, NORMAL, equilibration.rowScaleG, solution.s );
}

template<typename Real>
void UndoEquilibration
( const AffineLPSolution<DistMultiVec<Real>>& equilibratedSolution,
  const DistSparseAffineLPEquilibration<Real>& equilibration,
        AffineLPSolution<DistMultiVec<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;
    solution.x *= equilibration.sScale;
    solution.s *= equilibration.sScale;
    solution.y *= equilibration.zScale;
    solution.z *= equilibration.zScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.colScale,  solution.x );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScaleA, solution.y );
    DiagonalSolve( LEFT, NORMAL, equilibration.rowScaleG, solution.z );
    DiagonalScale( LEFT, NORMAL, equilibration.rowScaleG, solution.s );
}

template<typename Real>
void EquilibratedMehrotra
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Int degree = k;

    const Real bNrm2 = Nrm2( problem.b );
    const Real cNrm2 = Nrm2( problem.c );
    const Real hNrm2 = Nrm2( problem.h );
    if( ctrl.print )
    {
        const Real ANrm1 = OneNorm( problem.A );
        const Real GNrm1 = OneNorm( problem.G );
        Output("|| c ||_2 = ",cNrm2);
        Output("|| A ||_1 = ",ANrm1);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| G ||_1 = ",GNrm1);
        Output("|| h ||_2 = ",hNrm2);
    }

    Initialize
    ( problem, solution,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift );

    Real dimacsError = 1, dimacsErrorOld = 1;
    Real infeasError = 1;
    Matrix<Real> J, d;
    Matrix<Real> dSub;
    Permutation p;
    auto attemptToFactor = [&]()
      {
        EL_DEBUG_ONLY(auto callStack = CopyCallStack())
        try { LDL( J, dSub, p, false ); }
        catch(...)
        {
            if( dimacsError > ctrl.minTol )
                RuntimeError
                ("Unable to achieve minimum tolerance ",ctrl.minTol);
            EL_DEBUG_ONLY(SetCallStack(callStack))
            return false;
        }
        return true;
      };
    auto attemptToSolve = [&]( Matrix<Real>& rhs )
      {
        EL_DEBUG_ONLY(auto callStack = CopyCallStack())
        try { ldl::SolveAfter( J, dSub, p, rhs, false ); }
        catch(...)
        {
            if( dimacsError > ctrl.minTol )
                RuntimeError
                ("Unable to achieve minimum tolerance ",ctrl.minTol);
            EL_DEBUG_ONLY(SetCallStack(callStack))
            return false;
        }
        return true;
      };

    AffineLPSolution<Matrix<Real>> affineCorrection, correction;
    AffineLPResidual<Matrix<Real>> residual, error;
    const Int indent = PushIndent();
    Int numIts=0;
    for( ; numIts<=ctrl.maxIts; ++numIts, dimacsErrorOld=dimacsError )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = pos_orth::NumOutside( solution.s );
        const Int zNumNonPos = pos_orth::NumOutside( solution.z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure
        // ===========================
        const Real dualProd = Dot( solution.s, solution.z );
        const Real mu = dualProd / k;

        // Check for convergence
        // =====================
        const Real primObj = PrimalObjective<Real>( problem, solution );
        const Real dualObj = DualObjective<Real>( problem, solution );
        const Real relGap = RelativeDualityGap( primObj, dualObj, dualProd );
        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
        const Real primalInfeasNrm2 = Nrm2( residual.primalEquality );
        const Real primalInfeasNrm2Rel = primalInfeasNrm2 / (1+bNrm2);
        // || c + A^T y + G^T z ||_2 / (1 + || c ||_2)
        // -------------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        Gemv
        ( TRANSPOSE, Real(1), problem.G, solution.z,
          Real(1), residual.dualEquality );
        const Real dualInfeasNrm2 = Nrm2( residual.dualEquality );
        const Real dualInfeasNrm2Rel = dualInfeasNrm2 / (1+cNrm2);
        // || G x + s - h ||_2 / (1 + || h ||_2)
        // -------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Gemv
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        residual.primalConic += solution.s;
        const Real conicInfeasNrm2 = Nrm2( residual.primalConic );
        const Real conicInfeasNrm2Rel = conicInfeasNrm2 / (1+hNrm2);
        // Now check the pieces
        // --------------------
        const Real equalityError = Max(primalInfeasNrm2Rel,dualInfeasNrm2Rel);
        infeasError = Max(equalityError,conicInfeasNrm2Rel);
        dimacsError = Max(infeasError,relGap);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( solution.x );
            const Real yNrm2 = Nrm2( solution.y );
            const Real zNrm2 = Nrm2( solution.z );
            const Real sNrm2 = Nrm2( solution.s );
            Output
            ("iter ",numIts,":\n",Indent(),
             "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
             "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
             "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
             "  ||  s  ||_2 = ",sNrm2,"\n",Indent(),
             "  || primalInfeas ||_2 / (1 + || b ||_2) = ",
             primalInfeasNrm2Rel,"\n",Indent(),
             "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",
             dualInfeasNrm2Rel,"\n",Indent(),
             "  || conicInfeas  ||_2 / (1 + || h ||_2) = ",
             conicInfeasNrm2Rel,"\n",Indent(),
             "  scaled primal = ",primObj,"\n",Indent(),
             "  scaled dual   = ",dualObj,"\n",Indent(),
             "  scaled relative duality gap = ",relGap);
        }
        if( dimacsError <= ctrl.targetTol )
            break;
        // Exit if progress has stalled and we are sufficiently accurate.
        if( dimacsError <= ctrl.minTol &&
            dimacsError >= ctrl.minDimacsDecreaseRatio*dimacsErrorOld )
            break;
        if( numIts == ctrl.maxIts && dimacsError > ctrl.minTol )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := s o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

        // Construct the KKT system
        // ------------------------
        KKT( problem.A, problem.G, solution.s, solution.z, J );
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z,
          d );

        // Solve for the direction
        // -----------------------
        if( !attemptToFactor() )
            break;
        if( !attemptToSolve(d) )
            break;
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, solution.z,
          affineCorrection.x,
          affineCorrection.y,
          affineCorrection.z,
          affineCorrection.s );

        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Gemv
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            const Real dxErrorNrm2 = Nrm2( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Gemv
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Gemv
            ( TRANSPOSE, Real(1), problem.G, affineCorrection.z,
              Real(1), error.dualEquality );
            const Real dyErrorNrm2 = Nrm2( error.dualEquality );

            error.primalConic = residual.primalConic;
            Gemv
            ( NORMAL, Real(1), problem.G, affineCorrection.x,
              Real(1), error.primalConic );
            error.primalConic += affineCorrection.s;
            const Real dzErrorNrm2 = Nrm2( error.primalConic );

            // TODO(poulson): dmuError

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+primalInfeasNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+dualInfeasNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+conicInfeasNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri =
          pos_orth::MaxStep( solution.s, affineCorrection.s, Real(1) );
        Real alphaAffDual =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and ds are used as temporaries
        correction.s = solution.s;
        correction.z = solution.z;
        Axpy( alphaAffPri,  affineCorrection.s, correction.s );
        Axpy( alphaAffDual, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.s,correction.z) / degree;
        if( ctrl.print )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma =
          ctrl.centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dsAff o dzAff
            // ---------------------
            // NOTE: Using correction.z as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.s, correction.z );
            residual.dualConic += correction.z;
        }

        // Construct the new KKT RHS
        // -------------------------
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        // Solve for the direction
        // -----------------------
        if( !attemptToSolve(d) )
            break;
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, solution.z,
          correction.x, correction.y, correction.z, correction.s );
        // TODO(poulson): Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri =
          pos_orth::MaxStep( solution.s, correction.s, 1/ctrl.maxStepRatio );
        Real alphaDual =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaPri,  correction.s, solution.s );
        Axpy( alphaDual, correction.y, solution.y );
        Axpy( alphaDual, correction.z, solution.z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( dimacsError <= ctrl.minTol )
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
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.outerEquil )
    {
        AffineLPProblem<Matrix<Real>,Matrix<Real>> equilibratedProblem;
        AffineLPSolution<Matrix<Real>> equilibratedSolution;
        DenseAffineLPEquilibration<Real> equilibration;
        Equilibrate
        ( problem, solution,
          equilibratedProblem, equilibratedSolution,
          equilibration, ctrl );
        EquilibratedMehrotra( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedMehrotra( problem, solution, ctrl );
    }
}

template<typename Real>
void EquilibratedMehrotra
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    // Ensure that the problem is aligned.
    EL_DEBUG_ONLY(
      if( !SimpleAlignments(problem) )
          LogicError("Problem specification matrices were not simply aligned");
    )

    // Ensure that the solution vectors are aligned.
    EL_DEBUG_ONLY(
      if( !SimpleAlignments(solution) )
          LogicError("Solution matrices were not simply aligned");
    )

    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Int degree = k;
    const Grid& grid = problem.A.Grid();
    const int commRank = grid.Rank();

    const Real bNrm2 = Nrm2( problem.b );
    const Real cNrm2 = Nrm2( problem.c );
    const Real hNrm2 = Nrm2( problem.h );
    if( ctrl.print )
    {
        const Real ANrm1 = OneNorm( problem.A );
        const Real GNrm1 = OneNorm( problem.G );
        if( commRank == 0 )
        {
            Output("|| c ||_2 = ",cNrm2);
            Output("|| A ||_1 = ",ANrm1);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| G ||_1 = ",GNrm1);
            Output("|| h ||_2 = ",hNrm2);
        }
    }

    Initialize
    ( problem, solution,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift );

    Real dimacsError = 1, dimacsErrorOld = 1;
    Real infeasError = 1;
    DistMatrix<Real> J(grid), d(grid);
    DistMatrix<Real> dSub(grid);
    DistPermutation p(grid);
    auto attemptToFactor = [&]()
      {
        EL_DEBUG_ONLY(auto callStack = CopyCallStack())
        try { LDL( J, dSub, p, false ); }
        catch(...)
        {
            if( dimacsError > ctrl.minTol )
                RuntimeError
                ("Unable to achieve minimum tolerance ",ctrl.minTol);
            EL_DEBUG_ONLY(SetCallStack(callStack))
            return false;
        }
        return true;
      };
    auto attemptToSolve = [&]( DistMatrix<Real>& rhs )
      {
        EL_DEBUG_ONLY(auto callStack = CopyCallStack())
        try { ldl::SolveAfter( J, dSub, p, rhs, false ); }
        catch(...)
        {
            if( dimacsError > ctrl.minTol )
                RuntimeError
                ("Unable to achieve minimum tolerance ",ctrl.minTol);
            EL_DEBUG_ONLY(SetCallStack(callStack))
            return false;
        }
        return true;
      };

    AffineLPResidual<DistMatrix<Real>> residual, error;
    AffineLPSolution<DistMatrix<Real>> affineCorrection, correction;
    ForceSimpleAlignments( residual, grid );
    ForceSimpleAlignments( error, grid );
    ForceSimpleAlignments( affineCorrection, grid );
    ForceSimpleAlignments( correction, grid );

    const Int indent = PushIndent();
    Int numIts=0;
    for( ; numIts<=ctrl.maxIts; ++numIts, dimacsErrorOld=dimacsError )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = pos_orth::NumOutside( solution.s );
        const Int zNumNonPos = pos_orth::NumOutside( solution.z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure
        // ===========================
        const Real dualProd = Dot( solution.s, solution.z );
        const Real mu = dualProd / k;

        // Check for convergence
        // =====================
        const Real primObj = PrimalObjective<Real>( problem, solution );
        const Real dualObj = DualObjective<Real>( problem, solution );
        const Real relGap = RelativeDualityGap( primObj, dualObj, dualProd );
        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
        const Real primalInfeasNrm2 = Nrm2( residual.primalEquality );
        const Real primalInfeasNrm2Rel = primalInfeasNrm2 / (1+bNrm2);
        // || c + A^T y + G^T z ||_2 / (1 + || c ||_2)
        // -------------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        Gemv
        ( TRANSPOSE, Real(1), problem.G, solution.z,
          Real(1), residual.dualEquality );
        const Real dualInfeasNrm2 = Nrm2( residual.dualEquality );
        const Real dualInfeasNrm2Rel = dualInfeasNrm2 / (1+cNrm2);
        // || G x + s - h ||_2 / (1 + || h ||_2)
        // -------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Gemv
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        residual.primalConic += solution.s;
        const Real conicInfeasNrm2 = Nrm2( residual.primalConic );
        const Real conicInfeasNrm2Rel = conicInfeasNrm2 / (1+hNrm2);
        // Now check the pieces
        // --------------------
        const Real equalityError = Max(primalInfeasNrm2Rel,dualInfeasNrm2Rel);
        infeasError = Max(equalityError,conicInfeasNrm2Rel);
        dimacsError = Max(infeasError,relGap);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( solution.x );
            const Real yNrm2 = Nrm2( solution.y );
            const Real zNrm2 = Nrm2( solution.z );
            const Real sNrm2 = Nrm2( solution.s );
            if( commRank == 0 )
                Output
                ("iter ",numIts,":\n",Indent(),
                 "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
                 "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
                 "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
                 "  ||  s  ||_2 = ",sNrm2,"\n",Indent(),
                 "  || primalInfeas ||_2 / (1 + || b ||_2) = ",
                 primalInfeasNrm2Rel,"\n",Indent(),
                 "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",
                 dualInfeasNrm2Rel,"\n",Indent(),
                 "  || conicInfeas  ||_2 / (1 + || h ||_2) = ",
                 conicInfeasNrm2Rel,"\n",Indent(),
                 "  scaled primal = ",primObj,"\n",Indent(),
                 "  scaled dual   = ",dualObj,"\n",Indent(),
                 "  scaled relative gap = ",relGap);
        }
        if( dimacsError <= ctrl.targetTol )
            break;
        // Exit if progress has stalled and we are sufficiently accurate.
        if( dimacsError <= ctrl.minTol &&
            dimacsError >= ctrl.minDimacsDecreaseRatio*dimacsErrorOld )
            break;
        if( numIts == ctrl.maxIts && dimacsError > ctrl.minTol )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================
        // r_mu := s o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );
        // Construct the KKT system
        // ------------------------
        KKT( problem.A, problem.G, solution.s, solution.z, J );
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        // Solve for the proposed step
        // ---------------------------
        if( !attemptToFactor() )
            break;
        if( !attemptToSolve(d) )
            break;
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, solution.z,
          affineCorrection.x,
          affineCorrection.y,
          affineCorrection.z,
          affineCorrection.s );

        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Gemv
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            const Real dxErrorNrm2 = Nrm2( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Gemv
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Gemv
            ( TRANSPOSE, Real(1), problem.G, affineCorrection.z,
              Real(1), error.dualEquality );
            const Real dyErrorNrm2 = Nrm2( error.dualEquality );

            error.primalConic = residual.primalConic;
            Gemv
            ( NORMAL, Real(1), problem.G, affineCorrection.x,
              Real(1), error.primalConic );
            error.primalConic += affineCorrection.s;
            const Real dzErrorNrm2 = Nrm2( error.primalConic );

            // TODO(poulson): error.dualConic

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
                 dxErrorNrm2/(1+primalInfeasNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+dualInfeasNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+conicInfeasNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri =
          pos_orth::MaxStep( solution.s, affineCorrection.s, Real(1) );
        Real alphaAffDual =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and ds are used as temporaries
        correction.s = solution.s;
        correction.z = solution.z;
        Axpy( alphaAffPri,  affineCorrection.s, correction.s );
        Axpy( alphaAffDual, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.s,correction.z) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma =
          ctrl.centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dsAff o dzAff
            // ---------------------
            // NOTE: dz is used as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.s, correction.z );
            residual.dualConic += correction.z;
        }

        // Construct the new KKT RHS
        // -------------------------
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        // Solve for the direction
        // -----------------------
        if( !attemptToSolve(d) )
            break;
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, solution.z,
          correction.x, correction.y, correction.z, correction.s );
        // TODO(poulson): Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri =
          pos_orth::MaxStep( solution.s, correction.s, 1/ctrl.maxStepRatio );
        Real alphaDual =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaPri,  correction.s, solution.s );
        Axpy( alphaDual, correction.y, solution.y );
        Axpy( alphaDual, correction.z, solution.z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( dimacsError <= ctrl.minTol )
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
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Grid& grid = problem.A.Grid();
    if( ctrl.outerEquil )
    {
        AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>> equilibratedProblem;
        AffineLPSolution<DistMatrix<Real>> equilibratedSolution;
        DistDenseAffineLPEquilibration<Real> equilibration;
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
            AffineLPSolution<DistMatrix<Real>> alignedSolution;
            ForceSimpleAlignments( alignedSolution, grid );
            alignedSolution = solution;
            EquilibratedMehrotra( problem, alignedSolution, ctrl );
            solution = alignedSolution;
        }
        else if( SimpleAlignments(solution) )
        {
            AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>> alignedProblem;
            ForceSimpleAlignments( alignedProblem, grid );
            CopyOrViewHelper( problem.c, alignedProblem.c );
            CopyOrViewHelper( problem.A, alignedProblem.A );
            CopyOrViewHelper( problem.b, alignedProblem.b );
            CopyOrViewHelper( problem.G, alignedProblem.G );
            CopyOrViewHelper( problem.h, alignedProblem.h );
            EquilibratedMehrotra( alignedProblem, solution, ctrl );
        }
        else
        {
            AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>> alignedProblem;
            ForceSimpleAlignments( alignedProblem, grid );
            CopyOrViewHelper( problem.c, alignedProblem.c );
            CopyOrViewHelper( problem.A, alignedProblem.A );
            CopyOrViewHelper( problem.b, alignedProblem.b );
            CopyOrViewHelper( problem.G, alignedProblem.G );
            CopyOrViewHelper( problem.h, alignedProblem.h );
            AffineLPSolution<DistMatrix<Real>> alignedSolution;
            ForceSimpleAlignments( alignedSolution, grid );
            alignedSolution = solution;
            EquilibratedMehrotra( alignedProblem, alignedSolution, ctrl );
            solution = alignedSolution;
        }
    }
}

template<typename Real>
void EquilibratedMehrotra
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& origProblem,
        SparseAffineLPEquilibration<Real>& equilibration,
        AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Int degree = k;
    bool twoStage = ctrl.twoStage;

    if( ctrl.regIncreaseFactor <= Real(1) )
        LogicError("Regularization increase factor must be at least 1");

    const Real bNrm2 = Nrm2( problem.b );
    const Real cNrm2 = Nrm2( problem.c );
    const Real hNrm2 = Nrm2( problem.h );
    const Real twoNormEstA = TwoNormEstimate( problem.A, ctrl.basisSize );
    const Real twoNormEstG = TwoNormEstimate( problem.G, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstA + twoNormEstG + 1;
    if( ctrl.print )
    {
        Output("|| A ||_2 estimate: ",twoNormEstA);
        Output("|| G ||_2 estimate: ",twoNormEstG);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
        Output("|| h ||_2 = ",hNrm2);
    }

    const Real xRegSmall0 = origTwoNormEst*ctrl.xRegSmall;
    const Real yRegSmall0 = origTwoNormEst*ctrl.yRegSmall;
    const Real zRegSmall0 = origTwoNormEst*ctrl.zRegSmall;
    const Real xRegLarge0 = origTwoNormEst*ctrl.xRegLarge;
    const Real yRegLarge0 = origTwoNormEst*ctrl.yRegLarge;
    const Real zRegLarge0 = origTwoNormEst*ctrl.zRegLarge;
    Real xRegLarge = xRegLarge0;
    Real yRegLarge = yRegLarge0;
    Real zRegLarge = zRegLarge0;
    Real xRegSmall = xRegSmall0;
    Real yRegSmall = yRegSmall0;
    Real zRegSmall = zRegSmall0;
    // Once the permanent regularization is sufficiently large, we no longer
    // need to have separate 'temporary' regularization.

    Matrix<Real> regSmall;
    regSmall.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )        regSmall(i) =  xRegSmall;
        else if( i < n+m ) regSmall(i) = -yRegSmall;
        else               regSmall(i) = -zRegSmall;
    }

    Matrix<Real> regLarge;
    regLarge.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )        regLarge(i) =  xRegLarge;
        else if( i < n+m ) regLarge(i) = -yRegLarge;
        else               regLarge(i) = -zRegLarge;
    }

    // Initialize the static portion of the KKT system
    // ===============================================
    SparseMatrix<Real> JStatic;
    StaticKKT
    ( problem.A, problem.G, Sqrt(xRegSmall), Sqrt(yRegSmall), Sqrt(zRegSmall),
      JStatic, false );
    JStatic.FreezeSparsity();

    SparseLDLFactorization<Real> sparseLDLFact;
    Initialize
    ( problem, solution, JStatic, regLarge,
      sparseLDLFact,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift, ctrl.solveCtrl );

    Int numIts = 0;
    Real dimacsError = 1, dimacsErrorOld = 1;
    Real infeasError = 1;
    Matrix<Real> dInner;
    SparseMatrix<Real> J, JOrig;
    Matrix<Real> d, w;

    auto increaseRegularization = [&]() {
        if( twoStage )
        {
            //if( ctrl.print )
            //    Output("Falling back to single-stage strategy");
            //twoStage = false;
            if( ctrl.print )
                Output
                ("Increasing regularization by a factor of ",
             ctrl.regIncreaseFactor);

            UpdateDiagonal( JStatic, ctrl.regIncreaseFactor-Real(1), regSmall );
            regSmall *= ctrl.regIncreaseFactor;
            xRegSmall *= ctrl.regIncreaseFactor;
            yRegSmall *= ctrl.regIncreaseFactor;
            zRegSmall *= ctrl.regIncreaseFactor;

            regLarge *= ctrl.regIncreaseFactor;
            xRegLarge *= ctrl.regIncreaseFactor;
            yRegLarge *= ctrl.regIncreaseFactor;
            zRegLarge *= ctrl.regIncreaseFactor;
        }
        else
        {
            if( ctrl.print )
                Output
                ("Increasing regularization by a factor of ",
                 ctrl.regIncreaseFactor);
            regLarge *= ctrl.regIncreaseFactor;
            xRegLarge *= ctrl.regIncreaseFactor;
            yRegLarge *= ctrl.regIncreaseFactor;
            zRegLarge *= ctrl.regIncreaseFactor;
        }
    };

    auto attemptToFactor = [&]( const Real& wDynamicRange ) {
        EL_DEBUG_ONLY(auto callStack = CopyCallStack())
        try
        {
            // It seems that straight-forward equilibration can prevent the
            // iterative solver from converging.
            //
            // TODO(poulson): Determine if equilibration is safe on pilot87
            // if twoStage=false. (Yes, but adds a few more iterations.)
            //
            if( twoStage )
            {
                Ones( dInner, J.Height(), 1 );
            }
            else
            {
                if( wDynamicRange >= ctrl.ruizEquilTol )
                    SymmetricRuizEquil
                    ( J, dInner, ctrl.ruizMaxIter, ctrl.print );
                else if( wDynamicRange >= ctrl.diagEquilTol )
                    SymmetricDiagonalEquil( J, dInner, ctrl.print );
                else
                    Ones( dInner, J.Height(), 1 );
            }

            if( numIts == 0 && ctrl.primalInit && ctrl.dualInit )
            {
                const bool hermitian = true;
                const BisectCtrl bisectCtrl;
                sparseLDLFact.Initialize( J, hermitian, bisectCtrl );
            }
            else
            {
                sparseLDLFact.ChangeNonzeroValues( J );
            }
            sparseLDLFact.Factor();
        }
        catch(...)
        {
            if( dimacsError > ctrl.minTol )
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
            EL_DEBUG_ONLY(SetCallStack(callStack))
            return false;
        }
        return true;
      };

    auto attemptToSolve = [&]( Matrix<Real>& rhs )
      {
        if( twoStage )
        {
            EL_DEBUG_ONLY(auto callStack = CopyCallStack())
            RegSolveInfo<Real> solveInfo;
            auto origRhs( rhs );
            try
            {
                solveInfo =
                  reg_ldl::SolveAfter
                  ( JOrig, regLarge, dInner, sparseLDLFact, rhs,
                    ctrl.solveCtrl );
            }
            catch( const std::exception& except )
            {
                if( ctrl.print )
                    Output
                    ("reg_ldl::SolveAfter failed with error: ",except.what());
                EL_DEBUG_ONLY(SetCallStack(callStack))
                rhs = origRhs;
            }
            if( solveInfo.metRequestedTol )
            {
                return true;
            }
            else
            {
                if( dimacsError > ctrl.minTol )
                    Output
                    ("Could not solve with current two-stage regularization");
                return false;
            }
        }
        else
        {
            EL_DEBUG_ONLY(auto callStack = CopyCallStack())
            RegSolveInfo<Real> solveInfo;
            try
            {
                solveInfo =
                  reg_ldl::RegularizedSolveAfter
                  ( JOrig, regLarge, dInner, sparseLDLFact, rhs,
                    ctrl.solveCtrl.relTol,
                    ctrl.solveCtrl.maxRefineIts,
                    ctrl.solveCtrl.progress );
            }
            catch( const std::exception& except )
            {
                if( ctrl.print )
                    Output
                    ("WARNING: solveInfo failed with error: ",except.what());
                EL_DEBUG_ONLY(SetCallStack(callStack))
            }
            if( solveInfo.metRequestedTol )
            {
                return true;
            }
            else
            {
                if( dimacsError > ctrl.minTol )
                    Output
                    ("Could not solve with current one-stage regularization");
                return false;
            }
        }
      };

    AffineLPResidual<Matrix<Real>> residual, error;
    AffineLPSolution<Matrix<Real>> affineCorrection, correction;

    // TODO(poulson): Propagate this through more routines?
    const bool dynamicallyRescale = false;

    const Int indent = PushIndent();

    // We will monotonically drive down the barrier parameter
    // (with the exception of the initial value).
    Real muMin = 1.e6;

    Real mu = 0.1, muOld = 0.1;
    for( ; numIts<=ctrl.maxIts; ++numIts, muOld=mu, dimacsErrorOld=dimacsError )
    {
        if( dynamicallyRescale )
        {
            const Real sScaleNew = MaxNorm(solution.s);
            if( sScaleNew > Real(10) )
            {
                equilibration.sScale *= sScaleNew;
                if( ctrl.print )
                    Output
                    ("s /= ",sScaleNew," (total is ",equilibration.sScale,")");
                problem.b *= Real(1)/sScaleNew;
                problem.h *= Real(1)/sScaleNew;
                solution.x *= Real(1)/sScaleNew;
                solution.s *= Real(1)/sScaleNew;
            }

            // Rescale || c ||_max to roughly one.
            const Real zScaleNew = MaxNorm(solution.z);
            if( zScaleNew > Real(10) )
            {
                equilibration.zScale *= zScaleNew;
                if( ctrl.print )
                    Output
                    ("z /= ",zScaleNew," (total is ",equilibration.zScale,")");
                problem.c *= Real(1)/zScaleNew;
                solution.y *= Real(1)/zScaleNew;
                solution.z *= Real(1)/zScaleNew;
            }
        }

        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = pos_orth::NumOutside( solution.s );
        const Int zNumNonPos = pos_orth::NumOutside( solution.z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the scaling point
        // =========================
        pos_orth::NesterovTodd( solution.s, solution.z, w );
        Real lowerRatio=1; Real upperRatio=1;
        for( Int i=0; i<n; ++i )
        {
            if( solution.s(i) > solution.z(i) )
                upperRatio = Max( upperRatio, solution.s(i)/solution.z(i) );
            if( solution.z(i) > solution.s(i) )
                lowerRatio = Max( lowerRatio, solution.z(i)/solution.s(i) );
        }
        const Real wDynamicRange = lowerRatio*upperRatio;
        if( ctrl.print )
        {
            Output
            ("wLowerRatio=",lowerRatio,", wUpperRatio=",upperRatio,
             ", wDynamicRange=",wDynamicRange);
        }
        const Real compRatio =
          pos_orth::ComplementRatio( solution.s, solution.z );
        if( ctrl.print )
        {
            Output("Complement ratio: ",compRatio);
        }

        // Apply reg' sparingly to the KKT system's bottom-right corner
        auto sMod( solution.s );
        auto zMod( solution.z );
        /*
        const Real maxRatio = Pow(limits::Epsilon<Real>(),Real(-0.35));
        for( Int i=0; i<n; ++i )
        {
            if( solution.z(i)/solution.s(i) > maxRatio )
                zMod(i) = maxRatio*solution.s(i);
        }
        */

        // Check for convergence
        // =====================

        // Carefully compute a relative duality gap.
        const Real dualProd = Dot( solution.s, solution.z );
        const Real primObj = PrimalObjective<Real>( problem, solution );
        const Real dualObj = DualObjective<Real>( problem, solution );
        const Real relGap = RelativeDualityGap( primObj, dualObj, dualProd );

        // || A x - b - gamma_y y ||_2 / (1 + || b ||_2)
        // ---------------------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Multiply
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
        Axpy( -yRegSmall, solution.y, residual.primalEquality );
        if( !twoStage )
            Axpy( -yRegLarge, solution.y, residual.primalEquality );
        const Real primalInfeasNrm2 = Nrm2( residual.primalEquality );
        const Real primalInfeasNrm2Rel = primalInfeasNrm2 / (1+bNrm2);
        // || c + A^T y + G^T z + gamma_x x ||_2 / (1 + || c ||_2)
        // -------------------------------------------------------
        residual.dualEquality = problem.c;
        Multiply
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        Multiply
        ( TRANSPOSE, Real(1), problem.G, solution.z,
          Real(1), residual.dualEquality );
        Axpy( xRegSmall, solution.x, residual.dualEquality );
        if( !twoStage )
            Axpy( xRegLarge, solution.x, residual.dualEquality );
        const Real dualInfeasNrm2 = Nrm2( residual.dualEquality );
        const Real dualInfeasNrm2Rel = dualInfeasNrm2 / (1+cNrm2);
        // || G x + s - h - gamma_z z ||_2 / (1 + || h ||_2)
        // -------------------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Multiply
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        Axpy( -zRegSmall, solution.z, residual.primalConic );
        if( !twoStage )
            Axpy( -zRegLarge, solution.z, residual.primalConic );
        residual.primalConic += solution.s;
        const Real conicInfeasNrm2 = Nrm2( residual.primalConic );
        const Real conicInfeasNrm2Rel = conicInfeasNrm2 / (1+hNrm2);

        // Now check the pieces
        // --------------------
        const Real equalityError = Max(primalInfeasNrm2Rel,dualInfeasNrm2Rel);
        infeasError = Max(conicInfeasNrm2Rel,equalityError);
        dimacsError = Max(relGap,infeasError);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( solution.x );
            const Real yNrm2 = Nrm2( solution.y );
            const Real zNrm2 = Nrm2( solution.z );
            const Real sNrm2 = Nrm2( solution.s );
            Output
            ("iter ",numIts,":\n",Indent(),
             "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
             "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
             "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
             "  ||  s  ||_2 = ",sNrm2,"\n",Indent(),
             "  || primalInfeas ||_2 / (1 + || b ||_2) = ",
             primalInfeasNrm2Rel,"\n",Indent(),
             "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",
             dualInfeasNrm2Rel,"\n",Indent(),
             "  || conicInfeas  ||_2 / (1 + || h ||_2) = ",
             conicInfeasNrm2Rel,"\n",Indent(),
             "  scaled primal = ",primObj,"\n",Indent(),
             "  scaled dual   = ",dualObj,"\n",Indent(),
             "  scaled relative duality gap = ",relGap,"\n",Indent(),
             "  scaled dimacs error = ",dimacsError);

            // TODO(poulson): Move this into a subroutine.
            AffineLPSolution<Matrix<Real>> origSolution;
            origSolution.x = solution.x;
            origSolution.s = solution.s;
            origSolution.y = solution.y;
            origSolution.z = solution.z;
            origSolution.x *= equilibration.sScale;
            origSolution.s *= equilibration.sScale;
            origSolution.y *= equilibration.zScale;
            origSolution.z *= equilibration.zScale;
            DiagonalSolve
            ( LEFT, NORMAL, equilibration.colScale,  origSolution.x );
            DiagonalSolve
            ( LEFT, NORMAL, equilibration.rowScaleA, origSolution.y );
            DiagonalSolve
            ( LEFT, NORMAL, equilibration.rowScaleG, origSolution.z );

            const Real dualProdOrig = Dot( origSolution.s, origSolution.z );
            const Real primObjOrig =
              PrimalObjective<Real>( origProblem, origSolution );
            const Real dualObjOrig =
              DualObjective<Real>( origProblem, origSolution );
            const Real relGapOrig =
              RelativeDualityGap( primObjOrig, dualObjOrig, dualProdOrig );

            Output("  s^T z = ",dualProdOrig);
            Output("  primal: ",primObjOrig);
            Output("  dual: ",dualObjOrig);
            Output("  relative gap: ",relGapOrig);
        }
        if( dimacsError <= ctrl.targetTol )
            break;
        // Exit if progress has stalled and we are sufficiently accurate.
        if( dimacsError <= ctrl.minTol &&
            dimacsError >= ctrl.minDimacsDecreaseRatio*dimacsErrorOld )
            break;
        if( numIts == ctrl.maxIts && dimacsError > ctrl.minTol )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================
        const bool freezeBarrier = compRatio > 1000 || relGap < equalityError;
        if( freezeBarrier )
        {
            mu = muOld;
            if( ctrl.print )
                Output("Freezing barrier at ",mu);
        }
        else
        {
            mu = Dot(solution.s,solution.z) / k;
            mu = Min( mu, muMin );
            muMin = Min( mu, muMin );
            if( ctrl.print )
                Output("New barrier parameter is ",mu);
        }

        // r_mu := s o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

        // Construct the KKT system
        // ------------------------
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        FinishKKT( m, n, sMod, zMod, JOrig );
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        J = JOrig;
        J.FreezeSparsity();
        UpdateDiagonal( J, Real(1), regLarge );

        // Solve for the direction
        // -----------------------
        if( !attemptToFactor(wDynamicRange) )
        {
            increaseRegularization();
            continue;
        }
        if( !attemptToSolve(d) )
        {
            increaseRegularization();
            continue;
        }
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, solution.z,
          affineCorrection.x,
          affineCorrection.y,
          affineCorrection.z,
          affineCorrection.s );

        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Multiply
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            Axpy( -yRegSmall, affineCorrection.y, error.primalEquality );
            if( !twoStage )
                Axpy( -yRegLarge, affineCorrection.y, error.primalEquality );
            const Real dxErrorNrm2 = Nrm2( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Multiply
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Multiply
            ( TRANSPOSE, Real(1), problem.G, affineCorrection.z,
              Real(1), error.dualEquality );
            Axpy( xRegSmall, affineCorrection.x, error.dualEquality );
            if( !twoStage )
                Axpy( xRegLarge, affineCorrection.x, error.dualEquality );
            const Real dyErrorNrm2 = Nrm2( error.dualEquality );

            error.primalConic = residual.primalConic;
            Multiply
            ( NORMAL, Real(1), problem.G, affineCorrection.x,
              Real(1), error.primalConic );
            error.primalConic += affineCorrection.s;
            Axpy( -zRegSmall, affineCorrection.z, error.primalConic );
            if( !twoStage )
                Axpy( -zRegLarge, affineCorrection.z, error.primalConic );
            const Real dzErrorNrm2 = Nrm2( error.primalConic );

            // TODO(poulson): error.dualConic
            // TODO(poulson): Also compute and print the residuals with
            // regularization.

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+primalInfeasNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+dualInfeasNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+conicInfeasNrm2));
        }

        // Prevent the update from decreasing values of s or z below sqrt(eps)
        // ===================================================================
        Real sMin = MaxNorm( solution.s );
        Real zMin = MaxNorm( solution.z );
        for( Int i=0; i<n; ++i )
        {
            sMin = Min( sMin, solution.s(i) );
            zMin = Min( zMin, solution.z(i) );
            /*
            if( solution.s(i) < Pow(limits::Epsilon<Real>(),Real(0.75)) )
            {
                if( affineCorrection.s(i) < Real(0) )
                {
                    Output("Clipping affineCorrection.s(",i,")=",affineCorrection.s(i)," to zero since s(",i,")=",solution.s(i));
                    affineCorrection.s(i) = 0;
                }
            }
            */
            if( solution.z(i) < Pow(limits::Epsilon<Real>(),Real(0.75)) )
            {
                if( affineCorrection.z(i) < Real(0) )
                {
                    Output("Clipping affineCorrection.z(",i,")=",affineCorrection.z(i)," to zero since z(",i,")=",solution.z(i));
                    affineCorrection.z(i) = 0;
                }
            }
        }
        Output("sMin = ",sMin,", zMin = ",zMin);

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri =
          pos_orth::MaxStep( solution.s, affineCorrection.s, Real(1) );
        Real alphaAffDual =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and ds are used as temporaries
        correction.s = solution.s;
        correction.z = solution.z;
        Axpy( alphaAffPri,  affineCorrection.s, correction.s );
        Axpy( alphaAffDual, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.s,correction.z) / degree;
        if( ctrl.print )
            Output("muAff = ",muAff,", mu = ",mu);
        // TODO(poulson): Avoid unnecessary extra solve when freezing barrier
        Real sigma;
        if( freezeBarrier )
        {
            sigma = 1;
            if( ctrl.print )
                Output("Freezing sigma at one");
        }
        else
        {
            sigma = ctrl.centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
            if( ctrl.print )
                Output("Computed sigma=",sigma);
        }
        if( ctrl.print )
            Output("sigma=",sigma);

        // TODO(poulson): Experiment with several centrality choices?

        // Solve for the combined direction
        // ================================
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu := dsAff o dzAff
            // ---------------------
            // NOTE: dz is used as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.s, correction.z );
            residual.dualConic += correction.z;
        }

        // Construct the new KKT RHS
        // -------------------------
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        // Solve for the proposed step
        // ---------------------------
        if( !attemptToSolve(d) )
        {
            increaseRegularization();
            continue;
        }
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, solution.z,
          correction.x, correction.y, correction.z, correction.s );

        for( Int i=0; i<n; ++i )
        {
            /*
            if( solution.s(i) < Pow(limits::Epsilon<Real>(),Real(0.75)) )
            {
                if( correction.s(i) < Real(0) )
                {
                    Output("Clipping correction.s(",i,")=",correction.s(i)," to zero since s(",i,")=",solution.s(i));
                    correction.s(i) = 0;
                }
            }
            */
            if( solution.z(i) < Pow(limits::Epsilon<Real>(),Real(0.75)) )
            {
                if( correction.z(i) < Real(0) )
                {
                    Output("Clipping correction.z(",i,")=",correction.z(i)," to zero since z(",i,")=",solution.z(i));
                    correction.z(i) = 0;
                }
            }
        }

        // Update the current estimates
        // ============================
        Real alphaPri =
          pos_orth::MaxStep( solution.s, correction.s, 1/ctrl.maxStepRatio );
        Real alphaDual =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaPri,  correction.s, solution.s );
        Axpy( alphaDual, correction.y, solution.y );
        Axpy( alphaDual, correction.z, solution.z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( dimacsError <= ctrl.minTol )
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
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.outerEquil )
    {
        AffineLPProblem<SparseMatrix<Real>,Matrix<Real>> equilibratedProblem;
        AffineLPSolution<Matrix<Real>> equilibratedSolution;
        SparseAffineLPEquilibration<Real> equilibration;
        Equilibrate
        ( problem, solution,
          equilibratedProblem, equilibratedSolution,
          equilibration, ctrl );
        EquilibratedMehrotra
        ( problem, equilibration,
          equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        SparseAffineLPEquilibration<Real> equilibration;
        Ones( equilibration.rowScaleA, problem.A.Height(), 1 );
        Ones( equilibration.rowScaleG, problem.G.Height(), 1 );
        Ones( equilibration.colScale, problem.A.Width(), 1 );
        equilibration.sScale = 1;
        equilibration.zScale = 1;
        auto equilibratedProblem = problem;
        auto equilibratedSolution = solution;
        EquilibratedMehrotra
        ( problem, equilibration,
          equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
}

template<typename Real>
void EquilibratedMehrotra
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Int degree = k;
    const Grid& grid = problem.A.Grid();
    const int commRank = grid.Rank();
    Timer timer;

    const Real bNrm2 = Nrm2( problem.b );
    const Real cNrm2 = Nrm2( problem.c );
    const Real hNrm2 = Nrm2( problem.h );
    const Real twoNormEstA = TwoNormEstimate( problem.A, ctrl.basisSize );
    const Real twoNormEstG = TwoNormEstimate( problem.G, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstA + twoNormEstG + 1;
    if( ctrl.print )
    {
        const double imbalanceA = problem.A.Imbalance();
        const double imbalanceG = problem.G.Imbalance();
        if( commRank == 0 )
        {
            Output("|| c ||_2 = ",cNrm2);
            Output("|| A ||_2 estimate: ",twoNormEstA);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| G ||_2 estimate: ",twoNormEstG);
            Output("|| h ||_2 = ",hNrm2);
            Output("Imbalance factor of A: ",imbalanceA);
            Output("Imbalance factor of G: ",imbalanceG);
        }
    }

    DistMultiVec<Real> regLarge(grid);
    regLarge.Resize( n+m+k, 1 );
    for( Int iLoc=0; iLoc<regLarge.LocalHeight(); ++iLoc )
    {
        const Int i = regLarge.GlobalRow(iLoc);
        if( i < n )
          regLarge.SetLocal( iLoc, 0, ctrl.xRegLarge );
        else if( i < n+m )
          regLarge.SetLocal( iLoc, 0, -ctrl.yRegLarge );
        else
          regLarge.SetLocal( iLoc, 0, -ctrl.zRegLarge );
    }
    regLarge *= origTwoNormEst;

    // Construct the static part of the KKT system
    // ===========================================
    DistSparseMatrix<Real> JStatic(grid);
    StaticKKT
    ( problem.A, problem.G,
      Sqrt(ctrl.xRegSmall), Sqrt(ctrl.yRegSmall), Sqrt(ctrl.zRegSmall),
      JStatic, false );
    JStatic.FreezeSparsity();
    JStatic.InitializeMultMeta();
    if( ctrl.print )
    {
        const double imbalanceJ = JStatic.Imbalance();
        if( commRank == 0 )
            Output("Imbalance factor of J: ",imbalanceJ);
    }

    if( commRank == 0 && ctrl.time )
        timer.Start();
    DistSparseLDLFactorization<Real> sparseLDLFact;
    Initialize
    ( problem, solution, JStatic, regLarge,
      sparseLDLFact,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift, ctrl.solveCtrl );
    if( commRank == 0 && ctrl.time )
        Output("Init: ",timer.Stop()," secs");

    Int numIts = 0;
    Real dimacsError = 1, dimacsErrorOld = 1;
    Real infeasError = 1;
    DistSparseMatrix<Real> J(grid), JOrig(grid);
    DistMultiVec<Real> d(grid), w(grid), dInner(grid);
    auto attemptToFactor = [&]( const Real& wMaxNorm )
      {
        EL_DEBUG_ONLY(auto callStack = CopyCallStack())
        try
        {
            /*
            if( commRank == 0 && ctrl.time )
                timer.Start();
            if( wMaxNorm >= ctrl.ruizEquilTol )
                SymmetricRuizEquil( J, dInner, ctrl.ruizMaxIter, ctrl.print );
            else if( wMaxNorm >= ctrl.diagEquilTol )
                SymmetricDiagonalEquil( J, dInner, ctrl.print );
            else
                Ones( dInner, J.Height(), 1 );
            if( commRank == 0 && ctrl.time )
                Output("Equilibration: ",timer.Stop()," secs");
            */
            Ones( dInner, J.Height(), 1 );

            if( numIts == 0 && ctrl.primalInit && ctrl.dualInit )
            {
                const bool hermitian = true;
                const BisectCtrl bisectCtrl;
                sparseLDLFact.Initialize( J, hermitian, bisectCtrl );
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
        }
        catch(...)
        {
            if( dimacsError > ctrl.minTol )
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
            EL_DEBUG_ONLY(SetCallStack(callStack))
            return false;
        }
        return true;
      };
    auto attemptToSolve = [&]( DistMultiVec<Real>& rhs )
      {
        if( commRank == 0 && ctrl.time )
            timer.Start();
        if( ctrl.twoStage )
        {
            auto solveInfo =
              reg_ldl::SolveAfter
              ( JOrig, regLarge, dInner, sparseLDLFact, rhs, ctrl.solveCtrl );
            if( solveInfo.metRequestedTol )
            {
                if( commRank == 0 && ctrl.time )
                    Output("Affine: ",timer.Stop()," secs");
                return true;
            }
            else
            {
                if( commRank == 0 )
                    Output("WARNING: Could not resolve regularization");
                // TODO(poulson): twoStage = false
            }
        }
        auto solveInfo =
          reg_ldl::RegularizedSolveAfter
          ( JOrig, regLarge, dInner, sparseLDLFact, rhs,
            ctrl.solveCtrl.relTol,
            ctrl.solveCtrl.maxRefineIts,
            ctrl.solveCtrl.progress );
        if( commRank == 0 && ctrl.time )
            Output("Affine: ",timer.Stop()," secs");
        if( solveInfo.metRequestedTol )
        {
            return true;
        }
        else
        {
            if( dimacsError > ctrl.minTol )
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
            return false;
        }
      };

    AffineLPResidual<DistMultiVec<Real>> residual, error;
    AffineLPSolution<DistMultiVec<Real>> affineCorrection, correction;

    ForceSimpleAlignments( residual, grid );
    ForceSimpleAlignments( error, grid );
    ForceSimpleAlignments( affineCorrection, grid );
    ForceSimpleAlignments( correction, grid );

    const Int indent = PushIndent();
    for( ; numIts<=ctrl.maxIts; ++numIts, dimacsErrorOld=dimacsError )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = pos_orth::NumOutside( solution.s );
        const Int zNumNonPos = pos_orth::NumOutside( solution.z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure and scaling point
        // =============================================
        const Real mu = Dot(solution.s,solution.z) / k;
        pos_orth::NesterovTodd( solution.s, solution.z, w );
        const Real wMaxNorm = MaxNorm( w );

        // Check for convergence
        // =====================
        const Real dualProd = Dot( solution.s, solution.z );
        const Real primObj = PrimalObjective<Real>( problem, solution );
        const Real dualObj = DualObjective<Real>( problem, solution );
        const Real relGap = RelativeDualityGap( primObj, dualObj, dualProd );
        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Multiply
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
        Axpy( -ctrl.yRegSmall, solution.y, residual.primalEquality );
        const Real primalInfeasNrm2 = Nrm2( residual.primalEquality );
        const Real primalInfeasNrm2Rel = primalInfeasNrm2 / (1+bNrm2);

        // || c + A^T y + G^T z ||_2 / (1 + || c ||_2)
        // -------------------------------------------
        residual.dualEquality = problem.c;
        Multiply
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        Multiply
        ( TRANSPOSE, Real(1), problem.G, solution.z,
          Real(1), residual.dualEquality );
        Axpy( ctrl.xRegSmall, solution.x, residual.dualEquality );
        const Real dualInfeasNrm2 = Nrm2( residual.dualEquality );
        const Real dualInfeasNrm2Rel = dualInfeasNrm2 / (1+cNrm2);

        // || G x + s - h ||_2 / (1 + || h ||_2)
        // -------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Multiply
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        Axpy( -ctrl.zRegSmall, solution.z, residual.primalConic );
        residual.primalConic += solution.s;
        const Real conicInfeasNrm2 = Nrm2( residual.primalConic );
        const Real conicInfeasNrm2Rel = conicInfeasNrm2 / (1+hNrm2);

        // Now check the pieces
        // --------------------
        const Real equalityError = Max(primalInfeasNrm2Rel,dualInfeasNrm2Rel);
        infeasError = Max(equalityError,conicInfeasNrm2Rel);
        dimacsError = Max(infeasError,relGap);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( solution.x );
            const Real yNrm2 = Nrm2( solution.y );
            const Real zNrm2 = Nrm2( solution.z );
            const Real sNrm2 = Nrm2( solution.s );
            if( commRank == 0 )
                Output
                ("iter ",numIts,":\n",Indent(),
                 "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
                 "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
                 "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
                 "  ||  s  ||_2 = ",sNrm2,"\n",Indent(),
                 "  || primalInfeas ||_2 / (1 + || b ||_2) = ",
                 primalInfeasNrm2Rel,"\n",Indent(),
                 "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",
                 dualInfeasNrm2Rel,"\n",Indent(),
                 "  || conicInfeas  ||_2 / (1 + || h ||_2) = ",
                 conicInfeasNrm2Rel,"\n",Indent(),
                 "  scaled primal = ",primObj,"\n",Indent(),
                 "  scaled dual   = ",dualObj,"\n",Indent(),
                 "  scaled relative gap = ",relGap);
        }
        if( dimacsError <= ctrl.targetTol )
            break;
        // Exit if progress has stalled and we are sufficiently accurate.
        if( dimacsError <= ctrl.minTol &&
            dimacsError >= ctrl.minDimacsDecreaseRatio*dimacsErrorOld )
            break;
        if( numIts == ctrl.maxIts && dimacsError > ctrl.minTol )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := s o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

        // Construct the KKT system
        // ------------------------
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        JOrig.LockedDistGraph().multMeta = JStatic.LockedDistGraph().multMeta;
        FinishKKT( m, n, solution.s, solution.z, JOrig );
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        J = JOrig;
        J.FreezeSparsity();
        J.LockedDistGraph().multMeta = JStatic.LockedDistGraph().multMeta;
        UpdateDiagonal( J, Real(1), regLarge );

        // Solve for the direction
        // -----------------------
        if( !attemptToFactor(wMaxNorm) )
            break;
        if( !attemptToSolve(d) )
            break;
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, solution.z,
          affineCorrection.x,
          affineCorrection.y,
          affineCorrection.z,
          affineCorrection.s );

        if( ctrl.checkResiduals && ctrl.print )
        {
            error.primalEquality = residual.primalEquality;
            Multiply
            ( NORMAL, Real(1), problem.A, affineCorrection.x,
              Real(1), error.primalEquality );
            const Real dxErrorNrm2 = Nrm2( error.primalEquality );

            error.dualEquality = residual.dualEquality;
            Multiply
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Multiply
            ( TRANSPOSE, Real(1), problem.G, affineCorrection.z,
              Real(1), error.dualEquality );
            const Real dyErrorNrm2 = Nrm2( error.dualEquality );

            error.primalConic = residual.primalConic;
            Multiply
            ( NORMAL, Real(1), problem.G, affineCorrection.x,
              Real(1), error.primalConic );
            error.primalConic += affineCorrection.s;
            const Real dzErrorNrm2 = Nrm2( error.primalConic );

            // TODO(poulson): error.dualConic
            // TODO(poulson): Also compute and print the residuals with
            // regularization

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
                 dxErrorNrm2/(1+primalInfeasNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+dualInfeasNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+conicInfeasNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri =
          pos_orth::MaxStep( solution.s, affineCorrection.s, Real(1) );
        Real alphaAffDual =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: correction.z and correction.s are used as temporaries
        correction.s = solution.s;
        correction.z = solution.z;
        Axpy( alphaAffPri,  affineCorrection.s, correction.s );
        Axpy( alphaAffDual, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.s,correction.z) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma =
          ctrl.centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dsAff o dzAff
            // ---------------------
            // NOTE: correction.z is used as a temporary
            correction.z = affineCorrection.z;
            DiagonalScale( LEFT, NORMAL, affineCorrection.s, correction.z );
            residual.dualConic += correction.z;
        }

        // Construct the new KKT RHS
        // -------------------------
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        // Solve for the direction
        // -----------------------
        if( !attemptToSolve(d) )
            break;
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, solution.z,
          correction.x, correction.y, correction.z, correction.s );

        // Update the current estimates
        // ============================
        Real alphaPri =
          pos_orth::MaxStep( solution.s, correction.s, 1/ctrl.maxStepRatio );
        Real alphaDual =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaPri,  correction.s, solution.s );
        Axpy( alphaDual, correction.y, solution.y );
        Axpy( alphaDual, correction.z, solution.z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( dimacsError <= ctrl.minTol )
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
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const MehrotraCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    if( ctrl.outerEquil )
    {
        const Grid& grid = problem.A.Grid();
        AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>
          equilibratedProblem;
        AffineLPSolution<DistMultiVec<Real>> equilibratedSolution;
        DistSparseAffineLPEquilibration<Real> equilibration;

        ForceSimpleAlignments( equilibratedSolution, grid );
        ForceSimpleAlignments( equilibratedProblem, grid );

        Equilibrate
        ( problem, solution,
          equilibratedProblem, equilibratedSolution,
          equilibration, ctrl );
        EquilibratedMehrotra( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedMehrotra( problem, solution, ctrl );
    }
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
          AffineLPSolution<Matrix<Real>>& solution, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
          AffineLPSolution<DistMatrix<Real>>& solution, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
          AffineLPSolution<Matrix<Real>>& solution, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
          AffineLPSolution<DistMultiVec<Real>>& solution, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace affine
} // namespace lp
} // namespace El
