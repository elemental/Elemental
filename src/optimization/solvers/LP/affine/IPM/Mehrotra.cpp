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

namespace lp {
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
        Output("|| A ||_1 = ",ANrm1);
        Output("|| G ||_1 = ",GNrm1);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
        Output("|| h ||_2 = ",hNrm2);
    }

    Initialize
    ( problem, solution,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift );

    Real relError = 1;
    Matrix<Real> J, d;
    Matrix<Real> dSub;
    Permutation p;
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
    auto attemptToSolve = [&]( Matrix<Real>& rhs )
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

    AffineLPSolution<Matrix<Real>> affineCorrection, correction;
    AffineLPResidual<Matrix<Real>> residual, error;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
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
        const Real mu = Dot(solution.s,solution.z) / k;

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y) -
          Dot(problem.h,solution.z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
        const Real rbNrm2 = Nrm2( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        Gemv
        ( TRANSPOSE, Real(1), problem.G, solution.z,
          Real(1), residual.dualEquality );
        const Real rcNrm2 = Nrm2( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Gemv
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        residual.primalConic += solution.s;
        const Real rhNrm2 = Nrm2( residual.primalConic );
        const Real rhConv = rhNrm2 / (1+hNrm2);
        // Now check the pieces
        // --------------------
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
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
             "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
             "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
             "  || r_h ||_2 = ",rhNrm2,"\n",Indent(),
             "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
             "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
             "  || r_h ||_2 / (1 + || h ||_2) = ",rhConv,"\n",Indent(),
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
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+rhNrm2));
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
        residual.primalEquality *= 1-sigma;
        residual.primalConic *= 1-sigma;
        residual.dualEquality *= 1-sigma;
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
            Output("|| A ||_1 = ",ANrm1);
            Output("|| G ||_1 = ",GNrm1);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
            Output("|| h ||_2 = ",hNrm2);
        }
    }

    Initialize
    ( problem, solution,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift );

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

    AffineLPResidual<DistMatrix<Real>> residual, error;
    AffineLPSolution<DistMatrix<Real>> affineCorrection, correction;
    ForceSimpleAlignments( residual, grid );
    ForceSimpleAlignments( error, grid );
    ForceSimpleAlignments( affineCorrection, grid );
    ForceSimpleAlignments( correction, grid );

    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
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
        const Real mu = Dot(solution.s,solution.z) / k;

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y) -
          Dot(problem.h,solution.z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
        const Real rbNrm2 = Nrm2( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        residual.dualEquality = problem.c;
        Gemv
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        Gemv
        ( TRANSPOSE, Real(1), problem.G, solution.z,
          Real(1), residual.dualEquality );
        const Real rcNrm2 = Nrm2( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Gemv
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        residual.primalConic += solution.s;
        const Real rhNrm2 = Nrm2( residual.primalConic );
        const Real rhConv = rhNrm2 / (1+hNrm2);
        // Now check the pieces
        // --------------------
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
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
                 "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
                 "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
                 "  || r_h ||_2 = ",rhNrm2,"\n",Indent(),
                 "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
                 "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
                 "  || r_h ||_2 / (1 + || h ||_2) = ",rhConv,"\n",Indent(),
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
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rhNrm2));
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
        residual.primalEquality *= 1-sigma;
        residual.primalConic *= 1-sigma;
        residual.dualEquality *= 1-sigma;
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
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
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

    Matrix<Real> regTmp;
    regTmp.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )        regTmp(i) =  ctrl.reg0Tmp*ctrl.reg0Tmp;
        else if( i < n+m ) regTmp(i) = -ctrl.reg1Tmp*ctrl.reg1Tmp;
        else               regTmp(i) = -ctrl.reg2Tmp*ctrl.reg2Tmp;
    }
    regTmp *= origTwoNormEst;

    // Initialize the static portion of the KKT system
    // ===============================================
    SparseMatrix<Real> JStatic;
    StaticKKT
    ( problem.A, problem.G, ctrl.reg0Perm, ctrl.reg1Perm, ctrl.reg2Perm,
      JStatic, false );
    JStatic.FreezeSparsity();

    SparseLDLFactorization<Real> sparseLDLFact;
    Initialize
    ( problem, solution, JStatic, regTmp,
      sparseLDLFact,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift, ctrl.solveCtrl );

    Int numIts = 0;
    Real relError = 1;
    Matrix<Real> dInner;
    SparseMatrix<Real> J, JOrig;
    Matrix<Real> d, w;
    auto attemptToFactor = [&]( const Real& wMaxNorm )
      {
        try
        {
            if( wMaxNorm >= ctrl.ruizEquilTol )
                SymmetricRuizEquil( J, dInner, ctrl.ruizMaxIter, ctrl.print );
            else if( wMaxNorm >= ctrl.diagEquilTol )
                SymmetricDiagonalEquil( J, dInner, ctrl.print );
            else
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
            sparseLDLFact.Factor();
        }
        catch(...)
        {
            if( relError > ctrl.minTol )
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
            return false;
        }
        return true;
      };
    auto attemptToSolve = [&]( Matrix<Real>& rhs )
      {
        try
        {
            if( ctrl.resolveReg )
                reg_ldl::SolveAfter
                ( JOrig, regTmp, dInner, sparseLDLFact, rhs, ctrl.solveCtrl );
            else
                reg_ldl::RegularizedSolveAfter
                ( JOrig, regTmp, dInner, sparseLDLFact, rhs,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress );
        }
        catch(...)
        {
            if( relError > ctrl.minTol )
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
            return false;
        }
        return true;
      };

    AffineLPResidual<Matrix<Real>> residual, error;
    AffineLPSolution<Matrix<Real>> affineCorrection, correction;

    const Int indent = PushIndent();
    for( ; numIts<=ctrl.maxIts; ++numIts )
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
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y) -
          Dot(problem.h,solution.z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Multiply
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
        const Real rbNrm2 = Nrm2( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // TODO(poulson): Document this living *after* the norm computations
        Axpy
        ( -ctrl.reg1Perm*ctrl.reg1Perm, solution.y, residual.primalEquality );
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        residual.dualEquality = problem.c;
        Multiply
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        Multiply
        ( TRANSPOSE, Real(1), problem.G, solution.z,
          Real(1), residual.dualEquality );
        const Real rcNrm2 = Nrm2( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // TODO(poulson): Document this living *after* the norm computations
        Axpy( ctrl.reg0Perm*ctrl.reg0Perm, solution.x, residual.dualEquality );
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Multiply
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        residual.primalConic += solution.s;
        const Real rhNrm2 = Nrm2( residual.primalConic );
        const Real rhConv = rhNrm2 / (1+hNrm2);
        // TODO(poulson): Document this living *after* the norm computations
        Axpy( -ctrl.reg2Perm*ctrl.reg2Perm, solution.z, residual.primalConic );

        // Now check the pieces
        // --------------------
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
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
             "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
             "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
             "  || r_h ||_2 = ",rhNrm2,"\n",Indent(),
             "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
             "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
             "  || r_h ||_2 / (1 + || h ||_2) = ",rhConv,"\n",Indent(),
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

        // r_mu := s o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

        // Construct the KKT system
        // ------------------------
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        FinishKKT( m, n, solution.s, solution.z, JOrig );
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        J = JOrig;
        J.FreezeSparsity();
        UpdateDiagonal( J, Real(1), regTmp );

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
            // regularization.

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+rhNrm2));
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
        residual.primalEquality *= 1-sigma;
        residual.primalConic *= 1-sigma;
        residual.dualEquality *= 1-sigma;
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
        if( ctrl.print )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  correction.x, solution.x );
        Axpy( alphaPri,  correction.s, solution.s );
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
            Output("|| A ||_2 estimate: ",twoNormEstA);
            Output("|| G ||_2 estimate: ",twoNormEstG);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
            Output("|| h ||_2 = ",hNrm2);
            Output("Imbalance factor of A: ",imbalanceA);
            Output("Imbalance factor of G: ",imbalanceG);
        }
    }

    DistMultiVec<Real> regTmp(grid);
    regTmp.Resize( n+m+k, 1 );
    for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
    {
        const Int i = regTmp.GlobalRow(iLoc);
        if( i < n )
          regTmp.SetLocal( iLoc, 0, ctrl.reg0Tmp*ctrl.reg0Tmp );
        else if( i < n+m )
          regTmp.SetLocal( iLoc, 0, -ctrl.reg1Tmp*ctrl.reg1Tmp );
        else
          regTmp.SetLocal( iLoc, 0, -ctrl.reg2Tmp*ctrl.reg2Tmp );
    }
    regTmp *= origTwoNormEst;

    // Construct the static part of the KKT system
    // ===========================================
    DistSparseMatrix<Real> JStatic(grid);
    StaticKKT
    ( problem.A, problem.G, ctrl.reg0Perm, ctrl.reg1Perm, ctrl.reg2Perm,
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
    ( problem, solution, JStatic, regTmp,
      sparseLDLFact,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift, ctrl.solveCtrl );
    if( commRank == 0 && ctrl.time )
        Output("Init: ",timer.Stop()," secs");

    Int numIts = 0;
    Real relError = 1;
    DistSparseMatrix<Real> J(grid), JOrig(grid);
    DistMultiVec<Real> d(grid), w(grid), dInner(grid);
    auto attemptToFactor = [&]( const Real& wMaxNorm )
      {
        try
        {
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
            if( relError > ctrl.minTol )
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
            return false;
        }
        return true;
      };
    auto attemptToSolve = [&]( DistMultiVec<Real>& rhs )
      {
        try
        {
            if( commRank == 0 && ctrl.time )
                timer.Start();
            if( ctrl.resolveReg )
                reg_ldl::SolveAfter
                ( JOrig, regTmp, dInner, sparseLDLFact, rhs, ctrl.solveCtrl );
            else
                reg_ldl::RegularizedSolveAfter
                ( JOrig, regTmp, dInner, sparseLDLFact, rhs,
                  ctrl.solveCtrl.relTol,
                  ctrl.solveCtrl.maxRefineIts,
                  ctrl.solveCtrl.progress );
            if( commRank == 0 && ctrl.time )
                Output("Affine: ",timer.Stop()," secs");
        }
        catch(...)
        {
            if( relError > ctrl.minTol )
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
            return false;
        }
        return true;
      };

    AffineLPResidual<DistMultiVec<Real>> residual, error;
    AffineLPSolution<DistMultiVec<Real>> affineCorrection, correction;

    ForceSimpleAlignments( residual, grid );
    ForceSimpleAlignments( error, grid );
    ForceSimpleAlignments( affineCorrection, grid );
    ForceSimpleAlignments( correction, grid );

    const Int indent = PushIndent();
    for( ; numIts<=ctrl.maxIts; ++numIts )
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
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(problem.c,solution.x);
        const Real dualObj = -Dot(problem.b,solution.y) -
          Dot(problem.h,solution.z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Multiply
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
        const Real rbNrm2 = Nrm2( residual.primalEquality );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // TODO(poulson): Document this living *after* the norm computations
        Axpy
        ( -ctrl.reg1Perm*ctrl.reg1Perm, solution.y, residual.primalEquality );
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        residual.dualEquality = problem.c;
        Multiply
        ( TRANSPOSE, Real(1), problem.A, solution.y,
          Real(1), residual.dualEquality );
        Multiply
        ( TRANSPOSE, Real(1), problem.G, solution.z,
          Real(1), residual.dualEquality );
        const Real rcNrm2 = Nrm2( residual.dualEquality );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // TODO(poulson): Document this living *after* the norm computations
        Axpy( ctrl.reg0Perm*ctrl.reg0Perm, solution.x, residual.dualEquality );
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Multiply
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        residual.primalConic += solution.s;
        const Real rhNrm2 = Nrm2( residual.primalConic );
        const Real rhConv = rhNrm2 / (1+hNrm2);
        // TODO(poulson): Document this living *after* the norm computations
        Axpy( -ctrl.reg2Perm*ctrl.reg2Perm, solution.z, residual.primalConic );

        // Now check the pieces
        // --------------------
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
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
                 "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
                 "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
                 "  || r_h ||_2 = ",rhNrm2,"\n",Indent(),
                 "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
                 "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
                 "  || r_h ||_2 / (1 + || h ||_2) = ",rhConv,"\n",Indent(),
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
        UpdateDiagonal( J, Real(1), regTmp );

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
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rhNrm2));
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
        residual.primalEquality *= 1-sigma;
        residual.primalConic *= 1-sigma;
        residual.dualEquality *= 1-sigma;
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
