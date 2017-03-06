/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_SOLVERS_UTIL_HPP
#define EL_OPTIMIZATION_SOLVERS_UTIL_HPP

namespace El {

namespace KKTSystemNS {
enum KKTSystem {
  FULL_KKT,
  AUGMENTED_KKT,
  NORMAL_KKT
};
}
using namespace KKTSystemNS;

// Infeasible Interior Point Method
// ================================
template<typename Real>
Real StepLengthCentrality
( Real barrier, Real barrierAff, Real primalAffineStep, Real dualAffineStep )
{ return Pow(1-Min(primalAffineStep,dualAffineStep),Real(3)); }

template<typename Real>
Real PredictorCorrectorCentrality
( Real mu, Real muAff, Real alphaAffPri, Real alphaAffDual )
{ return Min(Pow(muAff/mu,Real(3)),Real(1)); }

template<typename Real>
struct IPMCtrl
{
    // Mark whether or not the primal and/or dual variables were
    // user-initialized.
    //
    // For problems with 'direct' cone constraints, i.e., x in K, the primal
    // variable is 'x' and the dual variables are 'y' and 'z'. For problems with
    // 'affine' cone constraints, i.e., (h - G x) in K, the primal variables are
    // 'x' and 's', while the dual variables are again 'y' and 'z'.
    //
    // NOTE: User initialization has not yet been tested, so there might be a
    //       trivial bug.
    bool primalInit = false, dualInit = false;

    // Demand that
    //
    //   max{ || A^T y + G^T z + c ||_2 / (1 + || c ||_2),
    //        || A x - b ||_2 / (1 + || b ||_2),
    //        || G x + s - h ||_2 / (1 + || h ||_2) }
    //
    // is less than this value.
    Real infeasibilityTolLogEps = Real(0.45);

    // Demand that
    //
    //   | primalObjective - dualObjective | /
    //   (max{ |primalObjective|, |dualObjective| } + 1)
    //
    // is less than this value. Unfortunately, some of the LP_data examples
    // (in particular, 80BAU3B) are challenging with symmetric indefinite KKT
    // system solvers and it seems challenging to achieve more than two digits
    // of accuracy in double precision with this metric. We therefore default
    // to only demanding that the relative complementarity gap is nontrivial
    // by default. (Although using larger amounts of regularization yields more
    // than three digits of relative objective gap accuracy.)
    Real relativeObjectiveGapTolLogEps = Real(0.05);

    // Demand that
    //
    //   s^T z / max{ -c^T x, -b^T y - h^T z }
    //
    // is less than this value. Note that this would be equivalent to the
    // relative objective gap if the iterates were feasible, as 'A x = b' and
    // 'G x + s = h' implies that
    //
    //   -b^T y - h^T z = -(A x)^T y - (G x + s)^T z
    //                  = -x^T (A^T y + G^T z) - s^T z,
    //
    // and 'A^T y + G^T z + c = 0' further implies that
    //
    //   -b^T y - h^T z = c^T x - s^T z.
    //
    Real relativeComplementarityGapTolLogEps = Real(0.3);

    // If the minimum tolerance has already been achieved, then exit if the
    // DIMACS error is multiplied by a factor larger than the following on any
    // subsequent iteration.
    Real minDimacsDecreaseRatio = Real(0.99);

    // The maximum number of iterations of the IPM. This should only be
    // activated in pathological circumstances, as even 100 iterations of an
    // IPM is usually excessive.
    Int maxIts = 1000;

    // If a step of length alpha would hit the boundary, instead step a distance
    // of 'maxStepRatio*alpha'.
    Real maxStepRatio = Real(0.99);

    // For configuring how many reductions of the first-order optimality
    // conditions should be performed before solving a linear system.
    // For problems with 'affine' cone constraints, i.e., (h - G x) in K,
    // this should remain as FULL_KKT. For sparse problems with 'direct' cone
    // constraints, i.e., x in K, both AUGMENTED_KKT (use a QSD solver) and
    // NORMAL_KKT (use a Cholesky solver) are also possible. The latter should
    // be avoided when the normal equations are sufficiently denser than the
    // (larger) augmented formulation.
    KKTSystem system = FULL_KKT;

    // Use a damped (level-1) composite Newton method?
    //
    // Cf. Tapia et al., "The Mehrotra Predictor-Corrector Interior-Point
    // Method as a Perturbed Composite Newton Method", July, 1990.
    // http://www.caam.rice.edu/tech_reports/1990/TR90-17.pdf
    bool compositeNewton = false;

    // If the composite Newton method is to be used, typical formulae assume
    // that the correction involves a strictly feasible point. But this is
    // often a gross error for challenging systems (e.g., netlib/LP_data).
    bool compositeNewtonAssumeFeasible = false;

    // For determining the ratio of the amount to balance the affine and
    // correction updates. The other common option is
    // 'PredictorCorrectorCentrality'.
    function<Real(Real,Real,Real,Real)> centralityRule =
      StepLengthCentrality<Real>;

    // Use a simple shift for forcing cone membership during initialization?
    bool standardInitShift = true;

    // Force the primal and dual step lengths to be the same size?
    // Allowing different step lengths can substantially decrease the number
    // of iterations for Linear Programming.
    bool forceSameStep = false;

    // The controls for quasi-(semi)definite solves
    RegSolveCtrl<Real> solveCtrl;

    // Wrap the Interior Point Method with an equilibration.
    // This should almost always be set to true.
    bool outerEquil = true;

    // The size of the Krylov subspace used for loosely estimating two-norms of
    // sparse matrices.
    Int twoNormKrylovBasisSize = 6;

    // Print the progress of the Interior Point Method?
    bool print = false;

    // Time the components of the Interior Point Method?
    bool time = false;

    // Rather than solving the prescribed conic programs, we could add in
    // "small" amounts of permanent regularization for the primal and dual
    // variables to yield an augmented Lagrangian of the form
    //
    //   L(x,s;y,z) = c^T x + y^T (A x - b) + z^T (G x + s - h) +
    //                + (1/2) xRegSmall || x - x_0 ||_2^2
    //                + (1/2)           || Sqrt(Gamma_s) (s - s_0) ||_2^2
    //                - (1/2) yRegSmall || y - y_0 ||_2^2
    //                - (1/2) zRegSmall || z - z_0 ||_2^2
    //                + mu Phi(s),
    //
    // where (x_0,y_0,z_0) will typically be set to the current estimate of the
    // solution, but Gamma_s and s_0 are only defined implicitly so that
    // 'z + Gamma_s (s - s_0)' is entrywise at least as large as
    // 'zMinPivotValue'.
    //
    // The regularization of (x,y,z) ensures that the reduced KKT system
    // is symmetric quasi-definite, whereas the regularization of s ensures
    // that pivoting on the (s,s) block of the full KKT system divides by
    // 'z + Gamma_s (s - s_0)' rather than 'z'.
    //
    Real xRegSmallLogEps = Real(0.8);
    Real yRegSmallLogEps = Real(0.8);
    Real zRegSmallLogEps = Real(0.8);
    Real xRegLargeLogEps = Real(0.6);
    Real yRegLargeLogEps = Real(0.6);
    Real zRegLargeLogEps = Real(0.6);
    // The convergence of 80BAU3B in double-precision (with the other
    // parameters held fixed) is severely damaged if this value is set to
    // something as strict as 1.0.
    Real zMinPivotValueLogEps = Real(1.5);

    // If it turns out that the "large" regularization cannot be resolved,'
    // we will increase it by the following factor at each iteration.
    Real regIncreaseFactorLogEps = Real(-0.02);

    // If the maximum ratio between the primary and dual variables exceeds this
    // value, a centering step is taken instead of a predictor-corrector
    // step.
    Real maxComplementRatio = Real(1000);
    bool softDualityTargets = true;
    Real lowerTargetRatioLogCompRatio = Real(-0.25);
    Real upperTargetRatioLogCompRatio = Real( 0.25);
    Real lowerTargetRatioLogMaxCompRatio = Real(-0.6);
    Real upperTargetRatioLogMaxCompRatio = Real( 0.6);

    // When solving for the combined direction, previous iterates can be used
    // to accurately estimate the norms of the x, y, and z portions of the
    // solution to aid in rescaling. This variable bounds the maximal rescaling
    // that we allow ourselves to perform there.
    Real maxRescaleRatioLogEps = Real(-0.5);

    // At each iteration, the primal and dual variables are rescaled so that
    // they live within the following two-norm ranges.
    Real primalNormLowerBound = Real(0.1);
    Real primalNormUpperBound = Real(10);
    Real dualNormLowerBound = Real(0.1);
    Real dualNormUpperBound = Real(10);
    Real backoffScalePower = Real(-0.5);

    // Initially attempt to solve with only the "small" regularization by
    // preconditioning with a factorization involving the "large"
    // regularization? If this two-stage procedure breaks down before the
    // minimum tolerance is met, we will not attempt to resolve down to the
    // "small" regularization level.
    //
    // DEPRECATED for LP's and QP's.
    bool twoStage = false;

    // A lower bound on the maximum entry in the Nesterov-Todd scaling point
    // before ad-hoc procedures to enforce the cone constraints should be
    // employed.
    //
    // DEPRECATED for LP's and QP's
    Real wSafeMaxNorm = Pow(limits::Epsilon<Real>(),Real(-0.15));

    // Equilibrating before factoring in the two-stage scheme can prevent the
    // iterative solver from converging due to small errors in the equilibrated
    // scale being very large (and of large rank) in the original scale.
    // Further, equilibration in the single-stage scheme seems to lead to a few
    // more iterations for PILOT87.
    //
    // DEPRECATED for LP's and QP's.
    bool equilibrateIfSingleStage = false;

    // If the Nesterov-Todd scaling point has an entry of magnitude greater than
    // the following and the minimum tolerance has been achieved, simply stop
    // and declare success. This is meant to prevent expensive (equilibrated)
    // further steps which are too polluted with floating-point error to
    // make substantial progress.
    //
    // DEPRECATED for LP's and QP's.
    Real wMaxLimit = Pow(limits::Epsilon<Real>(),Real(-0.4));

    // If the maximum entry in the NT scaling point is larger in magnitude than
    // this value, then use Ruiz equilibration on the KKT system (with the
    // specified limit on the number of iterations).
    //
    // DEPRECATED for LP's and QP's.
    Real ruizEquilTol = Pow(limits::Epsilon<Real>(),Real(-0.25));
    Int ruizMaxIter = 3;

    // If Ruiz equilibration was not performed, but the max norm of the NT
    // scaling point is larger than this value, then use diagonal equilibration
    // for solving the KKT system.
    //
    // DEPRECATED for LP's and QP's.
    Real diagEquilTol = Pow(limits::Epsilon<Real>(),Real(-0.15));
};

template<typename Real>
struct IPMInfo
{
    Real primalObjective;
    Real dualObjective;

    Real infeasibilityError;
    Real relativeObjectiveGap;
    Real relativeComplementarityGap;

    Int numIterations;
};

// Alternating Direction Method of Multipliers
// ===========================================
template<typename Real>
struct ADMMCtrl
{
    Real rho=Real(1);
    Real alpha=Real(1.2);
    Int maxIter=500;
    // TODO(poulson): Base upon machine epsilon?
    Real absTol=Real(1e-6);
    Real relTol=Real(1e-4);
    bool inv=true;
    bool print=true;
};

} // namespace El

#endif // ifndef EL_OPTIMIZATION_SOLVERS_UTIL_HPP
