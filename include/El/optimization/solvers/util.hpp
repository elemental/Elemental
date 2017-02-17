/*
   Copyright (c) 2009-2016, Jack Poulson
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

// Mehrotra's Predictor-Corrector Infeasible Interior Point Method
// ===============================================================
template<typename Real>
inline Real StepLengthCentrality
( Real mu, Real muAff, Real alphaAffPri, Real alphaAffDual )
{ return Pow(1-Min(alphaAffPri,alphaAffDual),Real(3)); }

template<typename Real>
inline Real MehrotraCentrality
( Real mu, Real muAff, Real alphaAffPri, Real alphaAffDual )
{ return Min(Pow(muAff/mu,Real(3)),Real(1)); }

template<typename Real>
struct MehrotraCtrl
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
    bool primalInit=false, dualInit=false;

    // Throw an exception if this tolerance could not be achieved.
    Real minTol=Pow(limits::Epsilon<Real>(),Real(0.3));

    // Exit the Interior Point Methods if this tolerance has been achieved.
    Real targetTol=Pow(limits::Epsilon<Real>(),Real(0.5));

    // The maximum number of iterations of the IPM. This should only be
    // activated in pathological circumstances, as even 100 iterations of an
    // IPM is excessive.
    Int maxIts=1000;

    // If a step of length alpha would hit the boundary, instead step a distance
    // of 'maxStepRatio*alpha'.
    Real maxStepRatio=Real(0.99);

    // For configuring how many reductions of the first-order optimality
    // conditions should be performed before solving a linear system.
    // For problems with 'affine' cone constraints, i.e., (h - G x) in K,
    // this should remain as FULL_KKT. For sparse problems with 'direct' cone
    // constraints, i.e., x in K, both AUGMENTED_KKT (use a QSD solver) and
    // NORMAL_KKT (use a Cholesky solver) are also possible. The latter should
    // be avoided when the normal equations are sufficiently denser than the
    // (larger) augmented formulation.
    KKTSystem system=FULL_KKT;

    // Use Mehrotra's second-order corrector?
    // TODO(poulson): Add support for Gondzio's correctors
    bool mehrotra=true;

    // For determining the ratio of the amount to balance the affine and 
    // correction updates. The other common option is 'MehrotraCentrality'.
    function<Real(Real,Real,Real,Real)>
      centralityRule=StepLengthCentrality<Real>;

    // Use a simple shift for forcing cone membership during initialization?
    bool standardInitShift=true;

    // If the maximum ratio between the primary and dual variables exceeds this
    // value, the barrier parameter is kept at its previous value to attempt to
    // increase the centrality.
    Real balanceTol=Pow(limits::Epsilon<Real>(),Real(-0.19));

    // Force the primal and dual step lengths to be the same size?
    bool forceSameStep=true;

    // The controls for quasi-(semi)definite solves
    RegSolveCtrl<Real> solveCtrl;

    // Always use an iterative solver to resolve the regularization?
    // TODO(poulson): Generalize this to a strategy (e.g., resolve if
    // stagnating), as this choice has been observed to substantially impact
    // both the cost and number of iterations.
    bool resolveReg=true;

    // Wrap the Interior Point Method with an equilibration.
    // This should almost always be set to true.
    bool outerEquil=true;

    // The size of the Krylov subspace used for loosely estimating two-norms of
    // sparse matrices.
    Int basisSize = 6;

    // Print the progress of the Interior Point Method?
    bool print=false;

    // Time the components of the Interior Point Method?
    bool time=false;

    // A lower bound on the maximum entry in the Nesterov-Todd scaling point
    // before ad-hoc procedures to enforce the cone constraints should be
    // employed.
    Real wSafeMaxNorm=Pow(limits::Epsilon<Real>(),Real(-0.15));

    // If the Nesterov-Todd scaling point has an entry of magnitude greater than
    // the following and the minimum tolerance has been achieved, simply stop
    // and declare success. This is meant to prevent expensive (equilibrated)
    // further steps which are too polluted with floating-point error to
    // make substantial progress.
    Real wMaxLimit=Pow(limits::Epsilon<Real>(),Real(-0.4));

    // If the maximum entry in the NT scaling point is larger in magnitude than
    // this value, then use Ruiz equilibration on the KKT system (with the
    // specified limit on the number of iterations).
    Real ruizEquilTol=Pow(limits::Epsilon<Real>(),Real(-0.25));
    Int ruizMaxIter=3;

    // If Ruiz equilibration was not performed, but the max norm of the NT
    // scaling point is larger than this value, then use diagonal equilibration
    // for solving the KKT system.
    Real diagEquilTol=Pow(limits::Epsilon<Real>(),Real(-0.15));

    // Whether or not additional matrix-vector multiplications should be
    // performed in order to check the accuracy of the solution to each
    // block row of the KKT system.
#ifdef EL_RELEASE
    bool checkResiduals=false;
#else
    bool checkResiduals=true;
#endif

    // Temporary regularization for the primal, dual, and dual slack variables
    Real reg0Tmp = Pow(limits::Epsilon<Real>(),Real(0.25));
    Real reg1Tmp = Pow(limits::Epsilon<Real>(),Real(0.25));
    Real reg2Tmp = Pow(limits::Epsilon<Real>(),Real(0.25));

    // Permanent regularization for the primal, dual, and dual slack variables
    Real reg0Perm = Pow(limits::Epsilon<Real>(),Real(0.35));
    Real reg1Perm = Pow(limits::Epsilon<Real>(),Real(0.35));
    Real reg2Perm = Pow(limits::Epsilon<Real>(),Real(0.35));

    // TODO(poulson): Add a user-definable (muAff,mu) -> sigma function to
    // replace the default, (muAff/mu)^3
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
