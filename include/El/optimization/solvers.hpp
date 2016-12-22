/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_OPTIMIZATION_SOLVERS_HPP
#define EL_OPTIMIZATION_SOLVERS_HPP

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

// Linear program
// ==============

namespace LPApproachNS {
enum LPApproach {
  LP_ADMM,
  LP_MEHROTRA
};
} // namespace LPApproachNS
using namespace LPApproachNS;

namespace lp {

namespace direct {

// Attempt to solve a pair of Linear Programs in "direct" conic form:
//
//   min c^T x,
//   s.t. A x = b, x >= 0
//
//   max -b^T y
//   s.t. A^T y -z + c = 0, z >= 0
//

// Control structure for the high-level "direct" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    LPApproach approach=LP_MEHROTRA;
    ADMMCtrl<Real> admmCtrl;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl( bool isSparse )
    { mehrotraCtrl.system = ( isSparse ? AUGMENTED_KKT : NORMAL_KKT ); }
};

} // namespace direct

namespace affine {

// Attempt to solve a pair of Linear Programs in "affine" conic form:
//
//   min c^T x,
//   s.t. A x = b, G x + s = h, s >= 0
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z >= 0
//

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    LPApproach approach=LP_MEHROTRA;
    MehrotraCtrl<Real> mehrotraCtrl;
};

} // namespace affine

} // namespace lp

// Direct conic form
// -----------------
template<typename MatrixType,typename VectorType>
struct DirectLPProblem
{
    // The objective is 'c^T x'.
    VectorType c;

    // The primal equality constraint is 'A x = b; x >= 0'.
    MatrixType A;
    VectorType b;
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const Grid& grid )
{
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    problem.c.Align(0,0);
    problem.A.Align(0,0);
    problem.b.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const Grid& grid )
{
    problem.c.SetComm( grid.Comm() );
    problem.A.SetComm( grid.Comm() );
    problem.b.SetComm( grid.Comm() );
}

template<typename Real>
bool SimpleAlignments
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem )
{
    return problem.c.ColAlign() == 0 && problem.c.RowAlign() == 0 &&
           problem.A.ColAlign() == 0 && problem.A.RowAlign() == 0 &&
           problem.b.ColAlign() == 0 && problem.b.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem )
{ return true; }

template<typename VectorType>
struct DirectLPSolution
{
    // The primal solution is 'x', where 'A x = b, x >= 0'.
    VectorType x;

    // The dual solution is '(y,z)', where 'A^T y - z + c = 0, z >= 0'.
    VectorType y;
    VectorType z;
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( DirectLPSolution<Matrix<Real>>& solution, const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( DirectLPSolution<DistMatrix<Real>>& solution, const Grid& grid )
{
    solution.x.SetGrid( grid );
    solution.y.SetGrid( grid );
    solution.z.SetGrid( grid );
    solution.x.Align(0,0);
    solution.y.Align(0,0);
    solution.z.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( DirectLPSolution<DistMultiVec<Real>>& solution, const Grid& grid )
{
    solution.x.SetComm( grid.Comm() );
    solution.y.SetComm( grid.Comm() );
    solution.z.SetComm( grid.Comm() );
}

// Quack...
template<typename Real>
bool SimpleAlignments
( const DirectLPSolution<Matrix<Real>>& solution )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const DirectLPSolution<DistMatrix<Real>>& solution )
{
    return solution.x.ColAlign() == 0 && solution.x.RowAlign() == 0 &&
           solution.y.ColAlign() == 0 && solution.y.RowAlign() == 0 &&
           solution.z.ColAlign() == 0 && solution.z.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const DirectLPSolution<DistSparseMatrix<Real>>& solution )
{ return true; }

template<typename VectorType>
struct DirectLPResidual
{
    VectorType primalEquality; // This residual is 'A x - b'.

    VectorType dualEquality; // This residual is 'A^T y - z + c'.
    VectorType dualConic; // This residual is 'x o z'.
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( DirectLPResidual<Matrix<Real>>& residual, const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( DirectLPResidual<DistMatrix<Real>>& residual, const Grid& grid )
{
    residual.primalEquality.SetGrid( grid );
    residual.dualEquality.SetGrid( grid );
    residual.dualConic.SetGrid( grid );
    residual.primalEquality.Align(0,0);
    residual.dualEquality.Align(0,0);
    residual.dualConic.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( DirectLPResidual<DistMultiVec<Real>>& residual, const Grid& grid )
{
    residual.primalEquality.SetComm( grid.Comm() );
    residual.dualEquality.SetComm( grid.Comm() );
    residual.dualConic.SetComm( grid.Comm() );
}

template<typename Real>
bool SimpleAlignments
( const DirectLPResidual<DistMatrix<Real>>& residual )
{
    return residual.primalEquality.ColAlign() == 0 &&
           residual.primalEquality.RowAlign() == 0 &&
           residual.dualEquality.ColAlign() == 0 &&
           residual.dualEquality.RowAlign() == 0 &&
           residual.dualConic.ColAlign() == 0 &&
           residual.dualConic.RowAlign() == 0;
}

template<typename Real>
void LP
( const DirectLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
void LP
( const DirectLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        DirectLPSolution<DistMatrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
void LP
( const DirectLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        DirectLPSolution<Matrix<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );
template<typename Real>
void LP
( const DirectLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        DirectLPSolution<DistMultiVec<Real>>& solution,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );

// These interfaces are now deprecated in favor of the above.
template<typename Real>
[[deprecated]]
void LP
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
[[deprecated]]
void LP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(false) );
template<typename Real>
[[deprecated]]
void LP
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );
template<typename Real>
[[deprecated]]
void LP
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const lp::direct::Ctrl<Real>& ctrl=lp::direct::Ctrl<Real>(true) );

// Affine conic form
// -----------------
template<typename MatrixType,typename VectorType>
struct AffineLPProblem
{
    // The objective is 'c^T x'.
    VectorType c;

    // The primal equality constraint is 'A x = b'.
    MatrixType A;
    VectorType b;

    // The primal cone constraint is 'G x + s = h, s >= 0'.
    MatrixType G;
    VectorType h;
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const Grid& grid )
{
    problem.c.SetGrid( grid );
    problem.A.SetGrid( grid );
    problem.b.SetGrid( grid );
    problem.G.SetGrid( grid );
    problem.h.SetGrid( grid );
    problem.c.Align(0,0);
    problem.A.Align(0,0);
    problem.b.Align(0,0);
    problem.G.Align(0,0);
    problem.h.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
  const Grid& grid )
{
    problem.c.SetComm( grid.Comm() );
    problem.A.SetComm( grid.Comm() );
    problem.b.SetComm( grid.Comm() );
    problem.G.SetComm( grid.Comm() );
    problem.h.SetComm( grid.Comm() );
}

// Quack...
template<typename Real>
bool SimpleAlignments
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem )
{ }
template<typename Real>
bool SimpleAlignments
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem )
{
    return problem.c.ColAlign() == 0 && problem.c.RowAlign() == 0 &&
           problem.A.ColAlign() == 0 && problem.A.RowAlign() == 0 &&
           problem.b.ColAlign() == 0 && problem.b.RowAlign() == 0 &&
           problem.G.ColAlign() == 0 && problem.G.RowAlign() == 0 &&
           problem.h.ColAlign() == 0 && problem.h.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem )
{ }
template<typename Real>
bool SimpleAlignments
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem )
{ }

template<typename VectorType>
struct AffineLPSolution
{
    // The primal solution is '(x,s)', where 'A x = b, G x + s = h, s >= 0'.
    VectorType x;
    VectorType s;

    // The dual solution is '(y,z)', where 'A^T y + G^T z + c = 0, z >= 0'.
    VectorType y;
    VectorType z;
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( AffineLPSolution<Matrix<Real>>& solution, const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( AffineLPSolution<DistMatrix<Real>>& solution, const Grid& grid )
{
    solution.x.SetGrid( grid );
    solution.s.SetGrid( grid );
    solution.y.SetGrid( grid );
    solution.z.SetGrid( grid );
    solution.x.Align(0,0);
    solution.s.Align(0,0);
    solution.y.Align(0,0);
    solution.z.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( AffineLPSolution<DistMultiVec<Real>>& solution, const Grid& grid )
{
    solution.x.SetComm( grid.Comm() );
    solution.s.SetComm( grid.Comm() );
    solution.y.SetComm( grid.Comm() );
    solution.z.SetComm( grid.Comm() );
}

// Quack...
template<typename Real>
bool SimpleAlignments
( const AffineLPSolution<Matrix<Real>>& solution )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const AffineLPSolution<DistMatrix<Real>>& solution )
{
    return solution.x.ColAlign() == 0 && solution.x.RowAlign() == 0 &&
           solution.y.ColAlign() == 0 && solution.y.RowAlign() == 0 &&
           solution.z.ColAlign() == 0 && solution.z.RowAlign() == 0 &&
           solution.s.ColAlign() == 0 && solution.s.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const AffineLPSolution<DistMultiVec<Real>>& solution )
{ return true; }

template<typename VectorType>
struct AffineLPResidual
{
    VectorType primalEquality; // This residual is 'A x - b'.
    VectorType primalConic; // This residual is 'G x + s - h'.

    VectorType dualEquality; // This residual is 'A^T y + G^T z + c'.
    VectorType dualConic; // This residual is 's o z'.
};

// Quack...
template<typename Real>
void ForceSimpleAlignments
( AffineLPResidual<Matrix<Real>>& residual, const Grid& grid )
{ }
template<typename Real>
void ForceSimpleAlignments
( AffineLPResidual<DistMatrix<Real>>& residual, const Grid& grid )
{
    residual.primalEquality.SetGrid( grid );
    residual.primalConic.SetGrid( grid );
    residual.dualEquality.SetGrid( grid );
    residual.dualConic.SetGrid( grid );
    residual.primalEquality.Align(0,0);
    residual.primalConic.Align(0,0);
    residual.dualEquality.Align(0,0);
    residual.dualConic.Align(0,0);
}
template<typename Real>
void ForceSimpleAlignments
( AffineLPResidual<DistMultiVec<Real>>& residual, const Grid& grid )
{
    residual.primalEquality.SetComm( grid.Comm() );
    residual.primalConic.SetComm( grid.Comm() );
    residual.dualEquality.SetComm( grid.Comm() );
    residual.dualConic.SetComm( grid.Comm() );
}

// Quack...
template<typename Real>
bool SimpleAlignments
( const AffineLPResidual<Matrix<Real>>& residual )
{ return true; }
template<typename Real>
bool SimpleAlignments
( const AffineLPResidual<DistMatrix<Real>>& residual )
{
    return residual.primalEquality.ColAlign() == 0 &&
           residual.primalEquality.RowAlign() == 0 &&
           residual.primalConic.ColAlign() == 0 &&
           residual.primalConic.RowAlign() == 0 &&
           residual.dualEquality.ColAlign() == 0 &&
           residual.dualEquality.RowAlign() == 0 &&
           residual.dualConic.ColAlign() == 0 &&
           residual.dualConic.RowAlign() == 0;
}
template<typename Real>
bool SimpleAlignments
( const AffineLPResidual<DistMultiVec<Real>>& residual )
{ return true; }

template<typename Real>
void LP
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
void LP
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );

// These interfaces are now deprecated in favor of the above.
template<typename Real>
[[deprecated]]
void LP
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
[[deprecated]]
void LP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
[[deprecated]]
void LP
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );
template<typename Real>
[[deprecated]]
void LP
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const lp::affine::Ctrl<Real>& ctrl=lp::affine::Ctrl<Real>() );

// Quadratic program
// =================

namespace QPApproachNS {
enum QPApproach {
  QP_ADMM,
  QP_MEHROTRA
};
} // namespace QPApproachNS
using namespace QPApproachNS;

namespace qp {

namespace direct {

// Attempt to solve a pair of Quadratic Programs in "direct" conic form:
//
//   min (1/2) x^T Q x + c^T x,
//   s.t. A x = b, x >= 0
//
//   max (1/2) (A^T y - z + c)^T pinv(Q) (A^T y -z + c) - b^T y
//   s.t. A^T y - z + c in range(Q), z >= 0
//

// Control structure for the high-level "direct" conic-form QP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    QPApproach approach=QP_MEHROTRA;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl() { mehrotraCtrl.system = AUGMENTED_KKT; }
};

} // namespace direct

namespace affine {

// Attempt to solve a pair of Quadratic Programs in "affine" conic form:
//
//   min (1/2) x^T Q x + c^T x,
//   s.t. A x = b, G x + s = h, s >= 0
//
//   max (1/2) (A^T y + G^T z + c)^T pinv(Q) (A^T y + G^T z + c)  -b^T y - h^T z
//   s.t. A^T y + G^T z + c in range(Q), z >= 0
//

// Control structure for the high-level "affine" conic-form QP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    QPApproach approach=QP_MEHROTRA;
    MehrotraCtrl<Real> mehrotraCtrl;
};

} // namespace affine

namespace box {

// Solve (a set of) quadratic programs of the form
//   min (1/2) x' Q x + c' x, subject to l_b <= x <= u_b
//    x

// Alternating Direction Method of Multipliers
// -------------------------------------------
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Int ADMM
( const Matrix<Real>& Q,
  const Matrix<Real>& C,
        Real lb,
        Real ub,
        Matrix<Real>& X,
  const ADMMCtrl<Real>& ctrl=ADMMCtrl<Real>() );
template<typename Real,
         typename=EnableIf<IsReal<Real>>>
Int ADMM
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& C,
        Real lb,
        Real ub,
        AbstractDistMatrix<Real>& X,
  const ADMMCtrl<Real>& ctrl=ADMMCtrl<Real>() );

} // namespace box

} // namespace qp

// Direct conic form
// -----------------
template<typename Real>
void QP
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );
template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const qp::direct::Ctrl<Real>& ctrl=qp::direct::Ctrl<Real>() );

// Affine conic form
// -----------------
template<typename Real>
void QP
( const Matrix<Real>& Q,
  const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const AbstractDistMatrix<Real>& Q,
  const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );
template<typename Real>
void QP
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const qp::affine::Ctrl<Real>& ctrl=qp::affine::Ctrl<Real>() );

// Second-order Cone Program
// =========================
namespace SOCPApproachNS {
enum SOCPApproach {
  SOCP_ADMM,     // NOTE: Not yet supported
  SOCP_MEHROTRA
};
} // namespace SOCPApproachNS
using namespace SOCPApproachNS;

namespace socp {
namespace direct {

// Attempt to solve a pair of Second-Order Cone Programs in "direct" conic form:
//
//   min c^T x,
//   s.t. A x = b, x in K,
//
//   max -b^T y
//   s.t. A^T y - z + c = 0, z in K,
//
// where the cone K is a product of second-order cones.
//

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    SOCPApproach approach=SOCP_MEHROTRA;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl()
    {
        mehrotraCtrl.system = AUGMENTED_KKT;
        mehrotraCtrl.minTol = Pow(limits::Epsilon<Real>(),Real(0.25));
        mehrotraCtrl.targetTol = Pow(limits::Epsilon<Real>(),Real(0.5));
    }
};

} // namespace direct

namespace affine {

// Attempt to solve a pair of Second-Order Cone Programs in "affine" conic form:
//
//   min c^T x,
//   s.t. A x = b, G x + s = h, s in K,
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z in K,
//
// where the cone K is a product of second-order cones.
//

// Control structure for the high-level "affine" conic-form LP solver
// ------------------------------------------------------------------
template<typename Real>
struct Ctrl
{
    SOCPApproach approach=SOCP_MEHROTRA;
    MehrotraCtrl<Real> mehrotraCtrl;

    Ctrl()
    {
        mehrotraCtrl.minTol = Pow(limits::Epsilon<Real>(),Real(0.25));
        mehrotraCtrl.targetTol = Pow(limits::Epsilon<Real>(),Real(0.5));
    }
};

} // namespace affine
} // namespace socp

// Direct conic form
// -----------------
template<typename Real>
void SOCP
( const Matrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const SparseMatrix<Real>& A,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );
template<typename Real>
void SOCP
( const DistSparseMatrix<Real>& A,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
  const socp::direct::Ctrl<Real>& ctrl=socp::direct::Ctrl<Real>() );

// Affine conic form
// -----------------
template<typename Real>
void SOCP
( const Matrix<Real>& A,
  const Matrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const AbstractDistMatrix<Real>& A,
  const AbstractDistMatrix<Real>& G,
  const AbstractDistMatrix<Real>& b,
  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
  const AbstractDistMatrix<Int>& orders,
  const AbstractDistMatrix<Int>& firstInds,
        AbstractDistMatrix<Real>& x,
        AbstractDistMatrix<Real>& y,
        AbstractDistMatrix<Real>& z,
        AbstractDistMatrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const Matrix<Real>& h,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );
template<typename Real>
void SOCP
( const DistSparseMatrix<Real>& A,
  const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const socp::affine::Ctrl<Real>& ctrl=socp::affine::Ctrl<Real>() );

} // namespace El

#endif // ifndef EL_OPTIMIZATION_SOLVERS_HPP
