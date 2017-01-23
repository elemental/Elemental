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

// TODO(poulson): Move this into a central location.
template<typename Real>
Real RelativeComplementarityGap
( const Real& primalObj, const Real& dualObj, const Real& dualityProduct )
{
    EL_DEBUG_CSE
    Real relCompGap;
    if( primalObj < Real(0) )
        relCompGap = dualityProduct / -primalObj;
    else if( dualObj > Real(0) )
        relCompGap = dualityProduct / dualObj;
    else
        relCompGap = 2; // 200% error if the signs differ inadmissibly.
    return relCompGap;
}

// TODO(poulson): Move this into a central location.
template<typename Real>
Real RelativeObjectiveGap
( const Real& primalObj, const Real& dualObj, const Real& dualityProduct )
{
    EL_DEBUG_CSE
    const Real relObjGap =
      Abs(primalObj-dualObj) / (Max(Abs(primalObj),Abs(dualObj))+1);
    return relObjGap;
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
// which corresponds to G = -I and h = 0.
//
// We make use of the regularized Lagrangian
//
//   L(x,s;y,z) = c^T x + y^T (A x - b) + z^T (G x + s - h)
//                + (1/2) gamma_x || x - x_0 ||_2^2
//                + (1/2)         || Sqrt(Gamma_s) (s - s_0) ||_2^2
//                - (1/2) gamma_y || y - y_0 ||_2^2
//                - (1/2) gamma_z || z - z_0 ||_2^2
//                + mu Phi(s),
//
// where we note that the two-norm regularization is positive for the primal
// variables x and s and *negative* for the dual variables y and z.
// The centering points (x_0,y_0,z_0) are typically chosen as the current
// estimates of the solution, whereas Gamma_s and s_0 are only implicitly
// defined and serve to force
//
//   z + Gamma_s (s - s_0) = LowerClip( z, zMinPivotValue ).
//
// The subsequent first-order optimality conditions for x, y, and z become
//
//   Nabla_x L = c + A^T y + G^T z + gamma_x (x - x_0) = 0,
//   Nabla_y L = A x - b - gamma_y (y - y_0) = 0,
//   Nabla_z L = G x + s - h - gamma_z (z - z_0) = 0.
//
// These can be arranged into the symmetric quasi-definite form
//
//   | gamma_x I,    A^T,         G^T     | | x | = | -c + gamma_x x_0  |
//   |     A,     -gamma_y I,      0      | | y |   |  b - gamma_y y_0  |
//   |     G,         0,       -gamma_z I | | z |   | h-s - gamma_z z_0 |.
//
//
// Rescalings of the problem can be understood in terms of the transformations
//
//   x = D_x x', y = D_y y', z = D_z z', s = D_s s',
//
// which lead to the rescaled Lagrangian
//
//   L(x',s';y',z') = (D_x c)^T x' + y'^T ((D_y A D_x) x' - (D_y b))
//     + z'^T ((D_z G D_x) x' + D_z D_s s' - (D_z h))
//     + (1/2) gamma_x || D_x x' - x_0 ||_2^2
//     + (1/2)         || Gamma_s (D_s s' - s_0) ||_2^2
//     - (1/2) gamma_y || D_y y' - y_0 ||_2^2
//     - (1/2) gamma_z || D_z z' - z_0 ||_2^2
//     + mu Phi(D_s s').
//
// One immediately notices that the oddity is the term "D_z D_s s'", which
// suggests that we must have that D_z D_s = ones, which would preserve the
// complementarity condition s o z = mu e. We therefore enforce the condition
// that D_s := inv(D_z) to arrive at
//
//   L(x',s';y',z') = (D_x c)^T x' + y'^T ((D_y A D_x) x' - (D_y b))
//     + z'^T ((D_z G D_x) x' + s' - (D_z h))
//     + (1/2) gamma_x ||    D_x   x' - x_0 ||_2^2
//     + (1/2)         || Gamma_s (inv(D_z) s' - s_0) ||_2^2
//     - (1/2) gamma_y ||    D_y   y' - y_0 ||_2^2
//     - (1/2) gamma_z ||    D_z   z' - z_0 ||_2^2
//     + mu Phi(inv(D_z) s').
//
// We can also make use of the fact that, for any primalScale != 0,
//
//     {A x = b, G x + s = h}  iff
//     {A (x/primalScale) = (b/primalScale),
//      G (x/primalScale) + (s/primalScale) = (h/primalScale)}
//
// and that, for any dualScale != 0,
//
//     {A^T y + G^T z + c = 0, z >= 0}  iff
//     {A^T (y/dualScale) + G^T (z/dualScale) + (c/dualScale) = 0,
//      (z/dualScale) >= 0}.
//
// Thus, we are free to insist that
//
//   max( || D_y b / primalScale ||_max, || D_z h / primalScale ||_max ) <= 1
//
// and
//
//   || D_x c / dualScale ||_max <= 1.
//
// Then, using the combined rescalings
//
//   x = primalScale D_x x~, s = primalScale inv(D_z) s~,
//   y = dualScale D_y y~, z = dualScale D_z z~,
//
// we arrive at the rescale, regularized Lagrangian
//
//   L(x~,s~;y~,z~) = (D_x c / dualScale)^T x~
//     + y~^T ((D_y A D_x) x~ - (D_y b / primalScale))
//     + z~^T ((D_z G D_x) x~ + s~ - (D_z h / primalScale))
//     + (1/2) gamma_x || primalScale D_x x~ - x_0 ||_2^2
//     + (1/2)         || Gamma_s (primalScale inv(D_z) s~ - s_0) ||_2^2
//     - (1/2) gamma_y || dualScale D_y y~ - y_0 ||_2^2
//     - (1/2) gamma_z || dualScale D_z z~ - z_0 ||_2^2
//     + mu Phi(primalScale inv(D_z) s').
//
// This can be simplified by setting
//
//   c~ = D_x c / dualScale,
//   A~ = D_y A D_x,
//   b~ = D_y b / primalScale,
//   G~ = D_z G D_x,
//   h~ = D_z h / primalScale,
//
// so that
//
//   L(x~,s~;y~,z~) = c~^T x~ + y~^T (A~ x~ - b~) + z~^T (G~ x~ + s~ - h~)
//     + (1/2) gamma_x || primalScale D_x x~ - x_0 ||_2^2
//     + (1/2)         || Gamma_s (primalScale inv(D_z) s~ - s_0) ||_2^2
//     - (1/2) gamma_y || dualScale D_y y~ - y_0 ||_2^2
//     - (1/2) gamma_z || dualScale D_z z~ - z_0 ||_2^2
//     + mu Phi(primalScale inv(D_z) s').
//
// In practice, we find (D_x,D_y,D_z) through (geometric) equilibration of the
// matrix [A; G] followed by rescaling the max-norms of the transformed b, c,
// and h down to one.
//
// And since any positive diagonal scaling is an automorphism of the positive
// orthant, we are free to remove it from, or compose it with, the barrier
// function. Further, since the two-norm regularization is added for stability
// reasons, it is acceptable to penalize the norms of the transformed
// variables rather than the original variables instead.
//
// Lastly, the gradient of the Lagrangian with respect to s being zero implies
//
//   Nabla_s L = z + Gamma_s (s - s_0) - mu inv(s) = 0.
//
// This can easily be arranged into the form
//
//   s o z + s o Gamma_s (s - s_0) = mu e,
//
// which implies that the 's' regularization about the point 's_0' modifies the
// four-by-four KKT system
//
//   | Z, 0,  0,   S  | | s | = | mu e |
//   | 0, 0, A^T, G^T | | x |   |  -c  |
//   | 0, A,  0,   0  | | y |   |   b  |
//   | I, G,  0,   0  | | z |   |  h-s |
//
// into
//
//   | Z + Gamma_s (S - S_0), 0,  0,   S  | | s | = | mu e |.
//   |          0,            0, A^T, G^T | | x |   |  -c  |
//   |          0,            A,  0,   0  | | y |   |   b  |
//   |          I,            G,  0,   0  | | z |   |  h-s |
//
// Since the upper-left block is typically an implicit initial pivot block,
// and many entries of z quickly converge towards zero, the complementarity
// condition implies that *either* 'z(i)' or 's(i)' is small, and so
// regularization of 's' about 's_0' should stabilize the solve as long as
// 's_0' is not exactly equal to 's_0'. 
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

template<typename Real>
struct DenseAffineLPEquilibration
{
    Real primalScale;
    Real dualScale;
    Matrix<Real> xScale;
    Matrix<Real> yScale;
    Matrix<Real> zScale;
};
template<typename Real>
struct DistDenseAffineLPEquilibration
{
    Real primalScale;
    Real dualScale;
    DistMatrix<Real,VC,STAR> xScale;
    DistMatrix<Real,VC,STAR> yScale;
    DistMatrix<Real,VC,STAR> zScale;
};
template<typename Real>
struct SparseAffineLPEquilibration
{
    Real primalScale;
    Real dualScale;
    Matrix<Real> xScale;
    Matrix<Real> yScale;
    Matrix<Real> zScale;
};
template<typename Real>
struct DistSparseAffineLPEquilibration
{
    Real primalScale;
    Real dualScale;
    DistMultiVec<Real> xScale;
    DistMultiVec<Real> yScale;
    DistMultiVec<Real> zScale;
};

template<typename Real>
void Equilibrate
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const AffineLPSolution<Matrix<Real>>& solution,
        AffineLPProblem<Matrix<Real>,Matrix<Real>>& equilibratedProblem,
        AffineLPSolution<Matrix<Real>>& equilibratedSolution,
        DenseAffineLPEquilibration<Real>& equilibration,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();

    equilibratedProblem = problem;
    equilibratedSolution = solution;

    // Equilibrate the LP by diagonally scaling [A;G].
    // Unfortunately, this routine returns the inverses of the desired scales,
    // so we manually invert them.
    Matrix<Real> rowScaleA, rowScaleG, colScale;
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      rowScaleA, rowScaleG, colScale,
      ctrl.print );
    Ones( equilibration.xScale, n, 1 );
    Ones( equilibration.yScale, m, 1 );
    Ones( equilibration.zScale, k, 1 );
    DiagonalSolve( LEFT, NORMAL, colScale,  equilibration.xScale );
    DiagonalSolve( LEFT, NORMAL, rowScaleA, equilibration.yScale );
    DiagonalSolve( LEFT, NORMAL, rowScaleG, equilibration.zScale );

    DiagonalScale
    ( LEFT, NORMAL, equilibration.xScale, equilibratedProblem.c );
    DiagonalScale
    ( LEFT, NORMAL, equilibration.yScale, equilibratedProblem.b );
    DiagonalScale
    ( LEFT, NORMAL, equilibration.zScale, equilibratedProblem.h );
    if( ctrl.primalInit )
    {
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.xScale, equilibratedSolution.x );
        DiagonalScale
        ( LEFT, NORMAL, equilibration.zScale, equilibratedSolution.s );
    }
    if( ctrl.dualInit )
    {
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.yScale, equilibratedSolution.y );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.zScale, equilibratedSolution.z );
    }

    // Rescale max(|| b ||_max,|| h ||_max) to roughly one.
    equilibration.primalScale =
      Max(Max(MaxNorm(equilibratedProblem.b),MaxNorm(equilibratedProblem.h)),
        Real(1e-6));
    equilibratedProblem.b *= Real(1)/equilibration.primalScale;
    equilibratedProblem.h *= Real(1)/equilibration.primalScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.primalScale;
        equilibratedSolution.s *= Real(1)/equilibration.primalScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.dualScale = Max(MaxNorm(equilibratedProblem.c),Real(1e-6));
    equilibratedProblem.c *= Real(1)/equilibration.dualScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.dualScale;
        equilibratedSolution.z *= Real(1)/equilibration.dualScale;
    }

    if( ctrl.print )
    {
        Output("primalScale = ",equilibration.primalScale);
        Output("dualScale   = ",equilibration.dualScale);
        Output("|| xScale ||_max = ",MaxNorm(equilibration.xScale));
        Output("|| yScale ||_max = ",MaxNorm(equilibration.yScale));
        Output("|| zScale ||_max = ",MaxNorm(equilibration.zScale));
    }
}

template<typename Real>
void Equilibrate
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
  const AffineLPSolution<DistMatrix<Real>>& solution,
        AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& equilibratedProblem,
        AffineLPSolution<DistMatrix<Real>>& equilibratedSolution,
        DistDenseAffineLPEquilibration<Real>& equilibration,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Grid& grid = problem.A.Grid();

    equilibratedProblem = problem;
    equilibratedSolution = solution;

    // Equilibrate the LP by diagonally scaling [A;G].
    // Unfortunately, this routine returns the inverses of the desired scales,
    // so we manually invert them.
    DistMatrix<Real,MR,STAR> colScale(grid);
    DistMatrix<Real,MC,STAR> rowScaleA(grid), rowScaleG(grid);
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      rowScaleA, rowScaleG, colScale,
      ctrl.print );
    Ones( equilibration.xScale, n, 1 );
    Ones( equilibration.yScale, m, 1 );
    Ones( equilibration.zScale, k, 1 );
    DiagonalSolve( LEFT, NORMAL, rowScaleA, equilibration.yScale );
    DiagonalSolve( LEFT, NORMAL, rowScaleG, equilibration.zScale );
    DiagonalSolve( LEFT, NORMAL, colScale,  equilibration.xScale );

    DiagonalScale
    ( LEFT, NORMAL, equilibration.xScale, equilibratedProblem.c );
    DiagonalScale
    ( LEFT, NORMAL, equilibration.yScale, equilibratedProblem.b );
    DiagonalScale
    ( LEFT, NORMAL, equilibration.zScale, equilibratedProblem.h );
    if( ctrl.primalInit )
    {
        DiagonalScale
        ( LEFT, NORMAL, equilibration.zScale, equilibratedSolution.s );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.xScale, equilibratedSolution.x );
    }
    if( ctrl.dualInit )
    {
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.yScale, equilibratedSolution.y );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.zScale, equilibratedSolution.z );
    }

    // Rescale max(|| b ||_max,|| h ||_max) to roughly one.
    equilibration.primalScale =
      Max(Max(MaxNorm(equilibratedProblem.b),MaxNorm(equilibratedProblem.h)),
        Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.primalScale;
    equilibratedProblem.h *= Real(1)/equilibration.primalScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.primalScale;
        equilibratedSolution.s *= Real(1)/equilibration.primalScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.dualScale = Max(MaxNorm(equilibratedProblem.c),Real(1));
    equilibratedProblem.c *= Real(1)/equilibration.dualScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.dualScale;
        equilibratedSolution.z *= Real(1)/equilibration.dualScale;
    }

    if( ctrl.print )
    {
        const Real xScaleMax = MaxNorm(equilibration.xScale);
        const Real yScaleMax = MaxNorm(equilibration.yScale);
        const Real zScaleMax = MaxNorm(equilibration.zScale);
        if( problem.A.Grid().Rank() == 0 )
        {
            Output("primalScale = ",equilibration.primalScale);
            Output("dualScale   = ",equilibration.dualScale);
            Output("|| xScale ||_max = ",xScaleMax);
            Output("|| yScale ||_max = ",yScaleMax);
            Output("|| zScale ||_max = ",zScaleMax);
        }
    }
}

template<typename Real>
void Equilibrate
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const AffineLPSolution<Matrix<Real>>& solution,
        AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& equilibratedProblem,
        AffineLPSolution<Matrix<Real>>& equilibratedSolution,
        SparseAffineLPEquilibration<Real>& equilibration,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();

    equilibratedProblem = problem;
    equilibratedSolution = solution;

    // Equilibrate the LP by diagonally scaling [A;G].
    // Unfortunately, this routine returns the inverses of the desired scales,
    // so we manually invert them.
    Matrix<Real> rowScaleA, rowScaleG, colScale;
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      rowScaleA, rowScaleG, colScale,
      ctrl.print );
    Ones( equilibration.xScale, n, 1 );
    Ones( equilibration.yScale, m, 1 );
    Ones( equilibration.zScale, k, 1 );
    DiagonalSolve( LEFT, NORMAL, colScale,  equilibration.xScale );
    DiagonalSolve( LEFT, NORMAL, rowScaleA, equilibration.yScale );
    DiagonalSolve( LEFT, NORMAL, rowScaleG, equilibration.zScale );

    DiagonalScale
    ( LEFT, NORMAL, equilibration.xScale, equilibratedProblem.c );
    DiagonalScale
    ( LEFT, NORMAL, equilibration.yScale, equilibratedProblem.b );
    DiagonalScale
    ( LEFT, NORMAL, equilibration.zScale, equilibratedProblem.h );
    if( ctrl.primalInit )
    {
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.xScale, equilibratedSolution.x );
        DiagonalScale
        ( LEFT, NORMAL, equilibration.zScale, equilibratedSolution.s );
    }
    if( ctrl.dualInit )
    {
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.yScale, equilibratedSolution.y );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.zScale, equilibratedSolution.z );
    }

    // Rescale max(|| b ||_max,|| h ||_max) to roughly one.
    equilibration.primalScale =
      Max(Max(MaxNorm(equilibratedProblem.b),MaxNorm(equilibratedProblem.h)),
        Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.primalScale;
    equilibratedProblem.h *= Real(1)/equilibration.primalScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.primalScale;
        equilibratedSolution.s *= Real(1)/equilibration.primalScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.dualScale = Max(MaxNorm(equilibratedProblem.c),Real(1));
    equilibratedProblem.c *= Real(1)/equilibration.dualScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.dualScale;
        equilibratedSolution.z *= Real(1)/equilibration.dualScale;
    }

    if( ctrl.print )
    {
        Output("primalScale = ",equilibration.primalScale);
        Output("dualScale   = ",equilibration.dualScale);
        Output("|| xcale  ||_2 = ",MaxNorm(equilibration.xScale));
        Output("|| yScale ||_2 = ",MaxNorm(equilibration.yScale));
        Output("|| zScale ||_2 = ",MaxNorm(equilibration.zScale));
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
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Grid& grid = problem.A.Grid();
    ForceSimpleAlignments( equilibratedProblem, grid );
    ForceSimpleAlignments( equilibratedSolution, grid );

    equilibratedProblem = problem;
    equilibratedSolution = solution;

    // Equilibrate the LP by diagonally scaling [A;G].
    // Unfortunately, this routine returns the inverses of the desired scales,
    // so we manually invert them.
    DistMultiVec<Real> rowScaleA(grid), rowScaleG(grid), colScale(grid);
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      rowScaleA, rowScaleG, colScale,
      ctrl.print );
    Ones( equilibration.xScale, n, 1 );
    Ones( equilibration.yScale, m, 1 );
    Ones( equilibration.zScale, k, 1 );
    DiagonalSolve( LEFT, NORMAL, colScale,  equilibration.xScale );
    DiagonalSolve( LEFT, NORMAL, rowScaleA, equilibration.yScale );
    DiagonalSolve( LEFT, NORMAL, rowScaleG, equilibration.zScale );

    DiagonalScale
    ( LEFT, NORMAL, equilibration.xScale, equilibratedProblem.c );
    DiagonalScale
    ( LEFT, NORMAL, equilibration.yScale, equilibratedProblem.b );
    DiagonalScale
    ( LEFT, NORMAL, equilibration.zScale, equilibratedProblem.h );
    if( ctrl.primalInit )
    {
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.xScale, equilibratedSolution.x );
        DiagonalScale
        ( LEFT, NORMAL, equilibration.zScale, equilibratedSolution.s );
    }
    if( ctrl.dualInit )
    {
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.yScale, equilibratedSolution.y );
        DiagonalSolve
        ( LEFT, NORMAL, equilibration.zScale, equilibratedSolution.z );
    }

    // Rescale max(|| b ||_max,|| h ||_max) to roughly one.
    equilibration.primalScale =
      Max(Max(MaxNorm(equilibratedProblem.b),MaxNorm(equilibratedProblem.h)),
        Real(1));
    equilibratedProblem.b *= Real(1)/equilibration.primalScale;
    equilibratedProblem.h *= Real(1)/equilibration.primalScale;
    if( ctrl.primalInit )
    {
        equilibratedSolution.x *= Real(1)/equilibration.primalScale;
        equilibratedSolution.s *= Real(1)/equilibration.primalScale;
    }

    // Rescale || c ||_max to roughly one.
    equilibration.dualScale = Max(MaxNorm(equilibratedProblem.c),Real(1));
    equilibratedProblem.c *= Real(1)/equilibration.dualScale;
    if( ctrl.dualInit )
    {
        equilibratedSolution.y *= Real(1)/equilibration.dualScale;
        equilibratedSolution.z *= Real(1)/equilibration.dualScale;
    }

    if( ctrl.print )
    {
        const Real xScaleMax = MaxNorm(equilibration.xScale);
        const Real yScaleMax = MaxNorm(equilibration.yScale);
        const Real zScaleMax = MaxNorm(equilibration.zScale);
        if( problem.A.Grid().Rank() == 0 )
        {
            Output("primalScale = ",equilibration.primalScale);
            Output("dualScale   = ",equilibration.dualScale);
            Output("|| xScale ||_max = ",xScaleMax);
            Output("|| yScale ||_max = ",yScaleMax);
            Output("|| zScale ||_max = ",zScaleMax);
        }
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

    // x = primalScale D_x x~
    solution.x *= equilibration.primalScale;
    DiagonalScale( LEFT, NORMAL, equilibration.xScale, solution.x );

    // s = primalScale inv(D_z) s~, as D_s = inv(D_z)
    solution.s *= equilibration.primalScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.zScale, solution.s );

    // y = dualScale D_y y~
    solution.y *= equilibration.dualScale;
    DiagonalScale( LEFT, NORMAL, equilibration.yScale, solution.y );

    // z = dualScale D_z z~
    solution.z *= equilibration.dualScale;
    DiagonalScale( LEFT, NORMAL, equilibration.zScale, solution.z );
}

template<typename Real>
void UndoEquilibration
( const AffineLPSolution<DistMatrix<Real>>& equilibratedSolution,
  const DistDenseAffineLPEquilibration<Real>& equilibration,
        AffineLPSolution<DistMatrix<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;

    // x = primalScale D_x x~
    solution.x *= equilibration.primalScale;
    DiagonalScale( LEFT, NORMAL, equilibration.xScale, solution.x );

    // s = primalScale inv(D_z) s~, as D_s = inv(D_z)
    solution.s *= equilibration.primalScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.zScale, solution.s );

    // y = dualScale D_y y~
    solution.y *= equilibration.dualScale;
    DiagonalScale( LEFT, NORMAL, equilibration.yScale, solution.y );

    // z = dualScale D_z z~
    solution.z *= equilibration.dualScale;
    DiagonalScale( LEFT, NORMAL, equilibration.zScale, solution.z );
}

template<typename Real>
void UndoEquilibration
( const AffineLPSolution<Matrix<Real>>& equilibratedSolution,
  const SparseAffineLPEquilibration<Real>& equilibration,
        AffineLPSolution<Matrix<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;

    // x = primalScale D_x x~
    solution.x *= equilibration.primalScale;
    DiagonalScale( LEFT, NORMAL, equilibration.xScale, solution.x );

    // s = primalScale inv(D_z) s~, as D_s = inv(D_z)
    solution.s *= equilibration.primalScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.zScale, solution.s );

    // y = dualScale D_y y~
    solution.y *= equilibration.dualScale;
    DiagonalScale( LEFT, NORMAL, equilibration.yScale, solution.y );

    // z = dualScale D_z z~
    solution.z *= equilibration.dualScale;
    DiagonalScale( LEFT, NORMAL, equilibration.zScale, solution.z );
}

template<typename Real>
void UndoEquilibration
( const AffineLPSolution<DistMultiVec<Real>>& equilibratedSolution,
  const DistSparseAffineLPEquilibration<Real>& equilibration,
        AffineLPSolution<DistMultiVec<Real>>& solution )
{
    EL_DEBUG_CSE
    solution = equilibratedSolution;

    // x = primalScale D_x x~
    solution.x *= equilibration.primalScale;
    DiagonalScale( LEFT, NORMAL, equilibration.xScale, solution.x );

    // s = primalScale inv(D_z) s~, as D_s = inv(D_z)
    solution.s *= equilibration.primalScale;
    DiagonalSolve( LEFT, NORMAL, equilibration.zScale, solution.s );

    // y = dualScale D_y y~
    solution.y *= equilibration.dualScale;
    DiagonalScale( LEFT, NORMAL, equilibration.yScale, solution.y );

    // z = dualScale D_z z~
    solution.z *= equilibration.dualScale;
    DiagonalScale( LEFT, NORMAL, equilibration.zScale, solution.z );
}

template<typename Real>
void EquilibratedIPM
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
            EL_DEBUG_ONLY(SetCallStack(callStack))
            return false;
        }
        return true;
      };

    // TODO(poulson): Set up 'increaseRegularization' function and call it
    // when factorizations or solves fail.

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
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real relObjGap =
          RelativeObjectiveGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relCompGap, relObjGap );

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
        dimacsError = Max(infeasError,maxRelGap);
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
             "  scaled relative duality gap = ",maxRelGap);
        }

        const bool metTolerances =
          infeasError <= ctrl.infeasibilityTol &&
          relCompGap <= ctrl.relativeComplementarityGapTol &&
          relObjGap <= ctrl.relativeObjectiveGapTol;
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

        // Factor the KKT system
        // =====================
        KKT( problem.A, problem.G, solution.s, solution.z, J );
        if( !attemptToFactor() )
        {
            // TODO(poulson): Increase regularization and continue instead.
            break;
        }

        // Compute the affine search direction
        // ===================================
        // r_mu := s o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z,
          d );
        if( !attemptToSolve(d) )
        {
            // TODO(poulson): Increase regularization and continue instead.
            break;
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
        {
            // TODO(poulson): Increase regularization and continue instead.
            break;
        }
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
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
        EquilibratedIPM( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedIPM( problem, solution, ctrl );
    }
}

template<typename Real>
void EquilibratedIPM
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real relObjGap =
          RelativeObjectiveGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relCompGap, relObjGap );
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
        dimacsError = Max(infeasError,maxRelGap);
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
                 "  scaled relative gap = ",maxRelGap);
        }

        const bool metTolerances =
          infeasError <= ctrl.infeasibilityTol &&
          relCompGap <= ctrl.relativeComplementarityGapTol &&
          relObjGap <= ctrl.relativeObjectiveGapTol;
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

        // Factor the KKT system
        // =====================
        KKT( problem.A, problem.G, solution.s, solution.z, J );
        if( !attemptToFactor() )
        {
            // TODO(poulson): Increase regularization and continue.
            break;
        }

        // Compute the affine search direction
        // ===================================
        // r_mu := s o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

        // Solve for the proposed step
        // ---------------------------
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        if( !attemptToSolve(d) )
        {
            // TODO(poulson): Increase regularization and continue.
            break;
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
        {
            // TODO(poulson): Increase regularization and continue.
            break;
        }
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
( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem,
        AffineLPSolution<DistMatrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
            AffineLPSolution<DistMatrix<Real>> alignedSolution;
            ForceSimpleAlignments( alignedSolution, grid );
            alignedSolution = solution;
            EquilibratedIPM( problem, alignedSolution, ctrl );
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
            EquilibratedIPM( alignedProblem, solution, ctrl );
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
            EquilibratedIPM( alignedProblem, alignedSolution, ctrl );
            solution = alignedSolution;
        }
    }
}

template<typename Real>
void EquilibratedIPM
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& origProblem,
        SparseAffineLPEquilibration<Real>& equilibration,
        AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
    const Real twoNormEstA =
      TwoNormEstimate( problem.A, ctrl.twoNormKrylovBasisSize );
    const Real twoNormEstG =
      TwoNormEstimate( problem.G, ctrl.twoNormKrylovBasisSize );
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
        const bool smallRegTooClose =
          ctrl.regIncreaseFactor*xRegSmall > xRegLarge ||
          ctrl.regIncreaseFactor*yRegSmall > yRegLarge ||
          ctrl.regIncreaseFactor*zRegSmall > zRegLarge;
        if( twoStage && smallRegTooClose )
        {
            twoStage = false;
            if( ctrl.print )
                Output("Disabling two-stage regularization");
        }
        if( twoStage )
        {
            if( ctrl.print )
                Output
                ("Increasing small regularization by a factor of ",
             ctrl.regIncreaseFactor);

            UpdateDiagonal( JStatic, ctrl.regIncreaseFactor-Real(1), regSmall );
            regSmall *= ctrl.regIncreaseFactor;
            xRegSmall *= ctrl.regIncreaseFactor;
            yRegSmall *= ctrl.regIncreaseFactor;
            zRegSmall *= ctrl.regIncreaseFactor;
        }
        else
        {
            if( ctrl.print )
                Output
                ("Increasing large regularization by a factor of ",
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
            // It seems that straight-forward equilibration in the two-stage
            // scheme can prevent the iterative solver from converging
            // (as small errors in the equilibrated scale can become enormous
            // in the original scale).
            if( !twoStage && ctrl.equilibrateIfSingleStage )
            {
                if( wDynamicRange >= ctrl.ruizEquilTol )
                    SymmetricRuizEquil
                    ( J, dInner, ctrl.ruizMaxIter, ctrl.print );
                else if( wDynamicRange >= ctrl.diagEquilTol )
                    SymmetricDiagonalEquil( J, dInner, ctrl.print );
                else
                    Ones( dInner, J.Height(), 1 );
            }
            else
            {
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
            return solveInfo.metRequestedTol;
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
            return solveInfo.metRequestedTol;
        }
      };

    AffineLPResidual<Matrix<Real>> residual, error;
    AffineLPSolution<Matrix<Real>> affineCorrection, correction;

    // matrix(row,:) *= scale
    auto scaleRow =
      [&]( Int row, const Real& scale, Matrix<Real>& matrix ) {
        matrix(row) *= scale;
      };
    auto scaleRowSparse =
      [&]( Int row, const Real& scale, SparseMatrix<Real>& matrix ) {
        Real* valBuf = matrix.ValueBuffer();
        const Int rowOffset = matrix.RowOffset(row);
        const Int nextRowOffset = matrix.RowOffset(row+1);
        for( Int index=rowOffset; index<nextRowOffset; ++index )
            valBuf[index] *= scale;
      };
    // matrix(scaleIndex,:) *= scale [except the diagonal]
    // matrix(:,scaleIndex) *= scale [except the diagonal]
    auto symmetricScaleG =
      [&]( const Matrix<Real>& zScaleNew, SparseMatrix<Real>& matrix ) {
        // | 0, A^T,  G^T |
        // | A,  0,    0  |
        // | G,  0,  -D^2 |
        const Int GHeight = zScaleNew.Height();
        const Int prevHeight = matrix.Height() - GHeight;
        Real* valBuf = matrix.ValueBuffer();
        const Int numEntries = matrix.NumEntries();
        for( Int index=0; index<numEntries; ++index )
        {
            const Int row = matrix.Row(index);
            const Int col = matrix.Col(index);
            if( row >= prevHeight && col >= prevHeight )
            {
                // We are in the nontrivial diagonal region.
            }
            else if( row >= prevHeight )
            {
                // We are in the row-prevHeight row of G
                valBuf[index] *= zScaleNew(row-prevHeight);
            }
            else if( col >= prevHeight )
            {
                // We are in the col-prevHeight column of G^T
                valBuf[index] *= zScaleNew(col-prevHeight);
            }
        }
      };

    // We will monotonically drive down the barrier parameter
    // (with the exception of the initial value).
    Real muMin = 1.e6;

    Real mu = 0.1, muOld = 0.1;
    const Int indent = PushIndent();
    for( ; numIts<=ctrl.maxIts; ++numIts, muOld=mu, dimacsErrorOld=dimacsError )
    {
        if( ctrl.dynamicallyRescale )
        {
            // Ensure that 0.01 <= || s ||_max <= 100
            const Real primalScaleNew = MaxNorm(solution.s);
            if( primalScaleNew < Real(0.1) || primalScaleNew > Real(10) )
            {
                equilibration.primalScale *= primalScaleNew;
                if( ctrl.print )
                    Output
                    ("s /= ",primalScaleNew," (total is ",
                     equilibration.primalScale,")");
                problem.b *= Real(1)/primalScaleNew;
                problem.h *= Real(1)/primalScaleNew;
                solution.x *= Real(1)/primalScaleNew;
                solution.s *= Real(1)/primalScaleNew;
            }

            // Ensure that 0.01 <= || z ||_max <= 100
            const Real dualScaleNew = MaxNorm(solution.z);
            if( dualScaleNew < Real(0.1) || dualScaleNew > Real(10) )
            {
                equilibration.dualScale *= dualScaleNew;
                if( ctrl.print )
                    Output
                    ("z /= ",dualScaleNew," (total is ",
                     equilibration.dualScale,")");
                problem.c *= Real(1)/dualScaleNew;
                solution.y *= Real(1)/dualScaleNew;
                solution.z *= Real(1)/dualScaleNew;
            }
        }
        const Real primalScale = MaxNorm(solution.s);
        const Real dualScale = MaxNorm(solution.z);

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
            Output("Complement ratio: ",compRatio);

        // Compute the regularized 'z' pivot.
        Matrix<Real> zPivot( solution.z );
        for( Int i=0; i<k; ++i )
            zPivot(i) = Max( zPivot(i), ctrl.zMinPivotValue );

        // Check for convergence
        // =====================

        // Carefully compute a relative duality gap.
        const Real dualProd =
          primalScale*dualScale*Dot( solution.s, solution.z );
        const Real primObj = PrimalObjective<Real>( problem, solution );
        const Real dualObj = DualObjective<Real>( problem, solution );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real relObjGap =
          RelativeObjectiveGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relCompGap, relObjGap );

        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Multiply
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
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
        const Real dualInfeasNrm2 = Nrm2( residual.dualEquality );
        const Real dualInfeasNrm2Rel = dualInfeasNrm2 / (1+cNrm2);

        // || G x + s - h ||_2 / (1 + || h ||_2)
        // -------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Multiply
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        residual.primalConic += solution.s;
        const Real conicInfeasNrm2 = Nrm2( residual.primalConic );
        const Real conicInfeasNrm2Rel = conicInfeasNrm2 / (1+hNrm2);

        // Now check the pieces
        // --------------------
        const Real equalityError = Max(primalInfeasNrm2Rel,dualInfeasNrm2Rel);
        infeasError = Max(conicInfeasNrm2Rel,equalityError);
        dimacsError = Max(maxRelGap,infeasError);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( solution.x );
            const Real yNrm2 = Nrm2( solution.y );
            const Real zNrm2 = Nrm2( solution.z );
            const Real sNrm2 = Nrm2( solution.s );
            Output
            ("iter ",numIts,":\n",Indent(),
             "  || x~ ||_2 = ",xNrm2,"\n",Indent(),
             "  || y~ ||_2 = ",yNrm2,"\n",Indent(),
             "  || z~ ||_2 = ",zNrm2,"\n",Indent(),
             "  || s~ ||_2 = ",sNrm2,"\n",Indent(),
             "  || primalInfeas ||_2 / (1 + || b ||_2) = ",
             primalInfeasNrm2Rel,"\n",Indent(),
             "  || dualInfeas   ||_2 / (1 + || c ||_2) = ",
             dualInfeasNrm2Rel,"\n",Indent(),
             "  || conicInfeas  ||_2 / (1 + || h ||_2) = ",
             conicInfeasNrm2Rel,"\n",Indent(),
             "  scaled primal = ",primObj,"\n",Indent(),
             "  scaled dual   = ",dualObj,"\n",Indent(),
             "  scaled rel obj gap = ",relObjGap,"\n",Indent(),
             "  scaled rel comp gap = ",relCompGap,"\n",Indent(),
             "  scaled dimacs error = ",dimacsError);

            AffineLPSolution<Matrix<Real>> origSolution;
            UndoEquilibration( solution, equilibration, origSolution );

            const Real dualProdOrig = Dot( origSolution.s, origSolution.z );
            const Real primObjOrig =
              PrimalObjective<Real>( origProblem, origSolution );
            const Real dualObjOrig =
              DualObjective<Real>( origProblem, origSolution );
            const Real relCompGapOrig =
              RelativeComplementarityGap
              ( primObjOrig, dualObjOrig, dualProdOrig );
            const Real relObjGapOrig =
              RelativeObjectiveGap( primObjOrig, dualObjOrig, dualProdOrig );
            const Real maxRelGapOrig = Max( relCompGapOrig, relObjGapOrig );

            Output("  s^T z = ",dualProdOrig);
            Output("  primal: ",primObjOrig);
            Output("  dual: ",dualObjOrig);
            Output("  relative gap: ",maxRelGapOrig);
        }

        const bool metTolerances =
          infeasError <= ctrl.infeasibilityTol &&
          relCompGap <= ctrl.relativeComplementarityGapTol &&
          relObjGap <= ctrl.relativeObjectiveGapTol;
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

        // Factor the KKT system
        // =====================
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        FinishKKT( m, n, solution.s, zPivot, JOrig );
        J = JOrig;
        J.FreezeSparsity();
        UpdateDiagonal( J, Real(1), regLarge );
        if( !attemptToFactor(wDynamicRange) )
        {
            increaseRegularization();
            continue;
        }

        // TODO(poulson): Pin the variables with corresponding sufficiently
        // large diagonal entries in 'J' to zero via modifying the right-hand
        // side and 'J'.

        // Compute the affine search direction
        // ===================================
        const bool largeCompRatio = compRatio > ctrl.maxComplementRatio;
        const bool freezeBarrier = largeCompRatio || maxRelGap < equalityError;
        const Real muClassical = Dot(solution.s,solution.z) / k;
        if( freezeBarrier )
        {
            mu = Max( muOld, muClassical );
            if( ctrl.print )
            {
                Output("Mu classical was ",muClassical);
                Output("Freezing barrier at ",mu);
            }
        }
        else
        {
            mu = Min( muClassical, muMin );
            muMin = Min( mu, muMin );
            if( ctrl.print )
                Output("New barrier parameter is ",mu);
        }

        // r_mu := s o zPivot
        // ------------------
        residual.dualConic = zPivot;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

        // Solve for the proposed step
        // ---------------------------
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          zPivot, d );
        if( !attemptToSolve(d) )
        {
            increaseRegularization();
            continue;
        }
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, zPivot,
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
            Axpy( xRegSmall, affineCorrection.x, error.dualEquality );
            if( !twoStage )
                Axpy( xRegLarge, affineCorrection.x, error.dualEquality );
            Multiply
            ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
              Real(1), error.dualEquality );
            Multiply
            ( TRANSPOSE, Real(1), problem.G, affineCorrection.z,
              Real(1), error.dualEquality );
            const Real dyErrorNrm2 = Nrm2( error.dualEquality );

            error.primalConic = residual.primalConic;
            error.primalConic += affineCorrection.s;
            Multiply
            ( NORMAL, Real(1), problem.G, affineCorrection.x,
              Real(1), error.primalConic );
            Axpy( -zRegSmall, affineCorrection.z, error.primalConic );
            const Real dzErrorNrm2 = Nrm2( error.primalConic );

            // TODO(poulson): error.dualConic

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+primalInfeasNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+dualInfeasNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+conicInfeasNrm2));
        }

        // Handle updates which push s or z close too close to zero
        // ========================================================
        const Real sAffineTwo = FrobeniusNorm( affineCorrection.s );
        const Real zAffineTwo = FrobeniusNorm( affineCorrection.z );
        Output("|| sAffine ||_2 = ",sAffineTwo);
        Output("|| zAffine ||_2 = ",zAffineTwo);
        Real sMin = MaxNorm( solution.s );
        Real zMin = MaxNorm( solution.z );
        for( Int i=0; i<n; ++i )
        {
            sMin = Min( sMin, solution.s(i) );
            zMin = Min( zMin, solution.z(i) );
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
        if( largeCompRatio && ctrl.softDualityTargets )
        {
            // Attempt to correct s o z entrywise into
            // [lowerTargetRatio,upperTargetRatio]*sigma*mu.
            Real lowerTargetRatio =
              Pow(compRatio,ctrl.lowerTargetRatioLogCompRatio);
            Real upperTargetRatio =
              Pow(compRatio,ctrl.upperTargetRatioLogCompRatio);
            Output
            ("compRatio=",compRatio,", lowerTargetRatio=",lowerTargetRatio,
             ", upperTargetRatio=",upperTargetRatio);
            lowerTargetRatio =
              Max( lowerTargetRatio, Pow(ctrl.maxComplementRatio,Real(-0.5)) );
            upperTargetRatio =
              Min( upperTargetRatio, Pow(ctrl.maxComplementRatio,Real(0.5)) );
            Output
            ("lowerTargetRatio=",lowerTargetRatio,
             ", upperTargetRatio=",upperTargetRatio);
            const Real lowerTarget = lowerTargetRatio*sigma*mu;
            const Real upperTarget = upperTargetRatio*sigma*mu;
            for( Int i=0; i<k; ++i )
            {
                const Real prod = solution.s(i)*zPivot(i);
                if( prod < lowerTarget )
                {
                    Output(i,": ",prod," was < lowerTarget=",lowerTarget);
                    residual.dualConic(i) -= lowerTarget;
                }
                else if( prod > upperTarget )
                {
                    Output(i,": ",prod," was > upperTarget=",upperTarget);
                    residual.dualConic(i) -= upperTarget;
                }
                else
                {
                    residual.dualConic(i) = 0;
                }
            }
        }
        else
        {
            Shift( residual.dualConic, -sigma*mu );
            if( ctrl.mehrotra )
            {
                // r_mu := dsAff o dzAff
                // ---------------------
                // NOTE: dz is used as a temporary
                correction.z = affineCorrection.z;
                DiagonalScale
                ( LEFT, NORMAL, affineCorrection.s, correction.z );
                residual.dualConic += correction.z;
            }
        }

        // Construct the new KKT RHS
        // -------------------------
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          zPivot, d );
        // Solve for the proposed step
        // ---------------------------
        if( !attemptToSolve(d) )
        {
            increaseRegularization();
            continue;
        }
        ExpandSolution
        ( m, n, d, residual.dualConic, solution.s, zPivot,
          correction.x, correction.y, correction.z, correction.s );

        const Real sCombinedTwo = FrobeniusNorm( correction.s );
        const Real zCombinedTwo = FrobeniusNorm( correction.z );
        Output("|| sCombined ||_2 = ",sCombinedTwo);
        Output("|| zCombined ||_2 = ",zCombinedTwo);

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
    if( ctrl.print )
    {
        Output
        ("Finished with:\n",Indent(),
         "  twoStage: ",twoStage,"\n",Indent(),
         "  xRegSmall: ",xRegSmall,"\n",Indent(),
         "  yRegSmall: ",yRegSmall,"\n",Indent(),
         "  zRegSmall: ",zRegSmall,"\n",Indent(),
         "  xRegLarge: ",xRegLarge,"\n",Indent(),
         "  yRegLarge: ",yRegLarge,"\n",Indent(),
         "  zRegLarge: ",zRegLarge,"\n",Indent());
    }
}

template<typename Real>
void RecenterAffineLPProblem
( const AffineLPSolution<Matrix<Real>>& center,
        AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem )
{
    EL_DEBUG_CSE
    // Recenter about (x,0;y,z), purposefully avoiding s. We put
    //   bHat := b - A x,
    //   hHat := h - G x,
    //   cHat := c + A^T y + G^T z.
    Multiply( NORMAL, Real(-1), problem.A, center.x, Real(1), problem.b );
    Multiply( NORMAL, Real(-1), problem.G, center.x, Real(1), problem.h );
    Multiply( TRANSPOSE, Real(1), problem.A, center.y, Real(1), problem.c );
    Multiply( TRANSPOSE, Real(1), problem.G, center.z, Real(1), problem.c );
}

template<typename Real>
void IPM
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
        EquilibratedIPM
        ( problem, equilibration,
          equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );

        // This seems to dramatically fail on several examples (e.g., dfl001).
        const bool tryRefinement = false;
        if( tryRefinement )
        {
            // Perform one refinement step using the modified RHS's.
            auto refinedProblem = equilibratedProblem;
            RecenterAffineLPProblem( equilibratedSolution, refinedProblem );

            // Compute the primal refinement scale.
            const Real bHatMax = MaxNorm( refinedProblem.b );
            const Real hHatMax = MaxNorm( refinedProblem.h );
            const Real primalRefineScale = Real(1)/Max(bHatMax,hHatMax);
            // Compute the dual refinement scale.
            const Real cHatMax = MaxNorm( refinedProblem.c );
            const Real dualRefineScale = Real(1)/cHatMax;
            Output("|| bHat ||_max = ",bHatMax);
            Output("|| hHat ||_max = ",hHatMax);
            Output("|| cHat ||_max = ",cHatMax);
            Output
            ("primalScale=",primalRefineScale,", dualScale=",dualRefineScale);
            refinedProblem.b *= primalRefineScale;
            refinedProblem.h *= primalRefineScale;
            refinedProblem.c *= dualRefineScale;

            AffineLPSolution<Matrix<Real>> refinedSolution,
              equilibratedRefinedSolution;
            SparseAffineLPEquilibration<Real> trivialEquilibration;
            trivialEquilibration.primalScale = 1;
            trivialEquilibration.dualScale = 1;
            Ones( trivialEquilibration.xScale, problem.A.Width(), 1 );
            Ones( trivialEquilibration.yScale, problem.A.Height(), 1 );
            Ones( trivialEquilibration.zScale, problem.G.Height(), 1 );
            auto refinedProblemCopy( refinedProblem );
            EquilibratedIPM
            ( refinedProblemCopy, trivialEquilibration,
              refinedProblem, equilibratedRefinedSolution, ctrl );
            UndoEquilibration
            ( equilibratedRefinedSolution, trivialEquilibration,
              refinedSolution );

            Axpy( Real(1)/primalRefineScale, refinedSolution.x,
              equilibratedSolution.x );
            Axpy( Real(1)/primalRefineScale, refinedSolution.s,
              equilibratedSolution.s );
            Axpy( Real(1)/dualRefineScale, refinedSolution.y,
              equilibratedSolution.y );
            Axpy( Real(1)/dualRefineScale, refinedSolution.z,
              equilibratedSolution.z );
            UndoEquilibration( equilibratedSolution, equilibration, solution );

            const Real newPrimal = Dot( problem.c, solution.x );
            const Real newDual = -Dot( problem.b, solution.y ) -
              Dot( problem.h, solution.z );
            Output("New primal: ",newPrimal,", new dual: ",newDual);
        }
    }
    else
    {
        SparseAffineLPEquilibration<Real> equilibration;
        equilibration.primalScale = 1;
        equilibration.dualScale = 1;
        Ones( equilibration.xScale, problem.A.Width(), 1 );
        Ones( equilibration.yScale, problem.A.Height(), 1 );
        Ones( equilibration.zScale, problem.G.Height(), 1 );
        auto equilibratedProblem = problem;
        auto equilibratedSolution = solution;
        EquilibratedIPM
        ( problem, equilibration,
          equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
}

template<typename Real>
void EquilibratedIPM
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Int degree = k;
    const Grid& grid = problem.A.Grid();
    const int commRank = grid.Rank();
    Timer timer;

    if( ctrl.regIncreaseFactor <= Real(1) )
        LogicError("Regularization increase factor must be at least 1");

    const Real bNrm2 = Nrm2( problem.b );
    const Real cNrm2 = Nrm2( problem.c );
    const Real hNrm2 = Nrm2( problem.h );
    const Real twoNormEstA =
      TwoNormEstimate( problem.A, ctrl.twoNormKrylovBasisSize );
    const Real twoNormEstG =
      TwoNormEstimate( problem.G, ctrl.twoNormKrylovBasisSize );
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
        return solveInfo.metRequestedTol;
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
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, dualProd );
        const Real relObjGap =
          RelativeObjectiveGap( primObj, dualObj, dualProd );
        const Real maxRelGap = Max( relCompGap, relObjGap );
        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        residual.primalEquality *= -1;
        Multiply
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(1), residual.primalEquality );
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
        const Real dualInfeasNrm2 = Nrm2( residual.dualEquality );
        const Real dualInfeasNrm2Rel = dualInfeasNrm2 / (1+cNrm2);

        // || G x + s - h ||_2 / (1 + || h ||_2)
        // -------------------------------------
        residual.primalConic = problem.h;
        residual.primalConic *= -1;
        Multiply
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(1), residual.primalConic );
        residual.primalConic += solution.s;
        const Real conicInfeasNrm2 = Nrm2( residual.primalConic );
        const Real conicInfeasNrm2Rel = conicInfeasNrm2 / (1+hNrm2);

        // Now check the pieces
        // --------------------
        const Real equalityError = Max(primalInfeasNrm2Rel,dualInfeasNrm2Rel);
        infeasError = Max(equalityError,conicInfeasNrm2Rel);
        dimacsError = Max(infeasError,maxRelGap);
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
                 "  scaled relative gap = ",maxRelGap);
        }

        const bool metTolerances =
          infeasError <= ctrl.infeasibilityTol &&
          relCompGap <= ctrl.relativeComplementarityGapTol &&
          relObjGap <= ctrl.relativeObjectiveGapTol;
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

        // Factor the KKT system
        // =====================
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        JOrig.LockedDistGraph().multMeta = JStatic.LockedDistGraph().multMeta;
        FinishKKT( m, n, solution.s, solution.z, JOrig );
        J = JOrig;
        J.FreezeSparsity();
        J.LockedDistGraph().multMeta = JStatic.LockedDistGraph().multMeta;
        UpdateDiagonal( J, Real(1), regLarge );
        if( !attemptToFactor(wMaxNorm) )
        {
            // TODO(poulson): Increase regularization and continue.
            break;
        }

        // Compute the affine search direction
        // ===================================

        // r_mu := s o z
        // -------------
        residual.dualConic = solution.z;
        DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

        // Solve for the proposed step
        // ---------------------------
        KKTRHS
        ( residual.dualEquality,
          residual.primalEquality,
          residual.primalConic,
          residual.dualConic,
          solution.z, d );
        if( !attemptToSolve(d) )
        {
            // TODO(poulson): Increase regularization and continue.
            break;
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
        {
            // TODO(poulson): Increase regularization and continue.
            break;
        }
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
( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem,
        AffineLPSolution<DistMultiVec<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
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
        EquilibratedIPM( equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    else
    {
        EquilibratedIPM( problem, solution, ctrl );
    }
}

#define PROTO(Real) \
  template void IPM \
  ( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem, \
          AffineLPSolution<Matrix<Real>>& solution, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const AffineLPProblem<DistMatrix<Real>,DistMatrix<Real>>& problem, \
          AffineLPSolution<DistMatrix<Real>>& solution, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem, \
          AffineLPSolution<Matrix<Real>>& solution, \
    const IPMCtrl<Real>& ctrl ); \
  template void IPM \
  ( const AffineLPProblem<DistSparseMatrix<Real>,DistMultiVec<Real>>& problem, \
          AffineLPSolution<DistMultiVec<Real>>& solution, \
    const IPMCtrl<Real>& ctrl );

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
