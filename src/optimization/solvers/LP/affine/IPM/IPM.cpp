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
//                - mu Phi(s),
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
//     - mu Phi(D_s s').
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
//     - mu Phi(inv(D_z) s').
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
//     - mu Phi(primalScale inv(D_z) s').
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
//     - mu Phi(primalScale inv(D_z) s').
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
struct AffineLPState
{
    Int numIts=0;

    Real barrier=Real(0.1);
    Real barrierOld=Real(0.1);

    Real xRegSmall;
    Real yRegSmall;
    Real zRegSmall;
    Real xRegLarge;
    Real yRegLarge;
    Real zRegLarge;

    Real cNrm2;
    Real bNrm2;
    Real hNrm2;

    Real primalEqualityInfeasNrm2=Real(1);
    Real primalEqualityInfeasNrm2Rel=Real(1);

    Real primalConicInfeasNrm2=Real(1);
    Real primalConicInfeasNrm2Rel=Real(1);

    Real dualEqualityInfeasNrm2=Real(1);
    Real dualEqualityInfeasNrm2Rel=Real(1);

    // The maximum of the primal and dual equality relative infeas. measures.
    Real equalityError=Real(1);

    // The maximum of the above relative infeasibility measures.
    Real infeasError=Real(1);

    // c^T x
    Real primalObj;
    Real primalObjOrig;

    // -b^T y - h^T z
    Real dualObj;
    Real dualObjOrig;

    // s^T z
    Real gap=Real(1);
    Real gapOrig=Real(1);

    // complementRatio := (max_i s_i * z_i) / (min_i s_i * z_i).
    Real complementRatio=Real(1);

    Real relObjGap=Real(1);
    Real relObjGapOrig=Real(1);

    Real relCompGap=Real(1);
    Real relCompGapOrig=Real(1);

    Real maxRelGap=Real(1);
    Real maxRelGapOrig=Real(1);

    Real dimacsError=Real(1);
    Real dimacsErrorOrig=Real(1);

    Real dimacsErrorOld=Real(1);
    Real dimacsErrorOrigOld=Real(1);

    bool metTolerances=false;
    bool metTolerancesOrig=false;

    bool backingOffEquilibration=false;
};

template<typename Real>
void RescalePrimal
( const Real& primalScale,
  DenseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  AffineLPSolution<Matrix<Real>>& solution,
  bool print )
{
    EL_DEBUG_CSE
    equilibration.primalScale *= primalScale;
    if( print )
        Output
        ("primal /= ",primalScale," (total is ",
         equilibration.primalScale,")");
    problem.b *= Real(1)/primalScale;
    problem.h *= Real(1)/primalScale;
    solution.x *= Real(1)/primalScale;
    solution.s *= Real(1)/primalScale;
}

template<typename Real>
void RescalePrimal
( const Real& primalScale,
  SparseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  AffineLPSolution<Matrix<Real>>& solution,
  bool print )
{
    EL_DEBUG_CSE
    equilibration.primalScale *= primalScale;
    if( print )
        Output
        ("primal /= ",primalScale," (total is ",
         equilibration.primalScale,")");
    problem.b *= Real(1)/primalScale;
    problem.h *= Real(1)/primalScale;
    solution.x *= Real(1)/primalScale;
    solution.s *= Real(1)/primalScale;
}

template<typename Real>
void RescaleDual
( const Real& dualScale,
  DenseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  AffineLPSolution<Matrix<Real>>& solution,
  bool print )
{
    equilibration.dualScale *= dualScale;
    if( print )
        Output
        ("dual /= ",dualScale," (total is ",equilibration.dualScale,")");
    problem.c *= Real(1)/dualScale;
    solution.y *= Real(1)/dualScale;
    solution.z *= Real(1)/dualScale;
}

template<typename Real>
void RescaleDual
( const Real& dualScale,
  SparseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  AffineLPSolution<Matrix<Real>>& solution,
  bool print )
{
    equilibration.dualScale *= dualScale;
    if( print )
        Output
        ("dual /= ",dualScale," (total is ",equilibration.dualScale,")");
    problem.c *= Real(1)/dualScale;
    solution.y *= Real(1)/dualScale;
    solution.z *= Real(1)/dualScale;
}

template<typename Real>
void RescaleKKTButNotRegularization
( const Matrix<Real>& xScale,
  const Matrix<Real>& yScale,
  const Matrix<Real>& zScale,
        Matrix<Real>& J )
{
    EL_DEBUG_CSE
    const Int n = xScale.Height();
    const Int m = yScale.Height();
    const Int k = zScale.Height();
    Real* JBuf = J.Buffer();
    const Int JLDim = J.LDim();
    for( Int col=0; col<n+m+k; ++col )
    {
        for( Int row=0; row<n+m+k; ++row )
        {
            Real& value = JBuf[row+col*JLDim];
            // | 0 A^T G^T |
            // | A  0   0  |
            // | G  0   0  |
            if( row < n )
            {
                if( col >= n && col < n+m )
                {
                    // A^T block
                    value *= yScale(col-n) * xScale(row);
                }
                else if( col >= n+m )
                {
                    // G^T block
                    value *= zScale(col-(n+m)) * xScale(row);
                }
            }
            else if( row < n+m )
            {
                if( col < n )
                {
                    // A block
                    value *= xScale(col) * yScale(row-n);
                }
            }
            else
            {
                if( col < n )
                {
                    // G block
                    value *= xScale(col) * zScale(row-(n+m));
                }
            }
        }
    }
}

template<typename Real>
void RescaleKKTButNotRegularization
( const Matrix<Real>& xScale,
  const Matrix<Real>& yScale,
  const Matrix<Real>& zScale,
        SparseMatrix<Real>& J )
{
    EL_DEBUG_CSE
    const Int n = xScale.Height();
    const Int m = yScale.Height();
    Real* valBuf = J.ValueBuffer();
    const Int numEntries = J.NumEntries();
    for( Int index=0; index<numEntries; ++index )
    {
        const Int row = J.Row(index);
        const Int col = J.Col(index);
        Real& value = valBuf[index];
        // | 0 A^T G^T |
        // | A  0   0  |
        // | G  0   0  |
        if( row < n )
        {
            if( col >= n && col < n+m )
            {
                // A^T block
                value *= yScale(col-n) * xScale(row);
            }
            else if( col >= n+m )
            {
                // G^T block
                value *= zScale(col-(n+m)) * xScale(row);
            }
        }
        else if( row < n+m )
        {
            if( col < n )
            {
                // A block
                value *= xScale(col) * yScale(row-n);
            }
        }
        else
        {
            if( col < n )
            {
                // G block
                value *= xScale(col) * zScale(row-(n+m));
            }
        }
    }
}

template<typename Real>
void BackOffEquilibration
( DenseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  AffineLPState<Real>& state,
  AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real xScaleMin = Min( equilibration.xScale );
    const Real yScaleMin = Min( equilibration.yScale );
    const Real zScaleMin = Min( equilibration.zScale );
    const Real xScaleMax = MaxNorm( equilibration.xScale );
    const Real yScaleMax = MaxNorm( equilibration.yScale );
    const Real zScaleMax = MaxNorm( equilibration.zScale );
    if( xScaleMin == xScaleMax &&
        yScaleMin == yScaleMax &&
        zScaleMin == zScaleMax )
        return;
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();

    if( ctrl.print )
    {
        Output
        ("Backing off with "
         "  primalScale = ",equilibration.primalScale,
         ", dualScale = ",equilibration.dualScale,
         "  || xScale ||_min = ",xScaleMin,
         ", || yScale ||_min = ",yScaleMin,
         ", || zScale ||_min = ",zScaleMin,
         "  || xScale ||_max = ",xScaleMax,
         ", || yScale ||_max = ",yScaleMax,
         ", || zScale ||_max = ",zScaleMax);
    }

    // Force || xScale ||_max =
    //       max(|| yScale ||_max,|| zScale ||_max) = 1
    // by absorbing into primalScale and dualScale.
    const Real primalRescaling = xScaleMax;
    const Real dualRescaling = Max(yScaleMax,zScaleMax);
    equilibration.xScale *= Real(1)/primalRescaling;
    equilibration.yScale *= Real(1)/dualRescaling;
    equilibration.zScale *= Real(1)/dualRescaling;
    equilibration.primalScale *= primalRescaling;
    equilibration.dualScale *= dualRescaling;

    Matrix<Real> xScaleMultiple, yScaleMultiple, zScaleMultiple;
    Ones( xScaleMultiple, n, 1 );
    Ones( yScaleMultiple, m, 1 );
    Ones( zScaleMultiple, k, 1 );
    for( Int i=0; i<n; ++i )
        xScaleMultiple(i) =
          Pow( equilibration.xScale(i), ctrl.backoffScalePower );
    for( Int i=0; i<m; ++i )
        yScaleMultiple(i) =
          Pow( equilibration.yScale(i), ctrl.backoffScalePower );
    for( Int i=0; i<k; ++i )
        zScaleMultiple(i) =
          Pow( equilibration.zScale(i), ctrl.backoffScalePower );

    // Handle x rescaling (except for first-order optimality matrix)
    DiagonalSolve( LEFT, NORMAL, xScaleMultiple, solution.x );
    DiagonalScale
    ( LEFT, NORMAL, xScaleMultiple, equilibration.xScale );
    DiagonalScale( LEFT, NORMAL, xScaleMultiple, problem.c );
    state.cNrm2 = FrobeniusNorm( problem.c );
    DiagonalScale( RIGHT, NORMAL, xScaleMultiple, problem.A );
    DiagonalScale( RIGHT, NORMAL, xScaleMultiple, problem.G );

    // Handle y rescaling (except for first-order optimality matrix)
    DiagonalSolve( LEFT, NORMAL, yScaleMultiple, solution.y );
    DiagonalScale
    ( LEFT, NORMAL, yScaleMultiple, equilibration.yScale );
    DiagonalScale( LEFT, NORMAL, yScaleMultiple, problem.A );
    DiagonalScale( LEFT, NORMAL, yScaleMultiple, problem.b );
    state.bNrm2 = FrobeniusNorm( problem.b );

    // Handle z rescaling (except for first-order optimality matrix)
    DiagonalScale( LEFT, NORMAL, zScaleMultiple, solution.s );
    DiagonalSolve( LEFT, NORMAL, zScaleMultiple, solution.z );
    DiagonalScale
    ( LEFT, NORMAL, zScaleMultiple, equilibration.zScale );
    DiagonalScale( LEFT, NORMAL, zScaleMultiple, problem.G );
    DiagonalScale( LEFT, NORMAL, zScaleMultiple, problem.h );
    state.hNrm2 = FrobeniusNorm( problem.h );
}

template<typename Real>
void BackOffEquilibration
( SparseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  AffineLPState<Real>& state,
  AffineLPSolution<Matrix<Real>>& solution,
  SparseMatrix<Real>& JStatic,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real xScaleMin = Min( equilibration.xScale );
    const Real yScaleMin = Min( equilibration.yScale );
    const Real zScaleMin = Min( equilibration.zScale );
    const Real xScaleMax = MaxNorm( equilibration.xScale );
    const Real yScaleMax = MaxNorm( equilibration.yScale );
    const Real zScaleMax = MaxNorm( equilibration.zScale );
    if( xScaleMin == xScaleMax &&
        yScaleMin == yScaleMax &&
        zScaleMin == zScaleMax )
        return;
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();

    if( ctrl.print )
    {
        Output
        ("Backing off with "
         "  primalScale = ",equilibration.primalScale,
         ", dualScale = ",equilibration.dualScale,
         "  || xScale ||_min = ",xScaleMin,
         ", || yScale ||_min = ",yScaleMin,
         ", || zScale ||_min = ",zScaleMin,
         "  || xScale ||_max = ",xScaleMax,
         ", || yScale ||_max = ",yScaleMax,
         ", || zScale ||_max = ",zScaleMax);
    }

    // Force || xScale ||_max =
    //       max(|| yScale ||_max,|| zScale ||_max) = 1
    // by absorbing into primalScale and dualScale.
    const Real primalRescaling = xScaleMax;
    const Real dualRescaling = Max(yScaleMax,zScaleMax);
    equilibration.xScale *= Real(1)/primalRescaling;
    equilibration.yScale *= Real(1)/dualRescaling;
    equilibration.zScale *= Real(1)/dualRescaling;
    equilibration.primalScale *= primalRescaling;
    equilibration.dualScale *= dualRescaling;

    Matrix<Real> xScaleMultiple, yScaleMultiple, zScaleMultiple;
    Ones( xScaleMultiple, n, 1 );
    Ones( yScaleMultiple, m, 1 );
    Ones( zScaleMultiple, k, 1 );
    for( Int i=0; i<n; ++i )
        xScaleMultiple(i) =
          Pow( equilibration.xScale(i), ctrl.backoffScalePower );
    for( Int i=0; i<m; ++i )
        yScaleMultiple(i) =
          Pow( equilibration.yScale(i), ctrl.backoffScalePower );
    for( Int i=0; i<k; ++i )
        zScaleMultiple(i) =
          Pow( equilibration.zScale(i), ctrl.backoffScalePower );

    // Handle x rescaling (except for first-order optimality matrix)
    DiagonalSolve( LEFT, NORMAL, xScaleMultiple, solution.x );
    DiagonalScale
    ( LEFT, NORMAL, xScaleMultiple, equilibration.xScale );
    DiagonalScale( LEFT, NORMAL, xScaleMultiple, problem.c );
    state.cNrm2 = FrobeniusNorm( problem.c );
    DiagonalScale( RIGHT, NORMAL, xScaleMultiple, problem.A );
    DiagonalScale( RIGHT, NORMAL, xScaleMultiple, problem.G );

    // Handle y rescaling (except for first-order optimality matrix)
    DiagonalSolve( LEFT, NORMAL, yScaleMultiple, solution.y );
    DiagonalScale
    ( LEFT, NORMAL, yScaleMultiple, equilibration.yScale );
    DiagonalScale( LEFT, NORMAL, yScaleMultiple, problem.A );
    DiagonalScale( LEFT, NORMAL, yScaleMultiple, problem.b );
    state.bNrm2 = FrobeniusNorm( problem.b );

    // Handle z rescaling (except for first-order optimality matrix)
    DiagonalScale( LEFT, NORMAL, zScaleMultiple, solution.s );
    DiagonalSolve( LEFT, NORMAL, zScaleMultiple, solution.z );
    DiagonalScale
    ( LEFT, NORMAL, zScaleMultiple, equilibration.zScale );
    DiagonalScale( LEFT, NORMAL, zScaleMultiple, problem.G );
    DiagonalScale( LEFT, NORMAL, zScaleMultiple, problem.h );
    state.hNrm2 = FrobeniusNorm( problem.h );

    RescaleKKTButNotRegularization
    ( xScaleMultiple, yScaleMultiple, zScaleMultiple, JStatic );
}

template<typename Real>
void NeutralizePrimalAndDualNorms
( DenseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  AffineLPState<Real>& state,
  AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    // Ensure that max( || x ||_2, || s ||_2 ) ~= 1.
    const Real primalNorm =
      Max( FrobeniusNorm(solution.x), FrobeniusNorm(solution.s) );
    if( primalNorm < ctrl.primalNormLowerBound ||
        primalNorm > ctrl.primalNormUpperBound )
    {
        RescalePrimal
        ( primalNorm, equilibration, problem, solution, ctrl.print );
        state.bNrm2 *= Real(1)/primalNorm;
        state.hNrm2 *= Real(1)/primalNorm;
    }

    // Ensure that max( || y ||_2, || z ||_2 ) ~= 1.
    const Real dualNorm =
      Max( FrobeniusNorm(solution.y), FrobeniusNorm(solution.z) );
    if( dualNorm < ctrl.dualNormLowerBound ||
        dualNorm > ctrl.dualNormUpperBound )
    {
        RescaleDual
        ( dualNorm, equilibration, problem, solution, ctrl.print );
        state.cNrm2 *= Real(1)/dualNorm;
    }
}

template<typename Real>
void NeutralizePrimalAndDualNorms
( SparseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  AffineLPState<Real>& state,
  AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE

    // Ensure that max( || x ||_2, || s ||_2 ) ~= 1.
    const Real primalNorm =
      Max( FrobeniusNorm(solution.x), FrobeniusNorm(solution.s) );
    if( primalNorm < ctrl.primalNormLowerBound ||
        primalNorm > ctrl.primalNormUpperBound )
    {
        RescalePrimal
        ( primalNorm, equilibration, problem, solution, ctrl.print );
        state.bNrm2 *= Real(1)/primalNorm;
        state.hNrm2 *= Real(1)/primalNorm;
    }

    // Ensure that max( || y ||_2, || z ||_2 ) ~= 1.
    const Real dualNorm =
      Max( FrobeniusNorm(solution.y), FrobeniusNorm(solution.z) );
    if( dualNorm < ctrl.dualNormLowerBound ||
        dualNorm > ctrl.dualNormUpperBound )
    {
        RescaleDual
        ( dualNorm, equilibration, problem, solution, ctrl.print );
        state.cNrm2 *= Real(1)/dualNorm;
    }
}

template<typename Real>
void AssertPositiveOrthantMembership
( const AffineLPSolution<Matrix<Real>>& solution )
{
    EL_DEBUG_CSE
    const Int sNumNonPos = pos_orth::NumOutside( solution.s );
    const Int zNumNonPos = pos_orth::NumOutside( solution.z );
    if( sNumNonPos > 0 || zNumNonPos > 0 )
        LogicError
        (sNumNonPos," entries of s were nonpositive and ",
         zNumNonPos," entries of z were nonpositive");
}

template<typename Real>
void UpdateCombinedResidualUsingAffine
( const Real& sigma,
  const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const AffineLPState<Real>& state,
  const AffineLPSolution<Matrix<Real>>& solution,
  const AffineLPSolution<Matrix<Real>>& affineCorrection,
        AffineLPResidual<Matrix<Real>>& residual,
  bool compositeNewton )
{
    EL_DEBUG_CSE
    Shift( residual.dualConic, -sigma*state.barrier );
    // Use the interpretation of Section 2 of [CITATION]
    //
    //   R.A. Tapia, Y. Zhang, M. Saltzman, and A. Weiser's
    //   "The Mehrotra Predictor-Corrector Interior-Point Method as a
    //   Perturbed Composite Newton Method", July, 1990
    //   http://www.caam.rice.edu/tech_reports/1990/TR90-17.pdf
    //
    // as a means of generalizing the classic Mehrotra predictor-corrector
    // to approximate feasibility by updating with
    // F(solution+affineCorrection) rather than (0,0,0,dsAff o dzAff).
    //
    // With that said, it appears that neither the original or generalized
    // algorithm is generally helpful on the netlib/LP_data examples
    // (at least, on top of the current regularization and dynamic scaling).
    //
    const bool infeasibleCompositeNewton = true;
    if( compositeNewton && infeasibleCompositeNewton )
    {
        Matrix<Real> input, update;

        // residual.primalEquality +=
        //   A (x + dxAff) - b - (yRegSmall + yRegLarge) dyAff
        Zeros( update, problem.A.Height(), 1 );
        input = solution.x;
        input += affineCorrection.x;
        Gemv( NORMAL, Real(1), problem.A, input, Real(1), update );
        Axpy( Real(-1), problem.b, update );
        Axpy( -state.yRegSmall-state.yRegLarge, affineCorrection.y, update );
        const Real primalEqualityUpdateNrm2 = FrobeniusNorm( update );
        residual.primalEquality += update;

        // residual.dualEquality +=
        //   c + A^T (y + dyAff) + G^T (z + dzAff) +
        //   (xRegSmall + xRegLarge) dxAff
        Zeros( update, problem.A.Width(), 1 );
        input += problem.c;
        input = solution.y;
        input += affineCorrection.y;
        Gemv( TRANSPOSE, Real(1), problem.A, input, Real(1), update );
        input = solution.z;
        input += affineCorrection.z;
        Gemv( TRANSPOSE, Real(1), problem.G, input, Real(1), update );
        Axpy( state.xRegSmall+state.xRegLarge, affineCorrection.x, update );
        const Real dualEqualityUpdateNrm2 = FrobeniusNorm( update );
        residual.dualEquality += update;

        // residual.primalConic +=
        //   G (x + dxAff) + (s + dsAff) - h - (zRegSmall + zRegLarge) dzAff
        Zeros( update, problem.G.Height(), 1 );
        input = solution.x;
        input += affineCorrection.x;
        Gemv( NORMAL, Real(1), problem.G, input, Real(1), update );
        input = solution.s;
        input += affineCorrection.s;
        update += input;
        Axpy( Real(-1), problem.h, update );
        Axpy( -state.zRegSmall-state.zRegLarge, affineCorrection.z, update );
        const Real primalConicUpdateNrm2 = FrobeniusNorm( update );
        residual.primalConic += update;

        // residual.dualConic += (s + dsAff) o (z + dzAff)
        // TODO(poulson): Incorporate zPivot instead of z?
        input = solution.s;
        input += affineCorrection.s;
        update = solution.z;
        update += affineCorrection.z;
        DiagonalScale( LEFT, NORMAL, input, update );
        const Real dualConicUpdateNrm2 = FrobeniusNorm( update );
        residual.dualConic += update;
    }
    else if( compositeNewton )
    {
        // Use the approximation of the composite Newton method that assumes
        // the update would provide strict feasibility.

        // r_mu := dsAff o dzAff
        // ---------------------
        Matrix<Real> update( affineCorrection.z );
        DiagonalScale( LEFT, NORMAL, affineCorrection.s, update );
        residual.dualConic += update;
    }
}

template<typename Real>
void UpdateCombinedResidualUsingAffine
( const Real& sigma,
  const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const AffineLPState<Real>& state,
  const AffineLPSolution<Matrix<Real>>& solution,
  const AffineLPSolution<Matrix<Real>>& affineCorrection,
        AffineLPResidual<Matrix<Real>>& residual,
  bool compositeNewton )
{
    EL_DEBUG_CSE
    Shift( residual.dualConic, -sigma*state.barrier );
    // Use the interpretation of Section 2 of [CITATION]
    //
    //   R.A. Tapia, Y. Zhang, M. Saltzman, and A. Weiser's
    //   "The Mehrotra Predictor-Corrector Interior-Point Method as a
    //   Perturbed Composite Newton Method", July, 1990
    //   http://www.caam.rice.edu/tech_reports/1990/TR90-17.pdf
    //
    // as a means of generalizing the classic Mehrotra predictor-corrector
    // to approximate feasibility by updating with
    // F(solution+affineCorrection) rather than (0,0,0,dsAff o dzAff).
    //
    // With that said, it appears that neither the original or generalized
    // algorithm is generally helpful on the netlib/LP_data examples
    // (at least, on top of the current regularization and dynamic scaling).
    //
    const bool infeasibleCompositeNewton = true;
    if( compositeNewton && infeasibleCompositeNewton )
    {
        Matrix<Real> input, update;

        // residual.primalEquality +=
        //   A (x + dxAff) - b - (yRegSmall + yRegLarge) dyAff
        Zeros( update, problem.A.Height(), 1 );
        input = solution.x;
        input += affineCorrection.x;
        Multiply( NORMAL, Real(1), problem.A, input, Real(1), update );
        Axpy( Real(-1), problem.b, update );
        Axpy( -state.yRegSmall-state.yRegLarge, affineCorrection.y, update );
        const Real primalEqualityUpdateNrm2 = FrobeniusNorm( update );
        residual.primalEquality += update;

        // residual.dualEquality +=
        //   c + A^T (y + dyAff) + G^T (z + dzAff) +
        //   (xRegSmall + xRegLarge) dxAff
        Zeros( update, problem.A.Width(), 1 );
        input += problem.c;
        input = solution.y;
        input += affineCorrection.y;
        Multiply( TRANSPOSE, Real(1), problem.A, input, Real(1), update );
        input = solution.z;
        input += affineCorrection.z;
        Multiply( TRANSPOSE, Real(1), problem.G, input, Real(1), update );
        Axpy( state.xRegSmall+state.xRegLarge, affineCorrection.x, update );
        const Real dualEqualityUpdateNrm2 = FrobeniusNorm( update );
        residual.dualEquality += update;

        // residual.primalConic +=
        //   G (x + dxAff) + (s + dsAff) - h - (zRegSmall + zRegLarge) dzAff
        Zeros( update, problem.G.Height(), 1 );
        input = solution.x;
        input += affineCorrection.x;
        Multiply( NORMAL, Real(1), problem.G, input, Real(1), update );
        input = solution.s;
        input += affineCorrection.s;
        update += input;
        Axpy( Real(-1), problem.h, update );
        Axpy( -state.zRegSmall-state.zRegLarge, affineCorrection.z, update );
        const Real primalConicUpdateNrm2 = FrobeniusNorm( update );
        residual.primalConic += update;

        // residual.dualConic += (s + dsAff) o (z + dzAff)
        // TODO(poulson): Incorporate zPivot instead of z?
        input = solution.s;
        input += affineCorrection.s;
        update = solution.z;
        update += affineCorrection.z;
        DiagonalScale( LEFT, NORMAL, input, update );
        const Real dualConicUpdateNrm2 = FrobeniusNorm( update );
        residual.dualConic += update;
    }
    else if( compositeNewton )
    {
        // Use the approximation of the composite Newton method that assumes
        // the update would provide strict feasibility.

        // r_mu := dsAff o dzAff
        // ---------------------
        Matrix<Real> update( affineCorrection.z );
        DiagonalScale( LEFT, NORMAL, affineCorrection.s, update );
        residual.dualConic += update;
    }
}

template<typename Real>
void UpdateCombinedDualConicResidualWithoutAffine
( const Real& sigma,
  const Real& mu,
  const Real& complementRatio,
  bool largeCompRatio,
  const AffineLPSolution<Matrix<Real>>& solution,
        AffineLPResidual<Matrix<Real>>& residual,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int k = solution.s.Height();
    if( largeCompRatio && ctrl.softDualityTargets )
    {
        // Attempt to correct s o z entrywise into
        // [lowerTargetRatio,upperTargetRatio]*sigma*mu.
        Real lowerTargetRatio =
          Pow(complementRatio,ctrl.lowerTargetRatioLogCompRatio);
        Real upperTargetRatio =
          Pow(complementRatio,ctrl.upperTargetRatioLogCompRatio);
        lowerTargetRatio =
          Max( lowerTargetRatio,
            Pow(ctrl.maxComplementRatio,ctrl.lowerTargetRatioLogMaxCompRatio) );
        upperTargetRatio =
          Min( upperTargetRatio,
            Pow(ctrl.maxComplementRatio,ctrl.upperTargetRatioLogMaxCompRatio) );
        if( ctrl.print )
            Output
            ("  lowerTargetRatio=",lowerTargetRatio,
             ", upperTargetRatio=",upperTargetRatio);
        const Real lowerTarget = lowerTargetRatio*sigma*mu;
        const Real upperTarget = upperTargetRatio*sigma*mu;
        Int numBelowLower=0, numInterior=0, numAboveUpper=0;
        for( Int i=0; i<k; ++i )
        {
            const Real prod = solution.s(i)*solution.z(i);
            if( prod < lowerTarget )
            {
                residual.dualConic(i) -= lowerTarget;
                ++numBelowLower;
            }
            else if( prod > upperTarget )
            {
                residual.dualConic(i) -= upperTarget;
                ++numAboveUpper;
            }
            else
            {
                residual.dualConic(i) -= prod;
                ++numInterior;
            }
        }
        if( ctrl.print )
            Output
            ("  num targets (belowLower,interior,aboveUpper)=(",
             numBelowLower,",",numInterior,",",numAboveUpper,")");
    }
    else
    {
        Shift( residual.dualConic, -sigma*mu );
    }
}

// Despite having a factorization of 'J', we will solve the linear system
//
//   (D J D) (inv(D) x) = D b,
//
// where 'D = diag(diagonalScale)', rather than the original system, 'J x = b'.
//
// Without the usage of an iterative method (such as iterative refinement or
// FGMRES(k)), this is does not change the solution process. Its primary effect
// is to reweight the norm used to detect convergence in a way that helps to
// promote relative accuracy in each of the subvectors of the solution.
//
template<typename Real>
bool AttemptToSolve
( const AffineLPState<Real>& state,
  const Matrix<Real>& JFactored,
  const Matrix<Real>& dSub,
  const Permutation& permutation,
  const Matrix<Real>& diagonalScale,
        Matrix<Real>& rhs,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_ONLY(auto callStack = CopyCallStack())
    try
    {
        ldl::SolveAfter( JFactored, dSub, permutation, rhs, false );
    }
    catch( const std::exception& except )
    {
        if( ctrl.print )
            Output
            ("WARNING: solve failed with error: ",except.what());
        EL_DEBUG_ONLY(SetCallStack(callStack))
        return false;
    }
    return true;
}

template<typename Real>
bool AttemptToSolve
( const AffineLPState<Real>& state,
  const SparseLDLFactorization<Real>& sparseLDLFact,
  const SparseMatrix<Real>& JOrigScaled,
  const Matrix<Real>& regLargeScaled,
  const Matrix<Real>& diagonalScale,
        Matrix<Real>& rhs,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_ONLY(auto callStack = CopyCallStack())
    RegSolveInfo<Real> solveInfo;
    try
    {
        solveInfo =
          reg_ldl::RegularizedSolveAfter
          ( JOrigScaled, regLargeScaled, diagonalScale,
            sparseLDLFact, rhs,
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

template<typename Real>
bool AttemptToSolveForUpdateUsingScaling
( const AffineLPSolution<Matrix<Real>>& solution,
  const AffineLPState<Real>& state,
        SparseMatrix<Real>& JOrigScaled,
        Matrix<Real>& regLargeScaled,
  const SparseLDLFactorization<Real>& sparseLDLFact,
  const Matrix<Real>& zDoublePivot,
  const AffineLPResidual<Matrix<Real>>& residual,
        Real& dxCombinedNrm2Last,
        Real& dyCombinedNrm2Last,
        Real& dzCombinedNrm2Last,
        Matrix<Real>& packedVector,
        AffineLPSolution<Matrix<Real>>& correction,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = solution.y.Height();
    const Int n = solution.x.Height();
    const Int k = solution.s.Height();
    const Real minNrm2Last =
      Min( Min(dxCombinedNrm2Last,dyCombinedNrm2Last), dzCombinedNrm2Last );
    const Real maxNrm2Last =
      Max( Max(dxCombinedNrm2Last,dyCombinedNrm2Last), dzCombinedNrm2Last );

    // Construct the new KKT RHS
    // -------------------------
    KKTRHS
    ( residual.dualEquality,
      residual.primalEquality,
      residual.primalConic,
      residual.dualConic,
      zDoublePivot, packedVector );

    // Solve for the proposed (scaled) step
    // ------------------------------------
    Matrix<Real> diagonalScale;
    Zeros( diagonalScale, n+m+k, 1 );
    const Real maxRescaleRatio =
      Pow(limits::Epsilon<Real>(),ctrl.maxRescaleRatioLogEps);
    const Real baselineNrm2 =
      maxNrm2Last > maxRescaleRatio*minNrm2Last ?
      maxNrm2Last / maxRescaleRatio :
      minNrm2Last;
    if( dxCombinedNrm2Last < baselineNrm2 )
    {
        for( Int i=0; i<n; ++i )
            diagonalScale(i) = Real(1) / maxRescaleRatio;
    }
    else
    {
        for( Int i=0; i<n; ++i )
            diagonalScale(i) = Real(1) /
              Min(maxNrm2Last/dxCombinedNrm2Last,maxRescaleRatio);
    }
    if( dyCombinedNrm2Last < baselineNrm2 )
    {
        for( Int i=0; i<m; ++i )
            diagonalScale(i+n) = Real(1) / maxRescaleRatio;
    }
    else
    {
        for( Int i=0; i<m; ++i )
            diagonalScale(i+n) = Real(1) /
              Min(maxNrm2Last/dyCombinedNrm2Last,maxRescaleRatio);
    }
    if( dzCombinedNrm2Last < baselineNrm2 )
    {
        for( Int i=0; i<k; ++i )
            diagonalScale(i+n+m) = Real(1) / maxRescaleRatio;
    }
    else
    {
        for( Int i=0; i<k; ++i )
            diagonalScale(i+n+m) = Real(1) /
              Min(maxNrm2Last/dzCombinedNrm2Last,maxRescaleRatio);
    }
    DiagonalScale( LEFT, NORMAL, diagonalScale, JOrigScaled );
    DiagonalScale( RIGHT, NORMAL, diagonalScale, JOrigScaled );
    DiagonalScale( LEFT, NORMAL, diagonalScale, regLargeScaled );
    DiagonalScale( LEFT, NORMAL, diagonalScale, regLargeScaled );
    DiagonalScale( LEFT, NORMAL, diagonalScale, packedVector );
    const bool solved =
      AttemptToSolve
      ( state, sparseLDLFact, JOrigScaled, regLargeScaled, diagonalScale,
        packedVector, ctrl );
    if( !solved )
        return false;

    DiagonalScale( LEFT, NORMAL, diagonalScale, packedVector );
    ExpandSolution
    ( m, n, packedVector, residual.dualConic, solution.s, zDoublePivot,
      correction.x, correction.y, correction.z, correction.s );

    dxCombinedNrm2Last = FrobeniusNorm( correction.x );
    dyCombinedNrm2Last = FrobeniusNorm( correction.y );
    dzCombinedNrm2Last = FrobeniusNorm( correction.z );

    return true;
}

template<typename Real>
bool AttemptToSolveForUpdate
( const AffineLPSolution<Matrix<Real>>& solution,
  const AffineLPState<Real>& state,
  const Matrix<Real>& JFactored,
  const Matrix<Real>& dSub,
  const Permutation& permutation,
  const Matrix<Real>& zDoublePivot,
  const AffineLPResidual<Matrix<Real>>& residual,
        Matrix<Real>& packedVector,
        AffineLPSolution<Matrix<Real>>& correction,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = solution.y.Height();
    const Int n = solution.x.Height();
    const Int k = solution.s.Height();

    // Construct the new KKT RHS
    // -------------------------
    KKTRHS
    ( residual.dualEquality,
      residual.primalEquality,
      residual.primalConic,
      residual.dualConic,
      zDoublePivot, packedVector );

    // Solve for the proposed (scaled) step
    // ------------------------------------
    Matrix<Real> diagonalScale;
    Ones( diagonalScale, n+m+k, 1 );
    const bool solved =
      AttemptToSolve
      ( state, JFactored, dSub, permutation, diagonalScale, packedVector,
        ctrl );
    if( !solved )
        return false;

    ExpandSolution
    ( m, n, packedVector, residual.dualConic, solution.s, zDoublePivot,
      correction.x, correction.y, correction.z, correction.s );
    return true;
}

template<typename Real>
bool AttemptToSolveForUpdate
( const AffineLPSolution<Matrix<Real>>& solution,
  const AffineLPState<Real>& state,
        SparseMatrix<Real>& JOrigScaled,
        Matrix<Real>& regLargeScaled,
  const SparseLDLFactorization<Real>& sparseLDLFact,
  const Matrix<Real>& zDoublePivot,
  const AffineLPResidual<Matrix<Real>>& residual,
        Matrix<Real>& packedVector,
        AffineLPSolution<Matrix<Real>>& correction,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = solution.y.Height();
    const Int n = solution.x.Height();
    const Int k = solution.s.Height();

    // Construct the new KKT RHS
    // -------------------------
    KKTRHS
    ( residual.dualEquality,
      residual.primalEquality,
      residual.primalConic,
      residual.dualConic,
      zDoublePivot, packedVector );

    // Solve for the proposed (scaled) step
    // ------------------------------------
    Matrix<Real> diagonalScale;
    Ones( diagonalScale, n+m+k, 1 );
    const bool solved =
      AttemptToSolve
      ( state, sparseLDLFact, JOrigScaled, regLargeScaled, diagonalScale,
        packedVector, ctrl );
    if( !solved )
        return false;

    ExpandSolution
    ( m, n, packedVector, residual.dualConic, solution.s, zDoublePivot,
      correction.x, correction.y, correction.z, correction.s );
    return true;
}

template<typename Real>
void CheckConvergence
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& origProblem,
  const DenseAffineLPEquilibration<Real>& equilibration,
  const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const AffineLPSolution<Matrix<Real>>& solution,
        Matrix<Real>& zPivot,
        Matrix<Real>& zPerturb,
        AffineLPResidual<Matrix<Real>>& residual,
        AffineLPState<Real>& state,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real epsilon = limits::Epsilon<Real>();
    const Real zMinPivotValue = Pow(epsilon,ctrl.zMinPivotValueLogEps);

    // Form the unequlibrated approximate solution.
    AffineLPSolution<Matrix<Real>> origSolution;
    UndoEquilibration( solution, equilibration, origSolution );

    // A x - b
    // -------
    residual.primalEquality = problem.b;
    Gemv
    ( NORMAL, Real(1), problem.A, solution.x,
      Real(-1), residual.primalEquality );

    // || A x - b ||_2
    // ---------------
    state.primalEqualityInfeasNrm2 = FrobeniusNorm( residual.primalEquality );

    // || A x - b ||_2 / (1 + || b ||_2)
    // ---------------------------------
    state.primalEqualityInfeasNrm2Rel =
      state.primalEqualityInfeasNrm2 / (1 + state.bNrm2);

    // c + A^T y + G^T z
    // -----------------
    residual.dualEquality = problem.c;
    Gemv
    ( TRANSPOSE, Real(1), problem.A, solution.y,
      Real(1), residual.dualEquality );
    Gemv
    ( TRANSPOSE, Real(1), problem.G, solution.z,
      Real(1), residual.dualEquality );

    // || c + A^T y + G^T z ||_2
    // -------------------------
    state.dualEqualityInfeasNrm2 = FrobeniusNorm( residual.dualEquality );

    // || c + A^T y + G^T z ||_2 / (1 + || c ||_2)
    // -------------------------------------------
    state.dualEqualityInfeasNrm2Rel =
      state.dualEqualityInfeasNrm2 / (1 + state.cNrm2);

    // G x + s - h
    // -----------
    residual.primalConic = problem.h;
    Gemv
    ( NORMAL, Real(1), problem.G, solution.x,
      Real(-1), residual.primalConic );
    residual.primalConic += solution.s;

    // || G x + s - h ||_2
    // -------------------
    state.primalConicInfeasNrm2 = FrobeniusNorm( residual.primalConic );
    state.primalConicInfeasNrm2Rel =
      state.primalConicInfeasNrm2 / (1 + state.hNrm2);

    state.equalityError =
      Max(state.primalEqualityInfeasNrm2Rel,
          state.dualEqualityInfeasNrm2Rel);
    state.infeasError = Max(state.primalConicInfeasNrm2Rel,state.equalityError);

    // Compute the regularized 'z' pivot.
    zPivot = solution.z;
    const Int k = solution.z.Height();
    Zeros( zPerturb, k, 1 );
    for( Int i=0; i<k; ++i )
    {
        if( zPivot(i) < zMinPivotValue )
        {
            zPerturb(i) = zMinPivotValue - zPivot(i);
            zPivot(i) = zMinPivotValue;
        }
    }
    const Real zPerturbNrm2 = FrobeniusNorm( zPerturb );
    // s o zPivot
    // ----------
    residual.dualConic = zPivot;
    DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

    // Form objective convergence information relative to the original scaling.
    state.gapOrig = Dot( origSolution.s, origSolution.z );
    state.primalObjOrig =
      PrimalObjective<Real>( origProblem, origSolution );
    state.dualObjOrig =
      DualObjective<Real>( origProblem, origSolution );
    state.relCompGapOrig =
      RelativeComplementarityGap
      ( state.primalObjOrig, state.dualObjOrig, state.gapOrig );
    state.relObjGapOrig =
      RelativeObjectiveGap( state.primalObjOrig, state.dualObjOrig );
    state.maxRelGapOrig = Max( state.relCompGapOrig, state.relObjGapOrig );

    // Form objective convergence information relative to the current scaling.
    state.gap = Dot( solution.s, solution.z );
    state.primalObj = PrimalObjective<Real>( problem, solution );
    state.dualObj = DualObjective<Real>( problem, solution );
    state.relCompGap =
      RelativeComplementarityGap( state.primalObj, state.dualObj, state.gap );
    state.relObjGap = RelativeObjectiveGap( state.primalObj, state.dualObj );
    state.maxRelGap = Max( state.relCompGap, state.relObjGap );

    state.dimacsError = Max(state.maxRelGap,state.infeasError);
    state.dimacsErrorOrig = Max(state.maxRelGapOrig,state.infeasError);

    state.metTolerances =
      state.infeasError <= Pow(epsilon,ctrl.infeasibilityTolLogEps) &&
      state.relCompGap <=
        Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
      state.relObjGap <= Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);

    state.metTolerancesOrig =
      state.metTolerances &&
      state.relCompGapOrig <=
        Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
      state.relObjGapOrig <=
        Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);

    state.complementRatio = pos_orth::ComplementRatio( solution.s, solution.z );

    if( ctrl.print )
    {
        Output
        ("iter ",state.numIts,":\n",Indent(),
         "  || primalEqualityInfeas ||_2 / (1 + || b ||_2) = ",
         state.primalEqualityInfeasNrm2Rel,"\n",Indent(),
         "  || dualEqualityInfeas   ||_2 / (1 + || c ||_2) = ",
         state.dualEqualityInfeasNrm2Rel,"\n",Indent(),
         "  || primalConicInfeas    ||_2 / (1 + || h ||_2) = ",
         state.primalConicInfeasNrm2Rel,"\n",Indent(),
         "\n",Indent(),
         "  orig   s^T z        = ",state.gapOrig,"\n",Indent(),
         "  orig   primal       = ",state.primalObjOrig,"\n",Indent(),
         "  orig   dual         = ",state.dualObjOrig,"\n",Indent(),
         "  orig   rel obj gap  = ",state.relObjGapOrig,"\n",Indent(),
         "  orig   rel comp gap = ",state.relCompGapOrig,"\n",Indent(),
         "  orig   dimacs error = ",state.dimacsErrorOrig,"\n",Indent(),
         "\n",Indent(),
         "  scaled s^T z        = ",state.gap,"\n",Indent(),
         "  scaled primal       = ",state.primalObj,"\n",Indent(),
         "  scaled dual         = ",state.dualObj,"\n",Indent(),
         "  scaled rel obj gap  = ",state.relObjGap,"\n",Indent(),
         "  scaled rel comp gap = ",state.relCompGap,"\n",Indent(),
         "  scaled dimacs error = ",state.dimacsError,"\n",Indent(),
         "\n",Indent(),
         "  complement ratio    = ",state.complementRatio);
        if( zPerturbNrm2 > Real(0) )
            Output
            ("\n",Indent(),
             "  || zPerturb ||_2 = ",zPerturbNrm2);
    }
}

template<typename Real>
void CheckConvergence
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& origProblem,
  const SparseAffineLPEquilibration<Real>& equilibration,
  const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const AffineLPSolution<Matrix<Real>>& solution,
        Matrix<Real>& zPivot,
        Matrix<Real>& zPerturb,
        AffineLPResidual<Matrix<Real>>& residual,
        AffineLPState<Real>& state,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Real epsilon = limits::Epsilon<Real>();
    const Real zMinPivotValue = Pow(epsilon,ctrl.zMinPivotValueLogEps);

    // Form the unequlibrated approximate solution.
    AffineLPSolution<Matrix<Real>> origSolution;
    UndoEquilibration( solution, equilibration, origSolution );

    // A x - b
    // -------
    residual.primalEquality = problem.b;
    Multiply
    ( NORMAL, Real(1), problem.A, solution.x,
      Real(-1), residual.primalEquality );

    // || A x - b ||_2
    // ---------------
    state.primalEqualityInfeasNrm2 = FrobeniusNorm( residual.primalEquality );

    // || A x - b ||_2 / (1 + || b ||_2)
    // ---------------------------------
    state.primalEqualityInfeasNrm2Rel =
      state.primalEqualityInfeasNrm2 / (1 + state.bNrm2);

    // c + A^T y + G^T z
    // -----------------
    residual.dualEquality = problem.c;
    Multiply
    ( TRANSPOSE, Real(1), problem.A, solution.y,
      Real(1), residual.dualEquality );
    Multiply
    ( TRANSPOSE, Real(1), problem.G, solution.z,
      Real(1), residual.dualEquality );

    // || c + A^T y + G^T z ||_2
    // -------------------------
    state.dualEqualityInfeasNrm2 = FrobeniusNorm( residual.dualEquality );

    // || c + A^T y + G^T z ||_2 / (1 + || c ||_2)
    // -------------------------------------------
    state.dualEqualityInfeasNrm2Rel =
      state.dualEqualityInfeasNrm2 / (1 + state.cNrm2);

    // G x + s - h
    // -----------
    residual.primalConic = problem.h;
    Multiply
    ( NORMAL, Real(1), problem.G, solution.x,
      Real(-1), residual.primalConic );
    residual.primalConic += solution.s;

    // || G x + s - h ||_2
    // -------------------
    state.primalConicInfeasNrm2 = FrobeniusNorm( residual.primalConic );
    state.primalConicInfeasNrm2Rel =
      state.primalConicInfeasNrm2 / (1 + state.hNrm2);

    state.equalityError =
      Max(state.primalEqualityInfeasNrm2Rel,
          state.dualEqualityInfeasNrm2Rel);
    state.infeasError = Max(state.primalConicInfeasNrm2Rel,state.equalityError);

    // Compute the regularized 'z' pivot.
    zPivot = solution.z;
    const Int k = solution.z.Height();
    Zeros( zPerturb, k, 1 );
    for( Int i=0; i<k; ++i )
    {
        if( zPivot(i) < zMinPivotValue )
        {
            zPerturb(i) = zMinPivotValue - zPivot(i);
            zPivot(i) = zMinPivotValue;
        }
    }
    const Real zPerturbNrm2 = FrobeniusNorm( zPerturb );
    // s o zPivot
    // ----------
    residual.dualConic = zPivot;
    DiagonalScale( LEFT, NORMAL, solution.s, residual.dualConic );

    // Form objective convergence information relative to the original scaling.
    state.gapOrig = Dot( origSolution.s, origSolution.z );
    state.primalObjOrig =
      PrimalObjective<Real>( origProblem, origSolution );
    state.dualObjOrig =
      DualObjective<Real>( origProblem, origSolution );
    state.relCompGapOrig =
      RelativeComplementarityGap
      ( state.primalObjOrig, state.dualObjOrig, state.gapOrig );
    state.relObjGapOrig =
      RelativeObjectiveGap( state.primalObjOrig, state.dualObjOrig );
    state.maxRelGapOrig = Max( state.relCompGapOrig, state.relObjGapOrig );

    // Form objective convergence information relative to the current scaling.
    state.gap = Dot( solution.s, solution.z );
    state.primalObj = PrimalObjective<Real>( problem, solution );
    state.dualObj = DualObjective<Real>( problem, solution );
    state.relCompGap =
      RelativeComplementarityGap( state.primalObj, state.dualObj, state.gap );
    state.relObjGap = RelativeObjectiveGap( state.primalObj, state.dualObj );
    state.maxRelGap = Max( state.relCompGap, state.relObjGap );

    state.dimacsError = Max(state.maxRelGap,state.infeasError);
    state.dimacsErrorOrig = Max(state.maxRelGapOrig,state.infeasError);

    state.metTolerances =
      state.infeasError <= Pow(epsilon,ctrl.infeasibilityTolLogEps) &&
      state.relCompGap <=
        Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
      state.relObjGap <= Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);

    state.metTolerancesOrig =
      state.metTolerances &&
      state.relCompGapOrig <=
        Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
      state.relObjGapOrig <=
        Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);

    state.complementRatio = pos_orth::ComplementRatio( solution.s, solution.z );

    if( ctrl.print )
    {
        Output
        ("iter ",state.numIts,":\n",Indent(),
         "  || primalEqualityInfeas ||_2 / (1 + || b ||_2) = ",
         state.primalEqualityInfeasNrm2Rel,"\n",Indent(),
         "  || dualEqualityInfeas   ||_2 / (1 + || c ||_2) = ",
         state.dualEqualityInfeasNrm2Rel,"\n",Indent(),
         "  || primalConicInfeas    ||_2 / (1 + || h ||_2) = ",
         state.primalConicInfeasNrm2Rel,"\n",Indent(),
         "\n",Indent(),
         "  orig   s^T z        = ",state.gapOrig,"\n",Indent(),
         "  orig   primal       = ",state.primalObjOrig,"\n",Indent(),
         "  orig   dual         = ",state.dualObjOrig,"\n",Indent(),
         "  orig   rel obj gap  = ",state.relObjGapOrig,"\n",Indent(),
         "  orig   rel comp gap = ",state.relCompGapOrig,"\n",Indent(),
         "  orig   dimacs error = ",state.dimacsErrorOrig,"\n",Indent(),
         "\n",Indent(),
         "  scaled s^T z        = ",state.gap,"\n",Indent(),
         "  scaled primal       = ",state.primalObj,"\n",Indent(),
         "  scaled dual         = ",state.dualObj,"\n",Indent(),
         "  scaled rel obj gap  = ",state.relObjGap,"\n",Indent(),
         "  scaled rel comp gap = ",state.relCompGap,"\n",Indent(),
         "  scaled dimacs error = ",state.dimacsError,"\n",Indent(),
         "\n",Indent(),
         "  complement ratio    = ",state.complementRatio);
        if( zPerturbNrm2 > Real(0) )
            Output
            ("\n",Indent(),
             "  || zPerturb ||_2 = ",zPerturbNrm2);
    }
}

template<typename Real>
void PrintErrorResiduals
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
  const AffineLPSolution<Matrix<Real>>& affineCorrection,
  const AffineLPResidual<Matrix<Real>>& residual,
  const AffineLPState<Real>& state )
{
    EL_DEBUG_CSE
    AffineLPResidual<Matrix<Real>> error;

    error.primalEquality = residual.primalEquality;
    Gemv
    ( NORMAL, Real(1), problem.A, affineCorrection.x,
      Real(1), error.primalEquality );
    Axpy( -state.yRegSmall, affineCorrection.y, error.primalEquality );
    Axpy( -state.yRegLarge, affineCorrection.y, error.primalEquality );
    const Real dxErrorNrm2 = Nrm2( error.primalEquality );

    error.dualEquality = residual.dualEquality;
    Axpy( state.xRegSmall, affineCorrection.x, error.dualEquality );
    Axpy( state.xRegLarge, affineCorrection.x, error.dualEquality );
    Gemv
    ( TRANSPOSE, Real(1), problem.A, affineCorrection.y,
      Real(1), error.dualEquality );
    Gemv
    ( TRANSPOSE, Real(1), problem.G, affineCorrection.z,
      Real(1), error.dualEquality );
    const Real dyErrorNrm2 = Nrm2( error.dualEquality );

    error.primalConic = residual.primalConic;
    error.primalConic += affineCorrection.s;
    Gemv
    ( NORMAL, Real(1), problem.G, affineCorrection.x,
      Real(1), error.primalConic );
    Axpy( -state.zRegSmall, affineCorrection.z, error.primalConic );
    Axpy( -state.zRegLarge, affineCorrection.z, error.primalConic );
    const Real dzErrorNrm2 = Nrm2( error.primalConic );

    // TODO(poulson): error.dualConic

    Output
    ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
     dxErrorNrm2/(1+state.primalEqualityInfeasNrm2),"\n",Indent(),
     "|| dyError ||_2 / (1 + || r_c ||_2) = ",
     dyErrorNrm2/(1+state.dualEqualityInfeasNrm2),"\n",Indent(),
     "|| dzError ||_2 / (1 + || r_h ||_2) = ",
     dzErrorNrm2/(1+state.primalConicInfeasNrm2));
}

template<typename Real>
void PrintErrorResiduals
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const AffineLPSolution<Matrix<Real>>& affineCorrection,
  const AffineLPResidual<Matrix<Real>>& residual,
  const AffineLPState<Real>& state )
{
    EL_DEBUG_CSE
    AffineLPResidual<Matrix<Real>> error;

    error.primalEquality = residual.primalEquality;
    Multiply
    ( NORMAL, Real(1), problem.A, affineCorrection.x,
      Real(1), error.primalEquality );
    Axpy( -state.yRegSmall, affineCorrection.y, error.primalEquality );
    Axpy( -state.yRegLarge, affineCorrection.y, error.primalEquality );
    const Real dxErrorNrm2 = Nrm2( error.primalEquality );

    error.dualEquality = residual.dualEquality;
    Axpy( state.xRegSmall, affineCorrection.x, error.dualEquality );
    Axpy( state.xRegLarge, affineCorrection.x, error.dualEquality );
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
    Axpy( -state.zRegSmall, affineCorrection.z, error.primalConic );
    Axpy( -state.zRegLarge, affineCorrection.z, error.primalConic );
    const Real dzErrorNrm2 = Nrm2( error.primalConic );

    // TODO(poulson): error.dualConic

    Output
    ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
     dxErrorNrm2/(1+state.primalEqualityInfeasNrm2),"\n",Indent(),
     "|| dyError ||_2 / (1 + || r_c ||_2) = ",
     dyErrorNrm2/(1+state.dualEqualityInfeasNrm2),"\n",Indent(),
     "|| dzError ||_2 / (1 + || r_h ||_2) = ",
     dzErrorNrm2/(1+state.primalConicInfeasNrm2));
}

template<typename Real>
std::pair<Real,Real>
ComputePrimalDualStepSizes
( const AffineLPSolution<Matrix<Real>>& solution,
  const AffineLPSolution<Matrix<Real>>& correction,
  const Real& upperBound,
  const Real& maxStepRatio,
  bool forceSameStep )
{
    EL_DEBUG_CSE
    Real primalStep = pos_orth::MaxStep( solution.s, correction.s, upperBound );
    Real dualStep = pos_orth::MaxStep( solution.z, correction.z, upperBound );
    if( forceSameStep )
        primalStep = dualStep = Min(primalStep,dualStep);
    primalStep = Min(maxStepRatio*primalStep,Real(1));
    dualStep = Min(maxStepRatio*dualStep,Real(1));
    return std::make_pair( primalStep, dualStep );
}

template<typename Real>
Real MedianBarrier
( const Matrix<Real>& s, const Matrix<Real>& z )
{
    EL_DEBUG_CSE
    auto dualProd( s );
    DiagonalScale( LEFT, NORMAL, z, dualProd );
    auto medianPair = Median( dualProd );
    const Real barrierMedian = medianPair.value;
    return barrierMedian;
}

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
        Output("|| xScale ||_2 = ",MaxNorm(equilibration.xScale));
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
void RecenterAffineLPProblem
( const AffineLPSolution<Matrix<Real>>& center,
        AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem )
{
    EL_DEBUG_CSE
    // Recenter about (x,0;y,z), purposefully avoiding s. We put
    //   bHat := b - A x,
    //   hHat := h - G x,
    //   cHat := c + A^T y + G^T z.
    Gemv( NORMAL, Real(-1), problem.A, center.x, Real(1), problem.b );
    Gemv( NORMAL, Real(-1), problem.G, center.x, Real(1), problem.h );
    Gemv( TRANSPOSE, Real(1), problem.A, center.y, Real(1), problem.c );
    Gemv( TRANSPOSE, Real(1), problem.G, center.z, Real(1), problem.c );
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
void EquilibratedIPM
( const AffineLPProblem<Matrix<Real>,Matrix<Real>>& origProblem,
        DenseAffineLPEquilibration<Real>& equilibration,
        AffineLPProblem<Matrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int m = problem.A.Height();
    const Int n = problem.A.Width();
    const Int k = problem.G.Height();
    const Int degree = k;
    const Real epsilon = limits::Epsilon<Real>();
    const Real regIncreaseFactor = Pow(epsilon,ctrl.regIncreaseFactorLogEps);
    if( regIncreaseFactor <= Real(1) )
        LogicError("Regularization increase factor must be > 1");

    const Real minStepSize = Pow(limits::Epsilon<Real>(),Real(0.33));

    AffineLPResidual<Matrix<Real>> residual, error;
    AffineLPSolution<Matrix<Real>> affineCorrection, correction;
    AffineLPState<Real> state;
    state.bNrm2 = Nrm2( problem.b );
    state.cNrm2 = Nrm2( problem.c );
    state.hNrm2 = Nrm2( problem.h );
    const Real ANrm1 = OneNorm( problem.A );
    const Real GNrm1 = OneNorm( problem.G );
    const Real origOneNormEst = ANrm1 + GNrm1 + 1;
    if( ctrl.print )
    {

        Output("|| c ||_2 = ",state.cNrm2);
        Output("|| A ||_1 = ",ANrm1);
        Output("|| b ||_2 = ",state.bNrm2);
        Output("|| G ||_1 = ",GNrm1);
        Output("|| h ||_2 = ",state.hNrm2);
    }
    const Real xRegSmall0 = origOneNormEst*Pow(epsilon,ctrl.xRegSmallLogEps);
    const Real yRegSmall0 = origOneNormEst*Pow(epsilon,ctrl.yRegSmallLogEps);
    const Real zRegSmall0 = origOneNormEst*Pow(epsilon,ctrl.zRegSmallLogEps);
    const Real xRegLarge0 = origOneNormEst*Pow(epsilon,ctrl.xRegLargeLogEps);
    const Real yRegLarge0 = origOneNormEst*Pow(epsilon,ctrl.yRegLargeLogEps);
    const Real zRegLarge0 = origOneNormEst*Pow(epsilon,ctrl.zRegLargeLogEps);
    state.xRegLarge = xRegLarge0;
    state.yRegLarge = yRegLarge0;
    state.zRegLarge = zRegLarge0;
    state.xRegSmall = xRegSmall0;
    state.yRegSmall = yRegSmall0;
    state.zRegSmall = zRegSmall0;

    Matrix<Real> regSmall;
    regSmall.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )        regSmall(i) =  state.xRegSmall;
        else if( i < n+m ) regSmall(i) = -state.yRegSmall;
        else               regSmall(i) = -state.zRegSmall;
    }

    Matrix<Real> regLarge, regLargeScaled;
    regLarge.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )        regLarge(i) =  state.xRegLarge;
        else if( i < n+m ) regLarge(i) = -state.yRegLarge;
        else               regLarge(i) = -state.zRegLarge;
    }

    // TODO(poulson): Incorporate regularization in the initialization.
    Initialize
    ( problem, solution,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift );

    Matrix<Real> diagonalScale;
    Matrix<Real> J;
    Matrix<Real> dSub;
    Permutation permutation;
    Matrix<Real> packedVector;

    auto increaseRegularization = [&]() {
        if( ctrl.print )
            Output
            ("Increasing large regularization by a factor of ",
             regIncreaseFactor);
        regLarge *= regIncreaseFactor;
        state.xRegLarge *= regIncreaseFactor;
        state.yRegLarge *= regIncreaseFactor;
        state.zRegLarge *= regIncreaseFactor;
    };

    auto attemptToFactor = [&]()
      {
        EL_DEBUG_ONLY(auto callStack = CopyCallStack())
        try { LDL( J, dSub, permutation, false ); }
        catch(...)
        {
            EL_DEBUG_ONLY(SetCallStack(callStack))
            return false;
        }
        return true;
      };

    const Int indent = PushIndent();
    auto lastSolution( solution );
    bool lastStepForcedBarrierIncrease = false;
    Int numIts=0;
    for( ; state.numIts<=ctrl.maxIts;
         ++state.numIts,
         state.barrierOld=state.barrier,
         state.dimacsErrorOld=state.dimacsError,
         state.dimacsErrorOrigOld=state.dimacsErrorOrig )
    {
        if( state.backingOffEquilibration )
        {
            BackOffEquilibration
            ( equilibration, problem, state, solution, ctrl );
        }
        NeutralizePrimalAndDualNorms
        ( equilibration, problem, state, solution, ctrl );
        AssertPositiveOrthantMembership( solution );

        Matrix<Real> zPivot, zPerturb;
        CheckConvergence
        ( origProblem, equilibration, problem, solution,
          zPivot, zPerturb, residual, state, ctrl );

        if( state.metTolerancesOrig )
        {
            if( state.dimacsErrorOrig >=
                ctrl.minDimacsDecreaseRatio*state.dimacsErrorOrigOld )
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
        else if( !lastStepForcedBarrierIncrease &&
                 !state.backingOffEquilibration &&
                 state.dimacsErrorOrig > Real(1.05)*state.dimacsErrorOrigOld &&
                 state.dimacsError > Real(1.05)*state.dimacsErrorOld &&
                 state.dimacsErrorOld < Real(0.001) )
        {
            if( ctrl.print )
                Output
                ("Both the original and scaled DIMACS errors increased");
            solution = lastSolution;
            increaseRegularization();
            continue;
        }

        if( state.metTolerances && !state.metTolerancesOrig )
        {
            state.backingOffEquilibration = true;
            if( ctrl.print )
                Output("  will back off from equilibration in next step");
        }
        else
        {
            state.backingOffEquilibration = false;
        }

        // Factor the KKT system
        // =====================
        auto zDoublePivot( zPivot );
        zDoublePivot += zPerturb;
        KKT( problem.A, problem.G, solution.s, zDoublePivot, J );
        UpdateDiagonal( J, Real(1), regLarge ); 
        if( !attemptToFactor() )
        {
            increaseRegularization();
            continue;
        }

        // Decide whether or not we will attempt a predictor-corrector primarily
        // based upon the complementarity ratio.
        const bool largeCompRatio =
          state.complementRatio > ctrl.maxComplementRatio;
        // It has been observed that models such as wood1w can continually try
        // to drive down the barrier parameter despite relative infeasibility
        // (e.g., a relative complementarity gap of 10^-14 while the relative
        // dual infeasibility is on the order of 10^-7).
        const Real relativeComplementarityLowerBound =
          Min(Pow(state.infeasError,Real(1.)),Real(1e-3));
        const bool relCompGapIsTooSmall =
          state.relCompGap < relativeComplementarityLowerBound;
        const bool centeringStep = largeCompRatio || relCompGapIsTooSmall;
        const Real barrierClassical = Dot(solution.s,solution.z) / k;
        const Real barrierMedian = MedianBarrier( solution.s, solution.z );
        if( ctrl.print )
            Output
            ("  barrierMedian = ",barrierMedian,
             ", barrierClassical=",barrierClassical);
        if( centeringStep )
        {
            if( ctrl.print )
                Output("  will take a centering step");
            if( relCompGapIsTooSmall )
            {
                state.barrier = Max( state.barrierOld, barrierMedian );
                if( ctrl.print )
                    Output
                    ("  using barrier = Max( barrierOld=",state.barrierOld,
                     ", barrierMedian=",barrierMedian,")");
                lastStepForcedBarrierIncrease = true;
            }
            else
            {
                state.barrier = barrierMedian;
                if( ctrl.print )
                    Output
                    ("  using barrier = barrierMedian = ",barrierMedian);
                lastStepForcedBarrierIncrease = false;
            }
        }
        else
        {
            if( ctrl.print )
                Output("  will take a predictor-corrector step.");
            if( state.backingOffEquilibration )
            {
                state.barrier = barrierMedian;
                if( ctrl.print )
                    Output
                    ("  backing off equilibration, using barrierMedian = ",
                     barrierMedian);
            }
            else
            {
                state.barrier = Min( state.barrierOld, barrierMedian );
                if( ctrl.print )
                    Output
                    ("  using barrier = Min( barrierOld=",state.barrierOld,
                     ", barrierMedian=",barrierMedian," )");
            }
            lastStepForcedBarrierIncrease = false;
        }

        Real sigma;
        if( centeringStep )
        {
            sigma = 1;
            if( ctrl.print )
                Output("  freezing sigma at one");

            // Solve for the combined direction
            // ================================
            UpdateCombinedDualConicResidualWithoutAffine
            ( sigma, state.barrier, state.complementRatio, largeCompRatio,
              solution, residual, ctrl );
            const bool solvedForCombined =
              AttemptToSolveForUpdate
              ( solution, state, J, dSub, permutation,
                zDoublePivot, residual, packedVector, correction, ctrl );
            if( !solvedForCombined )
            {
                increaseRegularization();
                continue;
            }
        }
        else
        {
            // Compute the affine search direction
            // ===================================

            // Solve for the proposed step
            // ---------------------------
            KKTRHS
            ( residual.dualEquality,
              residual.primalEquality,
              residual.primalConic,
              residual.dualConic,
              zDoublePivot,
              packedVector );
            Ones( diagonalScale, n+m+k, 1 );
            const bool solvedForAffine =
              AttemptToSolve
              ( state, J, dSub, permutation,
                diagonalScale, packedVector, ctrl );
            if( !solvedForAffine )
            {
                increaseRegularization();
                continue;
            }
            ExpandSolution
            ( m, n, packedVector, residual.dualConic, solution.s, zDoublePivot,
              affineCorrection.x,
              affineCorrection.y,
              affineCorrection.z,
              affineCorrection.s );

            // Compute a centrality parameter
            // ==============================
            auto affineStep =
              ComputePrimalDualStepSizes
              ( solution, affineCorrection, Real(1), Real(1),
                ctrl.forceSameStep );
            if( ctrl.print )
            {
                Output("  affine primal step size: ",affineStep.first);
                Output("  affine dual   step size: ",affineStep.second);
            }
            // TODO(poulson): Experiment with several centrality choices?
            auto sAffine( solution.s );
            Axpy( affineStep.first, affineCorrection.s, sAffine );
            auto zAffine( solution.z );
            Axpy( affineStep.second, affineCorrection.z, zAffine );
            //const Real barrierAff = Dot(sAffine,zAffine) / degree;
            const Real barrierAff = MedianBarrier( sAffine, zAffine );
            if( ctrl.print )
                Output
                ("  barrierAff = ",barrierAff,", barrier = ",state.barrier);
            const Real sigma =
              ctrl.centralityRule
              (state.barrier,barrierAff,affineStep.first,affineStep.second);
            if( ctrl.print )
                Output("  sigma=",sigma);

            // Solve for the combined direction
            // ================================
            const bool compositeNewton =
              ctrl.compositeNewton &&
              state.infeasError < Pow(epsilon,ctrl.infeasibilityTolLogEps);
            UpdateCombinedResidualUsingAffine
            ( sigma, problem, state, solution, affineCorrection, residual,
              compositeNewton );
            const bool solvedForCombined =
              AttemptToSolveForUpdate
              ( solution, state, J, dSub, permutation,
                zDoublePivot, residual, packedVector, correction, ctrl );
            if( !solvedForCombined )
            {
                increaseRegularization();
                continue;
            }
        }

        // Update the current estimates
        // ============================
        lastSolution = solution;
        auto combinedStep = ComputePrimalDualStepSizes
        ( solution, correction,
          Real(1)/ctrl.maxStepRatio, ctrl.maxStepRatio, ctrl.forceSameStep );
        if( ctrl.print )
        {
            Output("  combined primal step size: ",combinedStep.first);
            Output("  combined dual   step size: ",combinedStep.second);
        }
        Axpy( combinedStep.first, correction.x, solution.x );
        Axpy( combinedStep.first, correction.s, solution.s );
        Axpy( combinedStep.second, correction.y, solution.y );
        Axpy( combinedStep.second, correction.z, solution.z );
        if( combinedStep.first < minStepSize &&
            combinedStep.second < minStepSize )
        {
            if( state.metTolerancesOrig )
            {
                break;
            }
            else
            {
                RuntimeError("Attempted step sizes less than ",minStepSize);
            }
        }
        if( ctrl.print )
            Output("");
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
        EquilibratedIPM
        ( problem, equilibration,
          equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );

        // Refinement seems to lead to very challenging problems (e.g., for
        // 80bau3b).
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
            if( ctrl.print )
            {
                Output("|| bHat ||_max = ",bHatMax);
                Output("|| hHat ||_max = ",hHatMax);
                Output("|| cHat ||_max = ",cHatMax);
                Output
                ("primalScale=",primalRefineScale,
                 ", dualScale=",dualRefineScale);
            }
            refinedProblem.b *= primalRefineScale;
            refinedProblem.h *= primalRefineScale;
            refinedProblem.c *= dualRefineScale;

            AffineLPSolution<Matrix<Real>> refinedSolution,
              equilibratedRefinedSolution;
            DenseAffineLPEquilibration<Real> trivialEquilibration;
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
            if( ctrl.print )
                Output("New primal: ",newPrimal,", new dual: ",newDual);
        }
    }
    else
    {
        DenseAffineLPEquilibration<Real> equilibration;
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
    const Real epsilon = limits::Epsilon<Real>();
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
        const Real gap = Dot( solution.s, solution.z );
        const Real mu = gap / k;

        // Check for convergence
        // =====================
        const Real primObj = PrimalObjective<Real>( problem, solution );
        const Real dualObj = DualObjective<Real>( problem, solution );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, gap );
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real maxRelGap = Max( relCompGap, relObjGap );

        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        Gemv
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(-1), residual.primalEquality );
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
        Gemv
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(-1), residual.primalConic );
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

        /*
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
        */

        // Compute a centrality parameter
        // ==============================
        Real primalAffineStep =
          pos_orth::MaxStep( solution.s, affineCorrection.s, Real(1) );
        Real dualAffineStep =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            primalAffineStep = dualAffineStep =
              Min(primalAffineStep,dualAffineStep);
        if( ctrl.print && commRank == 0 )
            Output
            ("primalAffineStep = ",primalAffineStep,
             ", dualAffineStep = ",dualAffineStep);
        // NOTE: dz and ds are used as temporaries
        correction.s = solution.s;
        correction.z = solution.z;
        Axpy( primalAffineStep,  affineCorrection.s, correction.s );
        Axpy( dualAffineStep, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.s,correction.z) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma =
          ctrl.centralityRule(mu,muAff,primalAffineStep,dualAffineStep);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.compositeNewton )
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
        Real primalStep =
          pos_orth::MaxStep( solution.s, correction.s, 1/ctrl.maxStepRatio );
        Real dualStep =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        if( ctrl.forceSameStep )
            primalStep = dualStep = Min(primalStep,dualStep);
        primalStep = Min(ctrl.maxStepRatio*primalStep,Real(1));
        dualStep = Min(ctrl.maxStepRatio*dualStep,Real(1));
        if( ctrl.print && commRank == 0 )
            Output("primalStep = ",primalStep,", dualStep = ",dualStep);
        Axpy( primalStep,  correction.x, solution.x );
        Axpy( primalStep,  correction.s, solution.s );
        Axpy( dualStep, correction.y, solution.y );
        Axpy( dualStep, correction.z, solution.z );
        if( primalStep == Real(0) && dualStep == Real(0) )
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
    const Real epsilon = limits::Epsilon<Real>();
    const Real regIncreaseFactor = Pow(epsilon,ctrl.regIncreaseFactorLogEps);
    if( regIncreaseFactor <= Real(1) )
        LogicError("Regularization increase factor must be > 1");
    const Real minStepSize = Pow(limits::Epsilon<Real>(),Real(0.33));

    AffineLPResidual<Matrix<Real>> residual, error;
    AffineLPSolution<Matrix<Real>> affineCorrection, correction;
    AffineLPState<Real> state;
    state.bNrm2 = FrobeniusNorm( problem.b );
    state.cNrm2 = FrobeniusNorm( problem.c );
    state.hNrm2 = FrobeniusNorm( problem.h );
    const Real twoNormEstA =
      TwoNormEstimate( problem.A, ctrl.twoNormKrylovBasisSize );
    const Real twoNormEstG =
      TwoNormEstimate( problem.G, ctrl.twoNormKrylovBasisSize );
    const Real origTwoNormEst = twoNormEstA + twoNormEstG + 1;
    if( ctrl.print )
    {
        Output("|| A ||_2 estimate: ",twoNormEstA);
        Output("|| G ||_2 estimate: ",twoNormEstG);
        Output("|| b ||_2 = ",state.bNrm2);
        Output("|| c ||_2 = ",state.cNrm2);
        Output("|| h ||_2 = ",state.hNrm2);
    }
    const Real xRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.xRegSmallLogEps);
    const Real yRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.yRegSmallLogEps);
    const Real zRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.zRegSmallLogEps);
    const Real xRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.xRegLargeLogEps);
    const Real yRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.yRegLargeLogEps);
    const Real zRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.zRegLargeLogEps);
    state.xRegLarge = xRegLarge0;
    state.yRegLarge = yRegLarge0;
    state.zRegLarge = zRegLarge0;
    state.xRegSmall = xRegSmall0;
    state.yRegSmall = yRegSmall0;
    state.zRegSmall = zRegSmall0;

    Matrix<Real> regSmall;
    regSmall.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )        regSmall(i) =  state.xRegSmall;
        else if( i < n+m ) regSmall(i) = -state.yRegSmall;
        else               regSmall(i) = -state.zRegSmall;
    }

    Matrix<Real> regLarge, regLargeScaled;
    regLarge.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )        regLarge(i) =  state.xRegLarge;
        else if( i < n+m ) regLarge(i) = -state.yRegLarge;
        else               regLarge(i) = -state.zRegLarge;
    }

    // Initialize the static portion of the KKT system
    // ===============================================
    SparseMatrix<Real> JStatic;
    StaticKKT
    ( problem.A, problem.G,
      Sqrt(state.xRegSmall), Sqrt(state.yRegSmall), Sqrt(state.zRegSmall),
      JStatic, false );
    JStatic.FreezeSparsity();

    SparseLDLFactorization<Real> sparseLDLFact;
    Initialize
    ( problem, solution, JStatic, regLarge,
      sparseLDLFact,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift, ctrl.solveCtrl );

    Matrix<Real> diagonalScale;
    SparseMatrix<Real> J, JOrig, JOrigScaled;
    Matrix<Real> packedVector;

    auto increaseRegularization = [&]() {
        if( ctrl.print )
            Output
            ("Increasing large regularization by a factor of ",
             regIncreaseFactor);
        regLarge *= regIncreaseFactor;
        state.xRegLarge *= regIncreaseFactor;
        state.yRegLarge *= regIncreaseFactor;
        state.zRegLarge *= regIncreaseFactor;
    };

    auto attemptToFactor = [&]() {
        EL_DEBUG_ONLY(auto callStack = CopyCallStack())
        try
        {
            if( state.numIts == 0 && ctrl.primalInit && ctrl.dualInit )
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

    const Int indent = PushIndent();
    auto lastSolution( solution );
    bool lastStepForcedBarrierIncrease = false;
    Real dxCombinedNrm2Last=Real(1),
         dyCombinedNrm2Last=Real(1),
         dzCombinedNrm2Last=Real(1);
    for( ; state.numIts<=ctrl.maxIts;
         ++state.numIts,
         state.barrierOld=state.barrier,
         state.dimacsErrorOld=state.dimacsError,
         state.dimacsErrorOrigOld=state.dimacsErrorOrig )
    {
        if( state.backingOffEquilibration )
        {
            BackOffEquilibration
            ( equilibration, problem, state, solution, JStatic, ctrl );
        }
        NeutralizePrimalAndDualNorms
        ( equilibration, problem, state, solution, ctrl );
        AssertPositiveOrthantMembership( solution );

        Matrix<Real> zPivot, zPerturb;
        CheckConvergence
        ( origProblem, equilibration, problem, solution,
          zPivot, zPerturb, residual, state, ctrl );

        if( state.metTolerancesOrig )
        {
            if( state.dimacsErrorOrig >=
                ctrl.minDimacsDecreaseRatio*state.dimacsErrorOrigOld )
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
        else if( !lastStepForcedBarrierIncrease &&
                 !state.backingOffEquilibration &&
                 state.dimacsErrorOrig > Real(1.05)*state.dimacsErrorOrigOld &&
                 state.dimacsError > Real(1.05)*state.dimacsErrorOld &&
                 state.dimacsErrorOld < Real(0.001) )
        {
            if( ctrl.print )
                Output
                ("Both the original and scaled DIMACS errors increased");
            solution = lastSolution;
            increaseRegularization();
            continue;
        }

        if( state.metTolerances && !state.metTolerancesOrig )
        {
            state.backingOffEquilibration = true;
            if( ctrl.print )
                Output("  will back off from equilibration in next step");
        }
        else
        {
            state.backingOffEquilibration = false;
        }

        // Factor the (rescaled) KKT system
        // ================================
        Ones( diagonalScale, n+m+k, 1 );
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        auto zDoublePivot( zPivot );
        zDoublePivot += zPerturb;
        FinishKKT( m, n, solution.s, zDoublePivot, JOrig );
        JOrigScaled = JOrig;
        JOrigScaled.FreezeSparsity();
        regLargeScaled = regLarge;
        J = JOrig;
        J.FreezeSparsity();
        UpdateDiagonal( J, Real(1), regLarge );
        if( !attemptToFactor() )
        {
            increaseRegularization();
            continue;
        }

        // Decide whether or not we will attempt a predictor-corrector primarily
        // based upon the complementarity ratio.
        const bool largeCompRatio =
          state.complementRatio > ctrl.maxComplementRatio;
        // It has been observed that models such as wood1w can continually try
        // to drive down the barrier parameter despite relative infeasibility
        // (e.g., a relative complementarity gap of 10^-14 while the relative
        // dual infeasibility is on the order of 10^-7).
        const Real relativeComplementarityLowerBound =
          Min(Pow(state.infeasError,Real(1.)),Real(1e-3));
        const bool relCompGapIsTooSmall =
          state.relCompGap < relativeComplementarityLowerBound;
        const bool centeringStep = largeCompRatio || relCompGapIsTooSmall;
        const Real barrierClassical = Dot(solution.s,solution.z) / k;
        const Real barrierMedian = MedianBarrier( solution.s, solution.z );
        if( ctrl.print )
            Output
            ("  barrierMedian = ",barrierMedian,
             ", barrierClassical=",barrierClassical);
        if( centeringStep )
        {
            if( ctrl.print )
                Output("  will take a centering step");
            if( relCompGapIsTooSmall )
            {
                state.barrier = Max( state.barrierOld, barrierMedian );
                if( ctrl.print )
                    Output
                    ("  using barrier = Max( barrierOld=",state.barrierOld,
                     ", barrierMedian=",barrierMedian,")");
                lastStepForcedBarrierIncrease = true;
            }
            else
            {
                state.barrier = barrierMedian;
                if( ctrl.print )
                    Output
                    ("  using barrier = barrierMedian = ",barrierMedian);
                lastStepForcedBarrierIncrease = false;
            }
        }
        else
        {
            if( ctrl.print )
                Output("  will take a predictor-corrector step.");
            if( state.backingOffEquilibration )
            {
                state.barrier = barrierMedian;
                if( ctrl.print )
                    Output
                    ("  backing off equilibration, using barrierMedian = ",
                     barrierMedian);
            }
            else
            {
                state.barrier = Min( state.barrierOld, barrierMedian );
                if( ctrl.print )
                    Output
                    ("  using barrier = Min( barrierOld=",state.barrierOld,
                     ", barrierMedian=",barrierMedian," )");
            }
            lastStepForcedBarrierIncrease = false;
        }

        Real sigma;
        if( centeringStep )
        {
            sigma = 1;
            if( ctrl.print )
                Output("  freezing sigma at one");

            // Solve for the combined direction
            // ================================
            UpdateCombinedDualConicResidualWithoutAffine
            ( sigma, state.barrier, state.complementRatio, largeCompRatio,
              solution, residual, ctrl );
            const bool solvedForCombined =
              AttemptToSolveForUpdate
              ( solution, state, JOrigScaled, regLargeScaled, sparseLDLFact,
                zDoublePivot, residual, packedVector, correction, ctrl );
            if( !solvedForCombined )
            {
                increaseRegularization();
                continue;
            }
        }
        else
        {
            // Compute the affine search direction
            // ===================================

            // Solve for the proposed step
            // ---------------------------
            KKTRHS
            ( residual.dualEquality,
              residual.primalEquality,
              residual.primalConic,
              residual.dualConic,
              zDoublePivot,
              packedVector );
            Ones( diagonalScale, n+m+k, 1 );
            const bool solvedForAffine =
              AttemptToSolve
              ( state, sparseLDLFact, JOrigScaled, regLargeScaled,
                diagonalScale, packedVector, ctrl );
            if( !solvedForAffine )
            {
                increaseRegularization();
                continue;
            }
            ExpandSolution
            ( m, n, packedVector, residual.dualConic, solution.s, zDoublePivot,
              affineCorrection.x,
              affineCorrection.y,
              affineCorrection.z,
              affineCorrection.s );

            // Compute a centrality parameter
            // ==============================
            auto affineStep =
              ComputePrimalDualStepSizes
              ( solution, affineCorrection, Real(1), Real(1),
                ctrl.forceSameStep );
            if( ctrl.print )
            {
                Output("  affine primal step size: ",affineStep.first);
                Output("  affine dual   step size: ",affineStep.second);
            }
            // TODO(poulson): Experiment with several centrality choices?
            auto sAffine( solution.s );
            Axpy( affineStep.first, affineCorrection.s, sAffine );
            auto zAffine( solution.z );
            Axpy( affineStep.second, affineCorrection.z, zAffine );
            //const Real barrierAff = Dot(sAffine,zAffine) / degree;
            const Real barrierAff = MedianBarrier( sAffine, zAffine );
            if( ctrl.print )
                Output
                ("  barrierAff = ",barrierAff,", barrier = ",state.barrier);
            const Real sigma =
              ctrl.centralityRule
              (state.barrier,barrierAff,affineStep.first,affineStep.second);
            if( ctrl.print )
                Output("  sigma=",sigma);

            // Solve for the combined direction
            // ================================
            const bool compositeNewton =
              ctrl.compositeNewton &&
              state.infeasError < Pow(epsilon,ctrl.infeasibilityTolLogEps);
            UpdateCombinedResidualUsingAffine
            ( sigma, problem, state, solution, affineCorrection, residual,
              compositeNewton );
            const bool solvedForCombined =
              AttemptToSolveForUpdateUsingScaling
              ( solution, state, JOrigScaled, regLargeScaled, sparseLDLFact,
                zDoublePivot, residual,
                dxCombinedNrm2Last, dyCombinedNrm2Last, dzCombinedNrm2Last,
                packedVector, correction, ctrl );
            if( !solvedForCombined )
            {
                increaseRegularization();
                continue;
            }
        }

        // Update the current estimates
        // ============================
        lastSolution = solution;
        auto combinedStep = ComputePrimalDualStepSizes
        ( solution, correction,
          Real(1)/ctrl.maxStepRatio, ctrl.maxStepRatio, ctrl.forceSameStep );
        if( ctrl.print )
        {
            Output("  combined primal step size: ",combinedStep.first);
            Output("  combined dual   step size: ",combinedStep.second);
        }
        Axpy( combinedStep.first, correction.x, solution.x );
        Axpy( combinedStep.first, correction.s, solution.s );
        Axpy( combinedStep.second, correction.y, solution.y );
        Axpy( combinedStep.second, correction.z, solution.z );
        if( combinedStep.first < minStepSize &&
            combinedStep.second < minStepSize )
        {
            if( state.metTolerancesOrig )
            {
                break;
            }
            else
            {
                RuntimeError("Attempted step sizes less than ",minStepSize);
            }
        }
        if( ctrl.print )
            Output("");
    }
    SetIndent( indent );
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

        // Refinement seems to lead to very challenging problems (e.g., for
        // 80bau3b).
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
            if( ctrl.print )
            {
                Output("|| bHat ||_max = ",bHatMax);
                Output("|| hHat ||_max = ",hHatMax);
                Output("|| cHat ||_max = ",cHatMax);
                Output
                ("primalScale=",primalRefineScale,
                 ", dualScale=",dualRefineScale);
            }
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
            if( ctrl.print )
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
    const Real epsilon = limits::Epsilon<Real>();
    const Real regIncreaseFactor = Pow(epsilon,ctrl.regIncreaseFactorLogEps);
    if( regIncreaseFactor <= Real(1) )
        LogicError("Regularization increase factor must be > 1");

    AffineLPState<Real> state;
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

    const Real xRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.xRegSmallLogEps);
    const Real yRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.yRegSmallLogEps);
    const Real zRegSmall0 = origTwoNormEst*Pow(epsilon,ctrl.zRegSmallLogEps);
    const Real xRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.xRegLargeLogEps);
    const Real yRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.yRegLargeLogEps);
    const Real zRegLarge0 = origTwoNormEst*Pow(epsilon,ctrl.zRegLargeLogEps);
    state.xRegLarge = xRegLarge0;
    state.yRegLarge = yRegLarge0;
    state.zRegLarge = zRegLarge0;
    state.xRegSmall = xRegSmall0;
    state.yRegSmall = yRegSmall0;
    state.zRegSmall = zRegSmall0;

    DistMultiVec<Real> regLarge(grid);
    regLarge.Resize( n+m+k, 1 );
    for( Int iLoc=0; iLoc<regLarge.LocalHeight(); ++iLoc )
    {
        const Int i = regLarge.GlobalRow(iLoc);
        if( i < n )
          regLarge.SetLocal( iLoc, 0, state.xRegLarge );
        else if( i < n+m )
          regLarge.SetLocal( iLoc, 0, -state.yRegLarge );
        else
          regLarge.SetLocal( iLoc, 0, -state.zRegLarge );
    }
    regLarge *= origTwoNormEst;

    // Construct the static part of the KKT system
    // ===========================================
    DistSparseMatrix<Real> JStatic(grid);
    StaticKKT
    ( problem.A, problem.G,
      Sqrt(state.xRegSmall), Sqrt(state.yRegSmall), Sqrt(state.zRegSmall),
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
        const Real gap = Dot( solution.s, solution.z );
        const Real primObj = PrimalObjective<Real>( problem, solution );
        const Real dualObj = DualObjective<Real>( problem, solution );
        const Real relCompGap =
          RelativeComplementarityGap( primObj, dualObj, gap );
        const Real relObjGap = RelativeObjectiveGap( primObj, dualObj );
        const Real maxRelGap = Max( relCompGap, relObjGap );

        // || A x - b ||_2 / (1 + || b ||_2)
        // ---------------------------------
        residual.primalEquality = problem.b;
        Multiply
        ( NORMAL, Real(1), problem.A, solution.x,
          Real(-1), residual.primalEquality );
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
        Multiply
        ( NORMAL, Real(1), problem.G, solution.x,
          Real(-1), residual.primalConic );
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

        /*
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
        */

        // Compute a centrality parameter
        // ==============================
        Real primalAffineStep =
          pos_orth::MaxStep( solution.s, affineCorrection.s, Real(1) );
        Real dualAffineStep =
          pos_orth::MaxStep( solution.z, affineCorrection.z, Real(1) );
        if( ctrl.forceSameStep )
            primalAffineStep = dualAffineStep =
              Min(primalAffineStep,dualAffineStep);
        if( ctrl.print && commRank == 0 )
            Output
            ("primalAffineStep = ",primalAffineStep,
             ", dualAffineStep = ",dualAffineStep);
        // NOTE: correction.z and correction.s are used as temporaries
        correction.s = solution.s;
        correction.z = solution.z;
        Axpy( primalAffineStep,  affineCorrection.s, correction.s );
        Axpy( dualAffineStep, affineCorrection.z, correction.z );
        const Real muAff = Dot(correction.s,correction.z) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma =
          ctrl.centralityRule(mu,muAff,primalAffineStep,dualAffineStep);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        Shift( residual.dualConic, -sigma*mu );
        if( ctrl.compositeNewton )
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
        Real primalStep =
          pos_orth::MaxStep( solution.s, correction.s, 1/ctrl.maxStepRatio );
        Real dualStep =
          pos_orth::MaxStep( solution.z, correction.z, 1/ctrl.maxStepRatio );
        primalStep = Min(ctrl.maxStepRatio*primalStep,Real(1));
        dualStep = Min(ctrl.maxStepRatio*dualStep,Real(1));
        if( ctrl.forceSameStep )
            primalStep = dualStep = Min(primalStep,dualStep);
        if( ctrl.print && commRank == 0 )
            Output("primalStep = ",primalStep,", dualStep = ",dualStep);
        Axpy( primalStep, correction.x, solution.x );
        Axpy( primalStep, correction.s, solution.s );
        Axpy( dualStep, correction.y, solution.y );
        Axpy( dualStep, correction.z, solution.z );
        if( primalStep == Real(0) && dualStep == Real(0) )
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
