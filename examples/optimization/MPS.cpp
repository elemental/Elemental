/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>

namespace lp_hsd {
namespace affine {

using namespace El;
using El::lp::affine::Ctrl;

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
struct AffineLPState
{
    Int numIterations=0;

    Real barrier=Real(0.1);
    Real barrierOld=Real(0.1);
    bool centeringStep=false;

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
    Real lowerComplementRatio=Real(1);
    Real upperComplementRatio=Real(1);
    bool largeComplementRatio=false;

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

    bool approxFeasible=false;
    bool metTolerances=false;
    bool metTolerancesOrig=false;

    bool backingOffEquilibration=false;
    bool lastStepForcedBarrierIncrease=false;
};

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
void StaticKKT
( const SparseMatrix<Real>& A,
  const SparseMatrix<Real>& G,
        Real gamma,
        Real delta,
        Real beta,
        SparseMatrix<Real>& J,
  bool onlyLower )
{
    EL_DEBUG_CSE
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    Zeros( J, n+m+k, n+m+k );
    const Int numEntriesA = A.NumEntries();
    const Int numEntriesG = G.NumEntries();
    if( onlyLower )
        J.Reserve( numEntriesA + numEntriesG + n+m+k );
    else
        J.Reserve( 2*numEntriesA + 2*numEntriesG + n+m+k );

    // Jxx = gamma^2*I
    // ===================
    for( Int i=0; i<n; ++i )
        J.QueueUpdate( i, i, gamma*gamma );

    // Jyy = -delta^2*I
    // ================
    for( Int i=0; i<m; ++i )
        J.QueueUpdate( i+n, i+n, -delta*delta );

    // Jyx = A (and Jxy = A^T)
    // =======================
    for( Int e=0; e<numEntriesA; ++e )
    {
        J.QueueUpdate( n+A.Row(e), A.Col(e), A.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( A.Col(e), n+A.Row(e), A.Value(e) );
    }

    // Jzx = G (and Jxz = G^T)
    // =======================
    for( Int e=0; e<numEntriesG; ++e )
    {
        J.QueueUpdate( n+m+G.Row(e), G.Col(e), G.Value(e) );
        if( !onlyLower )
            J.QueueUpdate( G.Col(e), n+m+G.Row(e), G.Value(e) );
    }

    // Jzz = -beta^2*I
    // ===============
    for( Int i=0; i<k; ++i )
        J.QueueUpdate( i+m+n, i+m+n, -beta*beta );

    J.ProcessQueues();
    J.FreezeSparsity();
}

template<typename Real>
void FinishKKT
( Int m, Int n,
  const Matrix<Real>& s,
  const Matrix<Real>& z,
        SparseMatrix<Real>& J )
{
    EL_DEBUG_CSE
    const Int k = s.Height();

    // Jzz = -z <> s
    // =============
    if( !J.FrozenSparsity() )
        J.Reserve( J.NumEntries()+k );
    for( Int e=0; e<k; ++e )
        J.QueueUpdate( n+m+e, n+m+e, -s(e)/z(e) );
    J.ProcessQueues();
}

template<typename Real>
void KKTRHS
( const Matrix<Real>& rc,
  const Matrix<Real>& rb,
  const Matrix<Real>& rh,
  const Matrix<Real>& rmu,
  const Matrix<Real>& z,
        Matrix<Real>& d )
{
    EL_DEBUG_CSE
    const Int n = rc.Height();
    const Int m = rb.Height();
    const Int k = rh.Height();
    Zeros( d, n+m+k, 1 );

    auto dx = d(IR(0,n),ALL);
    dx = rc;
    dx *= -1;

    auto dy = d(IR(n,n+m),ALL);
    dy = rb;
    dy *= -1;

    auto dz = d(IR(n+m,n+m+k),ALL);
    dz = rmu;
    DiagonalSolve( LEFT, NORMAL, z, dz );
    dz -= rh;
}

template<typename Real>
void ExpandCoreSolution
( Int m, Int n, Int k,
  const Matrix<Real>& d,
        Matrix<Real>& dx,
        Matrix<Real>& dy,
        Matrix<Real>& dz )
{
    EL_DEBUG_CSE
    if( d.Height() != n+m+k || d.Width() != 1 )
        LogicError("Right-hand side was the wrong size");
    dx = d(IR(0,  n    ),ALL);
    dy = d(IR(n,  n+m  ),ALL);
    dz = d(IR(n+m,n+m+k),ALL);
}

template<typename Real>
void Initialize
( const SparseMatrix<Real>& JStatic,
  const Matrix<Real>& regLarge,
  const Matrix<Real>& b,
  const Matrix<Real>& c,
  const SparseMatrix<Real>& G,
  const Matrix<Real>& h,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
        SparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit,
  bool dualInit,
  bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE
    const Int m = b.Height();
    const Int n = c.Height();
    const Int k = h.Height();
    if( primalInit )
    {
        if( x.Height() != n || x.Width() != 1 )
            LogicError("x was of the wrong size");
        if( s.Height() != k || s.Width() != 1 )
            LogicError("s was of the wrong size");
    }
    if( dualInit )
    {
        if( y.Height() != m || y.Width() != 1 )
            LogicError("y was of the wrong size");
        if( z.Height() != k || z.Width() != 1 )
            LogicError("z was of the wrong size");
    }
    if( primalInit && dualInit )
    {
        // TODO(poulson): Perform a consistency check
        return;
    }

    // Form the KKT matrix
    // ===================
    auto JOrig = JStatic;
    JOrig.FreezeSparsity();
    Matrix<Real> ones;
    Ones( ones, k, 1 );
    FinishKKT( m, n, ones, ones, JOrig );
    auto J = JOrig;
    J.FreezeSparsity();
    UpdateRealPartOfDiagonal( J, Real(1), regLarge );

    // Analyze the sparsity pattern of the KKT system
    // ==============================================
    const bool hermitian = true;
    const BisectCtrl bisectCtrl;
    sparseLDLFact.Initialize( J, hermitian, bisectCtrl );

    // (Approximately) factor the KKT matrix
    // =====================================
    sparseLDLFact.Factor();

    // Compute the proposed step from the KKT system
    // ---------------------------------------------
    Matrix<Real> rc, rb, rh, rmu, u, d;
    Zeros( rmu, k, 1 );
    if( !primalInit )
    {
        // Minimize || G x - h ||^2, s.t. A x = b  by solving
        //
        //    | Q A^T G^T | |  x |   | 0 |
        //    | A  0   0  | |  u | = | b |,
        //    | G  0  -I  | | -s |   | h |
        //
        //   where 'u' is an unused dummy variable.
        Zeros( rc, n, 1 );
        rb = b;
        rb *= -1;
        rh = h;
        /*
        Output("Forcefully clipping 'h'");
        LowerClip( rh, Real(1e-2) );
        UpperClip( rh, Real(1e3) );
        */
        rh *= -1;
        KKTRHS( rc, rb, rh, rmu, ones, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, regLarge, sparseLDLFact, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );

        ExpandCoreSolution( m, n, k, d, x, u, s );
        // s := h - G x
        s = h;
        Multiply( NORMAL, Real(-1), G, x, Real(1), s );
    }
    if( !dualInit )
    {
        // Minimize || z ||^2, s.t. A^T y + G^T z + c in range(Q) by solving
        //
        //    | Q A^T G^T | | u |   | -c |
        //    | A  0   0  | | y | = |  0 |,
        //    | G  0  -I  | | z |   |  0 |
        //
        //    where 'u' is an unused dummy variable.
        rc = c;
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        KKTRHS( rc, rb, rh, rmu, ones, d );

        reg_ldl::RegularizedSolveAfter
        ( JOrig, regLarge, sparseLDLFact, d,
          solveCtrl.relTol, solveCtrl.maxRefineIts, solveCtrl.progress );

        ExpandCoreSolution( m, n, k, d, u, y, z );
    }

    /*
    const Real eps = limits::Epsilon<Real>();
    const Real sNorm = Nrm2( s );
    const Real zNorm = Nrm2( z );
    const Real gammaPrimal = Sqrt(eps)*Max(sNorm,Real(1));
    const Real gammaDual   = Sqrt(eps)*Max(zNorm,Real(1));
    if( standardShift )
    {
        // alpha_p := min { alpha : s + alpha*e >= 0 }
        // -------------------------------------------
        const auto sMinPair = VectorMinLoc( s );
        const Real alphaPrimal = -sMinPair.value;
        if( alphaPrimal >= Real(0) && primalInit )
            RuntimeError("initialized s was non-positive");

        // alpha_d := min { alpha : z + alpha*e >= 0 }
        // -------------------------------------------
        const auto zMinPair = VectorMinLoc( z );
        const Real alphaDual = -zMinPair.value;
        if( alphaDual >= Real(0) && dualInit )
            RuntimeError("initialized z was non-positive");

        if( alphaPrimal >= -gammaPrimal )
            Shift( s, alphaPrimal+1 );
        if( alphaDual >= -gammaDual )
            Shift( z, alphaDual+1 );
    }
    else
    {
        LowerClip( s, gammaPrimal );
        LowerClip( z, gammaDual   );
    }
    */
    // Use Mehrotra's approach:
    //
    //   Mehrotra S. (1992): On the Implementation of a Primal-Dual 
    //     Interior Point Method, SIAM Journal on Optimization 2, 
    //     No 4, pp. 575-601.
    //
    //  Mehrotra S. (1991): Higher Order Methods and their Performance,
    //     Technical Report 90-16R1, Department of Industrial Engineering
    //     and Management Sciences, Northwestern University, Evanston,
    //     Illinois 60208-3119, U.S.A.
    Output("Forcing Mehrotra's initialization strategy");
    Real sMin = MaxNorm( s );
    Real zMin = MaxNorm( z );
    for( Int i=0; i<k; ++i )
    {
        sMin = Min( sMin, s(i) );
        zMin = Min( zMin, z(i) );
    }
    const Real deltaPrimal = Max( Real(-1.5)*sMin, Real(0) );
    const Real deltaDual = Max( Real(-1.5)*zMin, Real(0) );
    Output("deltaPrimal=",deltaPrimal,", deltaDual=",deltaDual);
    Real simpleGap = Real(0);
    for( Int i=0; i<k; ++i )
        simpleGap += (s(i)+deltaPrimal)*(z(i)+deltaDual);
    Output("simpleGap=",simpleGap);
    Real sSum = Real(0);
    Real zSum = Real(0);
    for( Int i=0; i<k; ++i )
    {
        sSum += s(i) + deltaPrimal;
        zSum += z(i) + deltaDual;
    }
    Output("sSum=",sSum,", zSum=",zSum);
    for( Int i=0; i<k; ++i )
    {
        s(i) += deltaPrimal + simpleGap/(2*zSum);
        z(i) += deltaDual + simpleGap/(2*sSum);
    }
}

template<typename Real>
void Initialize
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const SparseMatrix<Real>& JStatic,
  const Matrix<Real>& regTmp,
        SparseLDLFactorization<Real>& sparseLDLFact,
  bool primalInit,
  bool dualInit,
  bool standardShift,
  const RegSolveCtrl<Real>& solveCtrl )
{
    EL_DEBUG_CSE
    Initialize
    ( JStatic, regTmp,
      problem.b, problem.c, problem.G, problem.h,
      solution.x, solution.y, solution.z, solution.s,
      sparseLDLFact,
      primalInit, dualInit, standardShift, solveCtrl );
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
void ExpandSolution
( Int m, Int n,
  const Matrix<Real>& d,
  const Matrix<Real>& rmu,
  const Matrix<Real>& s,
  const Matrix<Real>& zReg,
        Matrix<Real>& dx,
        Matrix<Real>& dy,
        Matrix<Real>& dz,
        Matrix<Real>& ds )
{
    EL_DEBUG_CSE
    const Int k = s.Height();
    ExpandCoreSolution( m, n, k, d, dx, dy, dz );
    // ds := - zReg <> ( rmu + s o dz )
    // ================================
    ds = dz;
    DiagonalScale( LEFT, NORMAL, s, ds );
    ds += rmu;
    DiagonalSolve( LEFT, NORMAL, zReg, ds );
    ds *= -1;
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
( SparseAffineLPEquilibration<Real>& equilibration,
  AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  AffineLPState<Real>& state,
  Matrix<Real>& regSmall,
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
( const AffineLPSolution<Matrix<Real>>& solution, bool print )
{
    EL_DEBUG_CSE
    const Int sNumNonPos = pos_orth::NumOutside( solution.s );
    const Int zNumNonPos = pos_orth::NumOutside( solution.z );
    if( sNumNonPos > 0 || zNumNonPos > 0 )
        LogicError
        (sNumNonPos," entries of s were nonpositive and ",
         zNumNonPos," entries of z were nonpositive");

    if( print )
    {
        const Int k = solution.s.Height();
        const Real sMax = MaxNorm( solution.s );
        const Real zMax = MaxNorm( solution.z );
        Real sMin = sMax;
        for( Int i=0; i<k; ++i )
            sMin = Min( sMin, solution.s(i) );
        Real zMin = zMax;
        for( Int i=0; i<k; ++i )
            zMin = Min( zMin, solution.z(i) );
        Output("  sMin = ",sMin,", zMin = ",zMin);
        Output("  sMax = ",sMax,", zMax = ",zMax);
    }
}

template<typename Real>
void UpdateCombinedResidualUsingAffine
( const Real& mu,
  const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
  const AffineLPState<Real>& state,
  const AffineLPSolution<Matrix<Real>>& solution,
  const AffineLPSolution<Matrix<Real>>& affineCorrection,
        AffineLPResidual<Matrix<Real>>& residual,
  const Real& complementRatio,
  bool compositeNewton,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int k = solution.s.Height();
    Shift( residual.dualConic, -mu );
    const bool infeasibleCompositeNewton = false;
    if( compositeNewton && infeasibleCompositeNewton )
    {
        Matrix<Real> input, update;

        // residual.primalEquality +=
        //   A (x + dxAff) - b - yRegSmall (y + dyAff) - yRegLarge dyAff
        Zeros( update, problem.A.Height(), 1 );
        input = solution.x;
        input += affineCorrection.x;
        Multiply( NORMAL, Real(1), problem.A, input, Real(1), update );
        Axpy( Real(-1), problem.b, update );
        input = solution.y;
        input += affineCorrection.y;
        Axpy( -state.yRegSmall, input, update );
        Axpy( -state.yRegLarge, affineCorrection.y, update );
        residual.primalEquality += update;

        // residual.dualEquality +=
        //   c + A^T (y + dyAff) + G^T (z + dzAff) + xRegSmall (x + dxAff) +
        //   xRegLarge dxAff
        Zeros( update, problem.A.Width(), 1 );
        input += problem.c;
        input = solution.y;
        input += affineCorrection.y;
        Multiply( TRANSPOSE, Real(1), problem.A, input, Real(1), update );
        input = solution.z;
        input += affineCorrection.z;
        Multiply( TRANSPOSE, Real(1), problem.G, input, Real(1), update );
        input = solution.x;
        input += affineCorrection.x;
        Axpy( state.xRegSmall, input, update );
        Axpy( state.xRegLarge, affineCorrection.x, update );
        residual.dualEquality += update;

        // residual.primalConic +=
        //   G (x + dxAff) + (s + dsAff) - h - zRegSmall (z + dzAff) -
        //   zRegLarge dzAff
        Zeros( update, problem.G.Height(), 1 );
        input = solution.x;
        input += affineCorrection.x;
        Multiply( NORMAL, Real(1), problem.G, input, Real(1), update );
        input = solution.s;
        input += affineCorrection.s;
        update += input;
        Axpy( Real(-1), problem.h, update );
        input = solution.z;
        input += affineCorrection.z;
        Axpy( -state.zRegSmall, input, update );
        Axpy( -state.zRegLarge, affineCorrection.z, update );
        residual.primalConic += update;

        // residual.dualConic += (s + dsAff) o (z + dzAff)
        // TODO(poulson): Incorporate zPivot instead of z?
        input = solution.s;
        input += affineCorrection.s;
        update = solution.z;
        update += affineCorrection.z;
        DiagonalScale( LEFT, NORMAL, input, update );
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
  const Real& lowerComplementRatio,
  const Real& upperComplementRatio,
  bool largeComplementRatio,
  const AffineLPSolution<Matrix<Real>>& solution,
        AffineLPResidual<Matrix<Real>>& residual,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE 
    const Int k = solution.s.Height();
    if( ctrl.softDualityTargets )
    {
        // Attempt to correct s o z entrywise into
        // [lowerTargetRatio,upperTargetRatio]*sigma*mu.
        /*
        Real lowerTargetRatio = 
          Pow(complementRatio,ctrl.lowerTargetRatioLogCompRatio);
        Real upperTargetRatio =
          Pow(complementRatio,ctrl.upperTargetRatioLogCompRatio);
        */
        /*
        lowerTargetRatio =
          Max( lowerTargetRatio,
            Pow(ctrl.maxComplementRatio,ctrl.lowerTargetRatioLogMaxCompRatio) );
        upperTargetRatio =
          Min( upperTargetRatio,
            Pow(ctrl.maxComplementRatio,ctrl.upperTargetRatioLogMaxCompRatio) );
        */
        Output
        ("lowerComplementRatio=",lowerComplementRatio,
         ", upperComplementRatio=",upperComplementRatio);
        Real lowerTargetRatio, upperTargetRatio;
        if( lowerComplementRatio < Real(0.01) )
        {
            lowerTargetRatio = lowerComplementRatio*Real(10.);
            upperTargetRatio = 10*upperComplementRatio;
        }
        else
        {
            lowerTargetRatio = Min( lowerComplementRatio*Real(10), Real(0.1) );
            upperTargetRatio = Max( upperComplementRatio/Real(10), Real(10.) );
        }
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

    if( ctrl.print )
    {
        auto dx = packedVector( IR(0,n), ALL );
        auto dy = packedVector( IR(n,n+m), ALL );
        auto dz = packedVector( IR(n+m,END), ALL );
        const Real dxScaledNrm2 = FrobeniusNorm( dx );
        const Real dyScaledNrm2 = FrobeniusNorm( dy );
        const Real dzScaledNrm2 = FrobeniusNorm( dz );
        DiagonalScale( LEFT, NORMAL, diagonalScale, packedVector );
        const Real dxNrm2 = FrobeniusNorm( dx );
        const Real dyNrm2 = FrobeniusNorm( dy );
        const Real dzNrm2 = FrobeniusNorm( dz );
        Output("  || dxScaled ||_2 = ",dxScaledNrm2);
        Output("  || dyScaled ||_2 = ",dyScaledNrm2);
        Output("  || dzScaled ||_2 = ",dzScaledNrm2);
        Output("  || dx       ||_2 = ",dxNrm2);
        Output("  || dy       ||_2 = ",dyNrm2);
        Output("  || dz       ||_2 = ",dzNrm2);
    }
    else
    {
        DiagonalScale( LEFT, NORMAL, diagonalScale, packedVector );
    }
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
bool CheckConvergence
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

    // primalEqualityResid := A x - b
    // ------------------------------
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

    // dualEqualityResid := c + A^T y + G^T z
    // --------------------------------------
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

    // primalConicResid := G x + s - h
    // -------------------------------
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

    // primalEqualityResid := A x - b - yRegSmall y
    // --------------------------------------------
    Axpy( -state.yRegSmall, solution.y, residual.primalEquality );

    // dualEqualityResid := c + A^T y + G^T z + xRegSmall x
    // ----------------------------------------------------
    Axpy( state.xRegSmall, solution.x, residual.dualEquality );

    // primalConicResid := G x + s - h - zRegSmall z
    // ---------------------------------------------
    Axpy( -state.zRegSmall, solution.z, residual.primalConic );

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

    state.approxFeasible =
      state.infeasError <= Pow(epsilon,ctrl.infeasibilityTolLogEps);

    state.metTolerances =
      state.approxFeasible &&
      state.relCompGap <=
        Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
      state.relObjGap <= Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);

    state.metTolerancesOrig =
      state.metTolerances &&
      state.relCompGapOrig <=
        Pow(epsilon,ctrl.relativeComplementarityGapTolLogEps) &&
      state.relObjGapOrig <=
        Pow(epsilon,ctrl.relativeObjectiveGapTolLogEps);

    const Real avgCompRatio = Dot( solution.s, solution.z ) / Real(k);
    state.lowerComplementRatio = Real(1);
    state.upperComplementRatio = Real(1);
    for( Int i=0; i<k; ++i )
    {
        const Real dualProd = solution.s(i) * solution.z(i);
        if( dualProd > avgCompRatio )
            state.upperComplementRatio =
              Max( state.upperComplementRatio, dualProd / avgCompRatio );
        else
            state.lowerComplementRatio =
              Min( state.lowerComplementRatio, dualProd / avgCompRatio );
    }
    state.complementRatio =
      state.upperComplementRatio / state.lowerComplementRatio;

    bool converged = false;
    if( state.metTolerancesOrig )
    {
        if( state.dimacsErrorOrig >=
            ctrl.minDimacsDecreaseRatio*state.dimacsErrorOrigOld )
        {
            // We have met the tolerances and progress in the last iteration
            // was not significant.
            converged = true;
        }
        else if( state.numIterations == ctrl.maxIts )
        {
            // We have hit the iteration limit but can declare success.
            converged = true;
        }
    }
    else if( state.metTolerances && !state.metTolerancesOrig )
    {
        state.backingOffEquilibration = true;
        if( ctrl.print )
            Output("  will back off from equilibration in next step");
    }
    else
    {
        state.backingOffEquilibration = false;
    }

    if( ctrl.print )
    {
        Output
        ("iter ",state.numIterations,":\n",Indent(),
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
        const Real xOrigNrm2 = FrobeniusNorm( origSolution.x );
        const Real yOrigNrm2 = FrobeniusNorm( origSolution.y );
        const Real zOrigNrm2 = FrobeniusNorm( origSolution.z );
        const Real sOrigNrm2 = FrobeniusNorm( origSolution.s );
        const Real xNrm2 = FrobeniusNorm( solution.x );
        const Real yNrm2 = FrobeniusNorm( solution.y );
        const Real zNrm2 = FrobeniusNorm( solution.z );
        const Real sNrm2 = FrobeniusNorm( solution.s );
        Output("  || xOrig ||_2 = ",xOrigNrm2);
        Output("  || yOrig ||_2 = ",yOrigNrm2);
        Output("  || zOrig ||_2 = ",zOrigNrm2);
        Output("  || sOrig ||_2 = ",sOrigNrm2);
        Output("  || x     ||_2 = ",xNrm2);
        Output("  || y     ||_2 = ",yNrm2);
        Output("  || z     ||_2 = ",zNrm2);
        Output("  || s     ||_2 = ",sNrm2);
    }

    return converged;
}

template<typename Real>
std::pair<Real,Real>
ComputePrimalDualStepSizes
( const AffineLPSolution<Matrix<Real>>& solution,
  const AffineLPSolution<Matrix<Real>>& correction,
  const Real& upperBound,
  const Real& maxStepRatio,
  const Real& maxPrimalStepImbalance,
  const Real& maxDualStepImbalance )
{
    EL_DEBUG_CSE
    Real primalStep = pos_orth::MaxStep( solution.s, correction.s, upperBound );
    Real dualStep = pos_orth::MaxStep( solution.z, correction.z, upperBound );
    primalStep = Min( primalStep, maxPrimalStepImbalance*dualStep );
    dualStep = Min( dualStep, maxDualStepImbalance*primalStep );
    primalStep = Min( maxStepRatio*primalStep, Real(1) );
    dualStep = Min( maxStepRatio*dualStep, Real(1) );
    return std::make_pair( primalStep, dualStep );
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
    /*
    Ones( rowScaleA, m, 1 );
    Ones( rowScaleG, k, 1 );
    Ones( colScale, n, 1 );
    */
    StackedGeomEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      rowScaleA, rowScaleG, colScale,
      ctrl.print );
    /*
    StackedRuizEquil
    ( equilibratedProblem.A,
      equilibratedProblem.G,
      rowScaleA, rowScaleG, colScale,
      ctrl.print );
    */
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
void DecideBarrierAndStepStrategy
(       AffineLPState<Real>& state,
  const AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    const Int k = solution.s.Height();
    const Real gap = Dot( solution.s, solution.z );
    const Real barrierClassical = gap / k;
    //const Real barrierMedian = MedianBarrier( solution.s, solution.z );
    if( ctrl.print )
        Output
        ("  barrierClassical=",barrierClassical);
        /*
        ("  barrierMedian = ",barrierMedian,
         ", barrierClassical=",barrierClassical);
        */

    // Decide whether or not we will attempt a predictor-corrector primarily
    // based upon the complementarity ratio.
    // TODO(poulson): Generalize these bounds.
    state.largeComplementRatio =
      state.lowerComplementRatio < Real(0.01) ||
      state.upperComplementRatio > Real(1.e4);
    // It has been observed that models such as wood1w can continually try
    // to drive down the barrier parameter despite relative infeasibility
    // (e.g., a relative complementarity gap of 10^-14 while the relative
    // dual infeasibility is on the order of 10^-7).
    const Real relativeComplementarityLowerBound =
      Min(Pow(state.infeasError,Real(1.)),Real(1e-3));
    const bool relCompGapIsTooSmall =
      state.relCompGap < relativeComplementarityLowerBound;
    const bool gapIsTooSmall = gap < state.infeasError;

    state.centeringStep =
      state.largeComplementRatio || relCompGapIsTooSmall || gapIsTooSmall;

    if( state.centeringStep )
    {
        if( ctrl.print )
            Output("  will take a centering step");
        state.barrier = barrierClassical;
        if( ctrl.print )
            Output
            ("  using barrier = barrierClassical = ",barrierClassical);
        state.lastStepForcedBarrierIncrease = false;
    }
    else
    {
        if( ctrl.print )
            Output("  will take a predictor-corrector step.");
        if( state.backingOffEquilibration )
        {
            state.barrier = barrierClassical;
            if( ctrl.print )
                Output
                ("  backing off equilibration, using barrierClassical = ",
                 barrierClassical);
        }
        else
        {
            state.barrier = barrierClassical;
            if( ctrl.print )
                Output("  using barrierClassical = ",barrierClassical);
        }
        state.lastStepForcedBarrierIncrease = false;
    }
}

template<typename Real>
IPMInfo<Real> EquilibratedIPM
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
    const Real epsilon = limits::Epsilon<Real>();
    const Real regIncreaseFactor = Pow(epsilon,ctrl.regIncreaseFactorLogEps);
    if( regIncreaseFactor <= Real(1) )
        LogicError("Regularization increase factor must be > 1");
    const Real minStepSize = Pow(limits::Epsilon<Real>(),Real(0.33));
    Real maxPredictorCorrectorStepImbalance =
      ctrl.forceSameStep ? Real(1) : Real(10);

    /*
    Print( problem.c, "cEquil" );
    Print( problem.A, "AEquil" );
    Print( problem.b, "bEquil" ); 
    Print( problem.G, "GEquil" );
    Print( problem.h, "hEquil" );
    */

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

    Matrix<Real> regSmall, regSmallScaled;
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
      //Sqrt(state.xRegSmall), Sqrt(state.yRegSmall), Sqrt(state.zRegSmall),
      Real(0), Real(0), Real(0),
      JStatic, false );
    JStatic.FreezeSparsity();

    SparseLDLFactorization<Real> sparseLDLFact;
    const bool hermitian = true;
    const BisectCtrl bisectCtrl;
    sparseLDLFact.Initialize( JStatic, hermitian, bisectCtrl );
    Ones( solution.x, n, 1 );
    solution.tau = Real(1);
    Zeros( solution.y, m, 1 );
    Ones( solution.z, k, 1 );
    solution.kappa = Real(1);
    /*
    Initialize
    ( problem, solution, JStatic, regLarge,
      sparseLDLFact,
      ctrl.primalInit, ctrl.dualInit, ctrl.standardInitShift, ctrl.solveCtrl );
    */

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
            if( state.numIterations == 0 && ctrl.primalInit && ctrl.dualInit )
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

    Int numTinySteps = 0;
    const Int indent = PushIndent();
    auto lastSolution( solution );
    Real dxCombinedNrm2Last=Real(1),
         dyCombinedNrm2Last=Real(1),
         dzCombinedNrm2Last=Real(1);
    for( ; state.numIterations<=ctrl.maxIts;
         ++state.numIterations,
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
        ( equilibration, problem, state, regSmall, solution, ctrl );
        AssertPositiveOrthantMembership( solution, ctrl.print );

        // Manage the duality products
        {
            const Real barrierClassical = Dot(solution.s,solution.z) / k;
            const Real primalObj = Dot(problem.c,solution.x);
            const Real dualObj = -Dot(problem.b,solution.y) -
              Dot(problem.h,solution.z);
            const Real dualityGap = primalObj - dualObj;
            const Real relDualityGap = Abs(dualityGap) /
              (Abs(primalObj) + Abs(dualObj) + Real(1));

            Real minDualProd = Real(1e-6)*barrierClassical;
            //minDualProd = Max( minDualProd, Real(1e-4)*Abs(dualityGap)/k );
            if( state.numIterations == 0 )
                minDualProd = Min( minDualProd, Real(1) );

            Real minPrimalConicValue = Real(0);
            Real minDualConicValue = Real(0);
            /*
            Real minPrimalConicValue = Real(1e-4)*relDualityGap;
            Real minDualConicValue = Real(1e-4)*relDualityGap;
            */
            minDualConicValue = Min( minDualConicValue, Real(1e-6) );

            Real maxProd = 0;
            for( Int i=0; i<k; ++i )
                maxProd = Max( maxProd, solution.s(i)*solution.z(i) );
            Real minProd = maxProd;
            for( Int i=0; i<k; ++i )
                minProd = Min( minProd, solution.s(i)*solution.z(i) );
            Output
            ("  Before: minProd=",minProd,", avgProd=",barrierClassical,
             ", maxProd=",maxProd);
            for( Int i=0; i<k; ++i )
            {
                const Real prod = solution.s(i)*solution.z(i);
                if( prod == minProd || prod == maxProd )
                    Output
                    ("  Extremal s(",i,")=",solution.s(i)," * z(",i,")=",
                     solution.z(i)," = ",prod);
            }

            for( Int i=0; i<k; ++i )
            {
                const Real dualProd = solution.s(i) * solution.z(i);
                if( dualProd < minDualProd )
                {
                    if( solution.s(i) <= solution.z(i) )
                    {
                        Output
                        ("  increasing s(",i,")=",solution.s(i)," to ",
                         minDualProd/solution.z(i));
                        solution.s(i) = minDualProd / solution.z(i);
                    }
                    else
                    {
                        Output
                        ("  increasing z(",i,")=",solution.z(i)," to ",
                         minDualProd/solution.s(i));
                        solution.z(i) = minDualProd / solution.s(i);
                    }
                }
                if( solution.s(i) < minPrimalConicValue )
                {
                    Output
                    ("  increasing s(",i,")=",solution.s(i)," to ",
                     minPrimalConicValue);
                    solution.s(i) = minPrimalConicValue;
                }
                if( solution.z(i) < minDualConicValue )
                {
                    Output
                    ("  increasing z(",i,")=",solution.z(i)," to ",
                     minDualConicValue);
                    solution.z(i) = minDualConicValue;
                }
            }

            maxProd = 0;
            for( Int i=0; i<k; ++i )
                maxProd = Max( maxProd, solution.s(i)*solution.z(i) );
            minProd = maxProd;
            for( Int i=0; i<k; ++i )
                minProd = Min( minProd, solution.s(i)*solution.z(i) );
            Output
            ("  Before: minProd=",minProd,", avgProd=",barrierClassical,
             ", maxProd=",maxProd);
            for( Int i=0; i<k; ++i )
            {
                const Real prod = solution.s(i)*solution.z(i);
                if( prod == minProd || prod == maxProd )
                    Output
                    ("  Extremal s(",i,")=",solution.s(i)," * z(",i,")=",
                     solution.z(i)," = ",prod);
            }
        }

        const Real maxAbsObj = Max( Abs(state.primalObj), Abs(state.dualObj) );
        Real regSmallValue;
        const Real barrierToRegSmallConstant = Pow(epsilon,Real(0));
        const Real barrierToRegSmallExponent = Real(0.8);
        regSmallValue = barrierToRegSmallConstant *
          Pow(state.barrier,barrierToRegSmallExponent);
        regSmallValue = Min( regSmallValue, Pow(epsilon,Real(0.8)) );

        /*
        Print( solution.x, "x"+std::to_string(state.numIterations) );
        Print( solution.s, "s"+std::to_string(state.numIterations) );
        Print( solution.y, "y"+std::to_string(state.numIterations) );
        Print( solution.z, "z"+std::to_string(state.numIterations) );
        */

        Output("barrier=",state.barrier,", regSmallValue=",regSmallValue);
        state.xRegSmall = state.yRegSmall = state.zRegSmall = regSmallValue;
        auto regSmall_x = regSmall( IR(0,n), ALL );
        auto regSmall_y = regSmall( IR(n,n+m), ALL );
        auto regSmall_z = regSmall( IR(n+m,END), ALL );
        Fill( regSmall_x, state.xRegSmall );
        Fill( regSmall_y, -state.yRegSmall );
        Fill( regSmall_z, -state.zRegSmall );

        Matrix<Real> zPivot, zPerturb;
        const bool converged =
          CheckConvergence
          ( origProblem, equilibration, problem, solution,
            zPivot, zPerturb, residual, state, ctrl );
        if( converged )
        {
            break;
        }
        else if( state.numIterations == ctrl.maxIts )
        {
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving tolerances");
        }
        else if( !state.lastStepForcedBarrierIncrease &&
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

        // Factor the (rescaled) KKT system
        // ================================
        Ones( diagonalScale, n+m+k, 1 );
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        auto zDoublePivot( zPivot );
        zDoublePivot += zPerturb;
        FinishKKT( m, n, solution.s, zDoublePivot, JOrig );
        UpdateDiagonal( JOrig, Real(1), regSmall );
        JOrigScaled = JOrig;
        JOrigScaled.FreezeSparsity();
        regSmallScaled = regSmall;
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
        DecideBarrierAndStepStrategy( state, solution, ctrl );

        Real sigma;
        if( state.centeringStep )
        {
            sigma = 1;
            if( ctrl.print )
                Output("  freezing sigma at one");

            // Solve for the combined direction
            // ================================
            UpdateCombinedDualConicResidualWithoutAffine
            ( sigma, state.barrier,
              state.lowerComplementRatio, state.upperComplementRatio,
              state.largeComplementRatio,
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
            const Real dxNrm2 = FrobeniusNorm( correction.x );
            const Real dsNrm2 = FrobeniusNorm( correction.s );
            const Real dyNrm2 = FrobeniusNorm( correction.y );
            const Real dzNrm2 = FrobeniusNorm( correction.z );
            Output("  || dx ||_2 = ",dxNrm2,", || ds ||_2 = ",dsNrm2);
            Output("  || dy ||_2 = ",dyNrm2,", || dz ||_2 = ",dzNrm2);
        }
        else
        {
            // Compute the affine search direction
            // ===================================

            // Set the dualConic residual to a small multiple of the current
            // median barrier parameter rather than the value of zero?
            const Real affineBarrierMultiple = Real(1e-3);
            Shift( residual.dualConic, -affineBarrierMultiple*state.barrier );
            Output("  affine barrier target: ",affineBarrierMultiple*state.barrier);

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

            // Undo the shift...
            Shift( residual.dualConic, affineBarrierMultiple*state.barrier );

            // Compute a centrality parameter
            // ==============================
            auto affineStep =
              ComputePrimalDualStepSizes
              ( solution, affineCorrection,
                Real(1) /* upperBound */,
                Real(1) /* maxStepRatio */,
                Real(1) /* maxPrimalStepImbalance */,
                Real(1) /* maxDualStepImbalance */ );
            if( ctrl.print )
            {
                Output("  affine primal step size: ",affineStep.first);
                Output("  affine dual   step size: ",affineStep.second);
            }
            const Real minAffineStep =
              Min( affineStep.first, affineStep.second );
            // TODO(poulson): Experiment with several centrality choices?
            auto sAffine( solution.s );
            Axpy( affineStep.first, affineCorrection.s, sAffine );
            auto zAffine( solution.z );
            Axpy( affineStep.second, affineCorrection.z, zAffine );
            //const Real barrierAff = MedianBarrier( sAffine, zAffine );
            const Real barrierAff = Dot( sAffine, zAffine ) / Real(k);
            if( ctrl.print )
                Output
                ("  barrierAff = ",barrierAff,", barrier = ",state.barrier);
            Real sigmaAff = Min(Pow(barrierAff/state.barrier,Real(2)),Real(1));
            Output("  sigmaAff before: ",sigmaAff);
            Real weightedLength = Real(0);
            for( Int i=0; i<k; ++i )
            {
                const Real theta = solution.z(i) / solution.s(i);
                weightedLength +=
                  affineCorrection.s(i)*affineCorrection.s(i)/theta +
                  affineCorrection.z(i)*affineCorrection.z(i)*theta;
            }
            if( ctrl.print )
                Output
                ("  || ds ||^2_{inv(theta)} + || dz ||^2_{theta} = ",
                 weightedLength);
            if( weightedLength < state.barrier )
            {
                if( minAffineStep >= Real(0.2) )
                    sigmaAff = Min( sigmaAff, Real(0.1) );
                else if( minAffineStep >= Real(0.1) )
                    sigmaAff = Min( sigmaAff, Real(0.2) );
                else
                    sigmaAff = Min( sigmaAff, Real(0.33) );
            }
            if( ctrl.print )
                Output("  sigmaAff=",sigmaAff);

            // Solve for the combined direction
            // ================================
            UpdateCombinedResidualUsingAffine
            ( sigmaAff*barrierAff, problem, state, solution, affineCorrection,
              residual, state.complementRatio, ctrl.compositeNewton, ctrl );
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
        Real primalCorrectionNrm2, dualCorrectionNrm2;
        {
            const Real xCorrectionNrm2 = FrobeniusNorm(correction.x);
            const Real sCorrectionNrm2 = FrobeniusNorm(correction.s);
            const Real yCorrectionNrm2 = FrobeniusNorm(correction.y);
            const Real zCorrectionNrm2 = FrobeniusNorm(correction.z);
            // TODO(poulson): Use the stable norm update routines instead.
            primalCorrectionNrm2 =
              Sqrt(xCorrectionNrm2*xCorrectionNrm2 +
                   sCorrectionNrm2*sCorrectionNrm2);
            dualCorrectionNrm2 =
              Sqrt(yCorrectionNrm2*yCorrectionNrm2 +
                   zCorrectionNrm2*zCorrectionNrm2);
        }
        Real maxPrimalStepImbalance, maxDualStepImbalance;
        if( state.centeringStep || !state.approxFeasible )
        {
            maxPrimalStepImbalance = Real(1);
            maxDualStepImbalance = Real(1);
        }
        else
        {
            const Real maxRatio = Pow(epsilon,Real(-0.1));

            if( primalCorrectionNrm2 > maxRatio*dualCorrectionNrm2 )
                maxDualStepImbalance = Real(1);
            else
                maxDualStepImbalance = maxPredictorCorrectorStepImbalance;

            if( dualCorrectionNrm2 > maxRatio*primalCorrectionNrm2 )
                maxPrimalStepImbalance = Real(1);
            else
                maxPrimalStepImbalance = maxPredictorCorrectorStepImbalance;
        }
        auto combinedStep =
          ComputePrimalDualStepSizes
          ( solution, correction,
            Real(1)/ctrl.maxStepRatio,
            ctrl.maxStepRatio,
            maxPrimalStepImbalance, maxDualStepImbalance );
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
                if( numTinySteps > 10 )
                {
                    Print( problem.c, "c" );
                    Print( problem.A, "A" );
                    Print( problem.b, "b" );
                    Print( problem.G, "G" );
                    Print( problem.h, "h" );
                    Print( solution.x, "x" );
                    Print( solution.s, "s" );
                    Print( solution.y, "y" );
                    Print( solution.z, "z" );
                    RuntimeError("Attempted too many tiny steps");
                }
                if( ctrl.print )
                  Output
                  ("Attempted step sizes less than ",minStepSize,
                   ", so increasing regularization");
                maxPredictorCorrectorStepImbalance = Real(1);
                increaseRegularization();
                ++numTinySteps;
                continue;
            }
        }
        if( ctrl.print )
            Output("");
    }
    SetIndent( indent );

    IPMInfo<Real> info;
    info.primalObjective = state.primalObjOrig;
    info.dualObjective = state.dualObjOrig;
    info.infeasibilityError = state.infeasError;
    info.relativeObjectiveGap = state.relObjGap;
    info.relativeComplementarityGap = state.relCompGap;
    info.numIterations = state.numIterations;
    return info;
}

template<typename Real>
IPMInfo<Real> IPM
( const AffineLPProblem<SparseMatrix<Real>,Matrix<Real>>& problem,
        AffineLPSolution<Matrix<Real>>& solution,
  const IPMCtrl<Real>& ctrl )
{
    EL_DEBUG_CSE
    IPMInfo<Real> info;
    if( ctrl.outerEquil )
    {
        AffineLPProblem<SparseMatrix<Real>,Matrix<Real>> equilibratedProblem;
        AffineLPSolution<Matrix<Real>> equilibratedSolution;
        SparseAffineLPEquilibration<Real> equilibration;
        Equilibrate
        ( problem, solution,
          equilibratedProblem, equilibratedSolution,
          equilibration, ctrl );
        try
        {
            info = EquilibratedIPM
            ( problem, equilibration,
              equilibratedProblem, equilibratedSolution, ctrl );
        }
        catch( const std::exception& except )
        {
            // Preserve the current approximate solution.
            UndoEquilibration( equilibratedSolution, equilibration, solution );
            throw except;
        }
        UndoEquilibration( equilibratedSolution, equilibration, solution );

        // Refinement seems to lead to very challenging problems (e.g., for
        // 80bau3b).
        /*
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
        */
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
        info = EquilibratedIPM
        ( problem, equilibration,
          equilibratedProblem, equilibratedSolution, ctrl );
        UndoEquilibration( equilibratedSolution, equilibration, solution );
    }
    return info;
}

} // namespace affine
} // namespace lp_hsd

template<typename Real>
El::LPInfo<Real> DenseLoadAndSolve
( const std::string& filename,
  bool metadataSummary,
  bool print,
  bool progress,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool forceSameStep,
  bool compositeNewton,
  El::Int maxIter )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::Matrix<",El::TypeName<Real>(),">");
    El::Timer timer;
    const bool compressed = false;
    const bool minimize = true;
    const bool keepNonnegativeWithZeroUpperBound=true;

    timer.Start();
    El::AffineLPProblem<El::Matrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
    if( metadataSummary )
        meta.PrintSummary();
    El::Output("Reading took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( problem.c, "c" );
        El::Print( problem.A, "A" );
        El::Print( problem.b, "b" );
        El::Print( problem.G, "G" );
        El::Print( problem.h, "h" );
    }


    El::AffineLPSolution<El::Matrix<Real>> solution;

    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = progress;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.forceSameStep = forceSameStep;
    ctrl.ipmCtrl.compositeNewton = compositeNewton;
    ctrl.ipmCtrl.maxIts = maxIter;
    ctrl.ipmCtrl.infeasibilityTolLogEps = Real(infeasibilityTolLogEps);
    ctrl.ipmCtrl.relativeObjectiveGapTolLogEps =
      Real(relativeObjectiveGapTolLogEps);
    ctrl.ipmCtrl.relativeComplementarityGapTolLogEps =
      Real(relativeComplementarityGapTolLogEps);
    ctrl.ipmCtrl.xRegLargeLogEps = Real(xRegLargeLogEps);
    ctrl.ipmCtrl.yRegLargeLogEps = Real(yRegLargeLogEps);
    ctrl.ipmCtrl.zRegLargeLogEps = Real(zRegLargeLogEps);
    ctrl.ipmCtrl.xRegSmallLogEps = Real(xRegSmallLogEps);
    ctrl.ipmCtrl.yRegSmallLogEps = Real(yRegSmallLogEps);
    ctrl.ipmCtrl.zRegSmallLogEps = Real(zRegSmallLogEps);
    ctrl.ipmCtrl.lowerTargetRatioLogCompRatio = lowerTargetRatioLogCompRatio;
    ctrl.ipmCtrl.upperTargetRatioLogCompRatio = upperTargetRatioLogCompRatio;

    timer.Start();
    auto info = El::LP( problem, solution, ctrl );
    El::Output("Solving took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( solution.x, "x" );
        El::Print( solution.s, "s" );
        El::Print( solution.y, "y" );
        El::Print( solution.z, "z" );
    }
    El::Output
    ("primal objective:       ",info.ipmInfo.primalObjective);
    El::Output
    ("dual   objective:       ",info.ipmInfo.dualObjective);
    El::Output
    ("infeasibility:          ",info.ipmInfo.infeasibilityError);
    El::Output
    ("relative objective gap: ",info.ipmInfo.relativeObjectiveGap);
    El::Output
    ("relative comp. gap:     ",info.ipmInfo.relativeComplementarityGap);
    El::Output
    ("num iterations:         ",info.ipmInfo.numIterations);
    return info;
}

template<typename Real>
El::LPInfo<Real> SparseLoadAndSolve
( const std::string& filename,
  bool metadataSummary,
  bool print,
  bool progress,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool forceSameStep,
  bool compositeNewton,
  El::Int maxIter )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;
    const bool compressed = false;
    const bool minimize = true;
    const bool keepNonnegativeWithZeroUpperBound=true;

    timer.Start();
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> problem;
    auto meta = El::ReadMPS
    ( problem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
    if( metadataSummary )
        meta.PrintSummary();
    El::Output("Reading took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( problem.c, "c" );
        El::Print( problem.A, "A" );
        El::Print( problem.b, "b" );
        El::Print( problem.G, "G" );
        El::Print( problem.h, "h" );
    }

    // We will wait to do any presolves until a fast mechanism exists for
    // deleting entries from a sparse matrix.

    El::AffineLPSolution<El::Matrix<Real>> solution;

    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = progress;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.forceSameStep = forceSameStep;
    ctrl.ipmCtrl.compositeNewton = compositeNewton;
    ctrl.ipmCtrl.maxIts = maxIter;
    ctrl.ipmCtrl.infeasibilityTolLogEps = Real(infeasibilityTolLogEps);
    ctrl.ipmCtrl.relativeObjectiveGapTolLogEps =
      Real(relativeObjectiveGapTolLogEps);
    ctrl.ipmCtrl.relativeComplementarityGapTolLogEps =
      Real(relativeComplementarityGapTolLogEps);
    ctrl.ipmCtrl.xRegLargeLogEps = Real(xRegLargeLogEps);
    ctrl.ipmCtrl.yRegLargeLogEps = Real(yRegLargeLogEps);
    ctrl.ipmCtrl.zRegLargeLogEps = Real(zRegLargeLogEps);
    ctrl.ipmCtrl.xRegSmallLogEps = Real(xRegSmallLogEps);
    ctrl.ipmCtrl.yRegSmallLogEps = Real(yRegSmallLogEps);
    ctrl.ipmCtrl.zRegSmallLogEps = Real(zRegSmallLogEps);
    ctrl.ipmCtrl.lowerTargetRatioLogCompRatio = lowerTargetRatioLogCompRatio;
    ctrl.ipmCtrl.upperTargetRatioLogCompRatio = upperTargetRatioLogCompRatio;
    ctrl.ipmCtrl.zMinPivotValueLogEps = Real(2.0);

    //ctrl.ipmCtrl.centralityRule = El::StepLengthCentrality<Real>;
    ctrl.ipmCtrl.centralityRule = El::PredictorCorrectorCentrality<Real>;

    timer.Start();
    El::LPInfo<Real> info;
    info.ipmInfo = lp_hsd::affine::IPM( problem, solution, ctrl.ipmCtrl );
    El::Output("Solving took ",timer.Stop()," seconds");
    if( print )
    {
        El::Print( solution.x, "x" );
        El::Print( solution.s, "s" );
        El::Print( solution.y, "y" );
        El::Print( solution.z, "z" );
    }
    El::Output("primal objective:       ",info.ipmInfo.primalObjective);
    El::Output("dual   objective:       ",info.ipmInfo.dualObjective);
    El::Output("infeasibility (scaled): ",info.ipmInfo.infeasibilityError);
    El::Output("relative objective gap: ",info.ipmInfo.relativeObjectiveGap);
    El::Output("relative comp. gap:     ",
      info.ipmInfo.relativeComplementarityGap);
    El::Output("num iterations:         ",info.ipmInfo.numIterations);

    // b - A x
    El::Matrix<Real> primalEqualityInfeas(problem.b);
    El::Multiply
    ( El::NORMAL, Real(1), problem.A, solution.x,
      Real(-1), primalEqualityInfeas );
    const Real primalEqualityInfeasNrm2 =
      El::FrobeniusNorm( primalEqualityInfeas );
    const Real bNrm2 = El::FrobeniusNorm( problem.b );
    El::Output
    ("|| b - A x ||_2 / (1 + || b ||_2) = ",primalEqualityInfeasNrm2/(1+bNrm2));

    // G x + s - h
    El::Matrix<Real> primalConicInfeas(problem.h);
    El::Multiply
    ( El::NORMAL, Real(1), problem.G, solution.x,
      Real(-1), primalConicInfeas );
    primalConicInfeas += solution.s;
    const Real primalConicInfeasNrm2 = El::FrobeniusNorm( primalConicInfeas );
    const Real hNrm2 = El::FrobeniusNorm( problem.h );
    El::Output
    ("|| G x + s - h ||_2 / (1 + || h ||_2) = ",
     primalConicInfeasNrm2/(1+hNrm2));

    // A^T y + G^T z + c
    El::Matrix<Real> dualEqualityInfeas(problem.c);
    El::Multiply
    ( El::TRANSPOSE, Real(1), problem.A, solution.y,
      Real(1), dualEqualityInfeas );
    El::Multiply
    ( El::TRANSPOSE, Real(1), problem.G, solution.z,
      Real(1), dualEqualityInfeas );
    const Real dualEqualityInfeasNrm2 = El::FrobeniusNorm( dualEqualityInfeas );
    const Real cNrm2 = El::FrobeniusNorm( problem.c );
    El::Output
    ("|| A^T y + G^T z + c ||_2 / (1 + || c ||_2) = ",
     dualEqualityInfeasNrm2/(1+cNrm2));

    return info;
}

template<typename Real>
void ConstructPrimalInfeasibilityTest
( const El::AffineLPProblem<El::SparseMatrix<Real>,
                            El::Matrix<Real>>& origProblem,
        El::AffineLPProblem<El::SparseMatrix<Real>,
                            El::Matrix<Real>>& primalInfeasTest,
  bool trivialObjective=true )
{
    EL_DEBUG_CSE
    const El::Int m = origProblem.A.Height();
    const El::Int n = origProblem.A.Width();
    const El::Int k = origProblem.G.Height();

    // Construct the LP
    //
    // min_{yCert,zCert} cCert^T [yCert; zCert]
    //
    // such that
    //
    // | A^T, G^T | | yCert | = |  0 |, | 0, -I | | yCert | <= 0.
    // | b^T, h^T | | zCert |   | -1 |            | zCert |
    //
    // If 'trivialObjective' is true, then 'cCert' is the all zeros vector.
    // Otherwise, 'cCert' is the all ones vector.
    //
    if( trivialObjective )
        El::Zeros( primalInfeasTest.c, m+k, 1 );
    else
        El::Ones( primalInfeasTest.c, m+k, 1 );
    El::Zeros( primalInfeasTest.A, n+1, m+k );
    primalInfeasTest.A.Reserve
    ( origProblem.A.NumEntries() + origProblem.G.NumEntries() + m + n );
    for( El::Int index=0; index<origProblem.A.NumEntries(); ++index )
    {
        const El::Int rowA = origProblem.A.Row(index);
        const El::Int colA = origProblem.A.Col(index);
        const Real value = origProblem.A.Value(index);
        primalInfeasTest.A.QueueUpdate( colA, rowA, value );
    }
    for( El::Int index=0; index<origProblem.G.NumEntries(); ++index )
    {
        const El::Int rowG = origProblem.G.Row(index);
        const El::Int colG = origProblem.G.Col(index);
        const Real value = origProblem.G.Value(index);
        primalInfeasTest.A.QueueUpdate( colG, rowG+m, value );
    }
    for( El::Int index=0; index<m; ++index )
        primalInfeasTest.A.QueueUpdate( n, index, origProblem.b(index) );
    for( El::Int index=0; index<k; ++index )
        primalInfeasTest.A.QueueUpdate( n, index+m, origProblem.h(index) );
    primalInfeasTest.A.ProcessQueues();
    El::Zeros( primalInfeasTest.b, n+1, 1 );
    primalInfeasTest.b(n) = Real(-1);
    El::Zeros( primalInfeasTest.G, k, m+k );
    primalInfeasTest.G.Reserve( k );
    for( El::Int index=0; index<k; ++index )
        primalInfeasTest.G.QueueUpdate( index, index+m, Real(-1) );
    primalInfeasTest.G.ProcessQueues();
    El::Zeros( primalInfeasTest.h, k, 1 );
}

template<typename Real>
void ConstructDualInfeasibilityTest
( const El::AffineLPProblem<El::SparseMatrix<Real>,
                            El::Matrix<Real>>& origProblem,
        El::AffineLPProblem<El::SparseMatrix<Real>,
                            El::Matrix<Real>>& primalInfeasTest,
  bool trivialObjective=true )
{
    EL_DEBUG_CSE
    const El::Int m = origProblem.A.Height();
    const El::Int n = origProblem.A.Width();
    const El::Int k = origProblem.G.Height();

    // Construct the LP
    //
    // min_{xCert,sCert} cCert^T xCert
    //
    // such that
    //
    // | A   | xCert = |  0 |, G xCert + sCert = 0, sCert >= 0.
    // | c^T |         | -1 |
    //
    // If 'trivialObjective' is true, then 'cCert' is the all zeros vector.
    // Otherwise, 'cCert' is the all ones vector.
    //
    if( trivialObjective )
        El::Zeros( primalInfeasTest.c, n, 1 );
    else
        El::Ones( primalInfeasTest.c, n, 1 );
    El::Zeros( primalInfeasTest.A, m+1, n );
    primalInfeasTest.A.Reserve( origProblem.A.NumEntries() + n );
    for( El::Int index=0; index<origProblem.A.NumEntries(); ++index )
    {
        const El::Int rowA = origProblem.A.Row(index);
        const El::Int colA = origProblem.A.Col(index);
        const Real value = origProblem.A.Value(index);
        primalInfeasTest.A.QueueUpdate( rowA, colA, value );
    }
    for( El::Int index=0; index<n; ++index )
        primalInfeasTest.A.QueueUpdate( m, index, origProblem.c(index) );
    primalInfeasTest.A.ProcessQueues();
    El::Zeros( primalInfeasTest.b, m+1, 1 );
    primalInfeasTest.b(m) = Real(-1);
    primalInfeasTest.G = origProblem.G;
    El::Zeros( primalInfeasTest.h, k, 1 );
}

template<typename Real>
void SparseLoadAndTestPrimalFeasibility
( const std::string& filename,
  bool metadataSummary,
  bool print,
  bool progress,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool forceSameStep,
  bool compositeNewton,
  El::Int maxIter )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;
    const bool compressed = false;
    const bool minimize = true;
    const bool keepNonnegativeWithZeroUpperBound=true;

    timer.Start();
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> origProblem;
    auto meta = El::ReadMPS
    ( origProblem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
    if( metadataSummary )
        meta.PrintSummary();
    El::Output("Reading took ",timer.Stop()," seconds");

    const bool trivialObjective = true;
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> infeasProblem;
    ConstructPrimalInfeasibilityTest
    ( origProblem, infeasProblem, trivialObjective );
    if( print )
    {
        El::Print( infeasProblem.c, "cInfeas" );
        El::Print( infeasProblem.A, "AInfeas" );
        El::Print( infeasProblem.b, "bInfeas" );
        El::Print( infeasProblem.G, "GInfeas" );
        El::Print( infeasProblem.h, "hInfeas" );
    }

    // We will wait to do any presolves until a fast mechanism exists for
    // deleting entries from a sparse matrix.

    El::AffineLPSolution<El::Matrix<Real>> solution;

    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = progress;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.forceSameStep = forceSameStep;
    ctrl.ipmCtrl.compositeNewton = compositeNewton;
    ctrl.ipmCtrl.maxIts = maxIter;
    ctrl.ipmCtrl.infeasibilityTolLogEps = Real(infeasibilityTolLogEps);
    ctrl.ipmCtrl.relativeObjectiveGapTolLogEps =
      Real(relativeObjectiveGapTolLogEps);
    ctrl.ipmCtrl.relativeComplementarityGapTolLogEps =
      Real(relativeComplementarityGapTolLogEps);
    ctrl.ipmCtrl.xRegLargeLogEps = Real(xRegLargeLogEps);
    ctrl.ipmCtrl.yRegLargeLogEps = Real(yRegLargeLogEps);
    ctrl.ipmCtrl.zRegLargeLogEps = Real(zRegLargeLogEps);
    ctrl.ipmCtrl.xRegSmallLogEps = Real(xRegSmallLogEps);
    ctrl.ipmCtrl.yRegSmallLogEps = Real(yRegSmallLogEps);
    ctrl.ipmCtrl.zRegSmallLogEps = Real(zRegSmallLogEps);
    ctrl.ipmCtrl.lowerTargetRatioLogCompRatio = lowerTargetRatioLogCompRatio;
    ctrl.ipmCtrl.upperTargetRatioLogCompRatio = upperTargetRatioLogCompRatio;
    ctrl.ipmCtrl.zMinPivotValueLogEps = Real(2.0);

    timer.Start();
    try
    {
        auto info = El::LP( infeasProblem, solution, ctrl );
        El::Output("Solving took ",timer.Stop()," seconds");
        if( print )
        {
            El::Print( solution.x, "xInfeas" );
            El::Print( solution.s, "sInfeas" );
            El::Print( solution.y, "yInfeas" );
            El::Print( solution.z, "zInfeas" );
        }
    }
    catch( const std::exception& except )
    {
        El::Output("LP did not converge.");
        if( print )
        {
            El::Print( solution.x, "xInfeas" );
            El::Print( solution.s, "sInfeas" );
            El::Print( solution.y, "yInfeas" );
            El::Print( solution.z, "zInfeas" );
        }
    }

    const El::Int m = origProblem.A.Height();
    const El::Int n = origProblem.A.Width();
    const El::Int k = origProblem.G.Height();
    auto yCert = solution.x( El::IR(0,m), El::ALL );
    auto zCert = solution.x( El::IR(m,m+k), El::ALL );
    const Real dualObjective = -El::Dot(origProblem.b,yCert) -
      El::Dot(origProblem.h,zCert);
    El::Output("-b^T y - h^T z = ",dualObjective);

    // Test dual feasibility
    El::Matrix<Real> dualEqualityResid;
    El::Zeros( dualEqualityResid, n, 1 );
    El::Multiply
    ( El::TRANSPOSE, Real(1), origProblem.A, yCert,
      Real(1), dualEqualityResid );
    El::Multiply
    ( El::TRANSPOSE, Real(1), origProblem.G, zCert,
      Real(1), dualEqualityResid );
    const Real dualEqualityResidNrm2 = El::FrobeniusNorm( dualEqualityResid );
    El::Output("|| A^T yCert + G^T zCert ||_2 = ",dualEqualityResidNrm2);
}

template<typename Real>
void SparseLoadAndTestDualFeasibility
( const std::string& filename,
  bool metadataSummary,
  bool print,
  bool progress,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool forceSameStep,
  bool compositeNewton,
  El::Int maxIter )
{
    EL_DEBUG_CSE
    El::Output("Will load into El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;
    const bool compressed = false;
    const bool minimize = true;
    const bool keepNonnegativeWithZeroUpperBound=true;

    timer.Start();
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> origProblem;
    auto meta = El::ReadMPS
    ( origProblem, filename, compressed,
      minimize, keepNonnegativeWithZeroUpperBound );
    if( metadataSummary )
        meta.PrintSummary();
    El::Output("Reading took ",timer.Stop()," seconds");

    const bool trivialObjective = true;
    El::AffineLPProblem<El::SparseMatrix<Real>,El::Matrix<Real>> infeasProblem;
    ConstructDualInfeasibilityTest
    ( origProblem, infeasProblem, trivialObjective );
    if( print )
    {
        El::Print( infeasProblem.c, "cInfeas" );
        El::Print( infeasProblem.A, "AInfeas" );
        El::Print( infeasProblem.b, "bInfeas" );
        El::Print( infeasProblem.G, "GInfeas" );
        El::Print( infeasProblem.h, "hInfeas" );
    }

    // We will wait to do any presolves until a fast mechanism exists for
    // deleting entries from a sparse matrix.

    El::AffineLPSolution<El::Matrix<Real>> solution;

    El::lp::affine::Ctrl<Real> ctrl;
    ctrl.ipmCtrl.print = progress;
    ctrl.ipmCtrl.outerEquil = outerEquil;
    ctrl.ipmCtrl.forceSameStep = forceSameStep;
    ctrl.ipmCtrl.compositeNewton = compositeNewton;
    ctrl.ipmCtrl.maxIts = maxIter;
    ctrl.ipmCtrl.infeasibilityTolLogEps = Real(infeasibilityTolLogEps);
    ctrl.ipmCtrl.relativeObjectiveGapTolLogEps =
      Real(relativeObjectiveGapTolLogEps);
    ctrl.ipmCtrl.relativeComplementarityGapTolLogEps =
      Real(relativeComplementarityGapTolLogEps);
    ctrl.ipmCtrl.xRegLargeLogEps = Real(xRegLargeLogEps);
    ctrl.ipmCtrl.yRegLargeLogEps = Real(yRegLargeLogEps);
    ctrl.ipmCtrl.zRegLargeLogEps = Real(zRegLargeLogEps);
    ctrl.ipmCtrl.xRegSmallLogEps = Real(xRegSmallLogEps);
    ctrl.ipmCtrl.yRegSmallLogEps = Real(yRegSmallLogEps);
    ctrl.ipmCtrl.zRegSmallLogEps = Real(zRegSmallLogEps);
    ctrl.ipmCtrl.lowerTargetRatioLogCompRatio = lowerTargetRatioLogCompRatio;
    ctrl.ipmCtrl.upperTargetRatioLogCompRatio = upperTargetRatioLogCompRatio;
    ctrl.ipmCtrl.zMinPivotValueLogEps = Real(2.0);

    timer.Start();
    try
    {
        auto info = El::LP( infeasProblem, solution, ctrl );
        El::Output("Solving took ",timer.Stop()," seconds");
        if( print )
        {
            El::Print( solution.x, "xInfeas" );
            El::Print( solution.s, "sInfeas" );
            El::Print( solution.y, "yInfeas" );
            El::Print( solution.z, "zInfeas" );
        }
    }
    catch( const std::exception& except )
    {
        El::Output("LP did not converge.");
        if( print )
        {
            El::Print( solution.x, "xInfeas" );
            El::Print( solution.s, "sInfeas" );
            El::Print( solution.y, "yInfeas" );
            El::Print( solution.z, "zInfeas" );
        }
    }

    const El::Int m = origProblem.A.Height();
    const El::Int n = origProblem.A.Width();
    const El::Int k = origProblem.G.Height();
    const Real primalObjective = El::Dot(origProblem.c,solution.x);
    El::Output("c^T x = ",primalObjective);

    // Test auxiliary primal feasibility
    // A x
    El::Matrix<Real> primalEqualityResid;
    El::Zeros( primalEqualityResid, m, 1 );
    El::Multiply
    ( El::NORMAL, Real(1), origProblem.A, solution.x,
      Real(1), primalEqualityResid );
    const Real primalEqualityResidNrm2 =
      El::FrobeniusNorm( primalEqualityResid );
    El::Output("|| A xCert ||_2 = ",primalEqualityResidNrm2);
    // G x + s
    El::Matrix<Real> primalConicResid;
    El::Zeros( primalConicResid, k, 1 );
    El::Multiply
    ( El::NORMAL, Real(1), origProblem.G, solution.x,
      Real(1), primalConicResid );
    primalConicResid += solution.s;
    const Real dualEqualityResidNrm2 = El::FrobeniusNorm( primalConicResid );
    El::Output("|| G xCert + sCert ||_2 = ",dualEqualityResidNrm2);
}

template<typename Real>
void SparseNetlibLPData
( const std::string& directory,
  bool metadataSummary,
  bool print,
  bool progress,
  bool outerEquil,
  double infeasibilityTolLogEps,
  double relativeObjectiveGapTolLogEps,
  double relativeComplementarityGapTolLogEps,
  double xRegSmallLogEps,
  double yRegSmallLogEps,
  double zRegSmallLogEps,
  double xRegLargeLogEps,
  double yRegLargeLogEps,
  double zRegLargeLogEps,
  double lowerTargetRatioLogCompRatio,
  double upperTargetRatioLogCompRatio,
  bool forceSameStep,
  bool compositeNewton,
  El::Int maxIter )
{
    EL_DEBUG_CSE
    El::Output
    ("Will test netlib LP_data suite at ",directory,
     " using El::SparseMatrix<",El::TypeName<Real>(),">");
    El::Timer timer;

    std::vector<std::pair<std::string,Real>> problemObjectives;
    // The netlib-documented value for 25fv47 is +5.5018458883e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "25fv47",
      Real(+5.501845888286744794581232588391644178145257e+03) );
    // The netlib-documented value for 80bau3b is +9.8723216072e+05, which
    // appears to only be accurate in its first four digits.
    problemObjectives.emplace_back( "80bau3b",
      Real(+9.872241924090902575791171173879905331234157e+05) );
    // The netlib-documented value for adlittle is +2.2549496316e+05, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "adlittle",
      Real(+2.254949631623803822810117662149174838132293e+05) );
    // The netlib-documented value for afiro is -4.6475314286e+02, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "afiro",
      Real(-4.647531428571428571428571428571428571428571e+02) );
    // The netlib-documented value for agg is -3.5991767287e+07, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "agg",
      Real(-3.599176728657650671264082431963582753710339e+07) );
    // The netlib-documented value for agg2 is -2.0239252356e+07, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "agg2",
      Real(-2.023925235597710902431766192613304279598203e+07) );
    // The netlib-documented value for agg3 is +1.0312115935e+07, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "agg3",
      Real(+1.031211593508922557906105879621508395548832e+07) );
    // The netlib-documented value for bandm is -1.5862801845e+02, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "bandm",
      Real(-1.586280184501206405217412376873610552983025e+02) );
    // The netlib-documented value for beaconfd is +3.3592485807e+04, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "beaconfd",
      Real(+3.359248580720000000000000000000000000000000e+04) );
    // The netlib-documented value for blend is -3.0812149846e+01, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "blend",
      Real(-3.081214984582822017377435612498398053134658e+01) );
    // The netlib-documented value for bnl1 is +1.9776292856e+03, which appears
    // to only be accurate in its first seven digits.
    problemObjectives.emplace_back( "bnl1",
      Real(+1.977629561522889243956439833182149143045779e+03) );
    // The netlib-documented value for bnl2 is +1.8112365404e+03, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "bnl2",
      Real(+1.8112365403585451704484137e+03) );
    // The netlib-documented value for boeing1 is -3.3521356751e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "boeing1",
      Real(-3.35213567507126621842969732e+02) );
    // The netlib-documented value for boeing2 is -3.1501872802e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "boeing2",
      Real(-3.150187280152028787046220e+02) );
    // The netlib-documented value for bore3d is +1.3730803942e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "bore3d",
      Real(+1.3730803942084927215581987e+03) );
    // The netlib-documented value for brandy is +1.5185098965e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "brandy",
      Real(+1.5185098964881283835426751e+03) );
    // The netlib-documented value for capri is +2.6900129138e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "capri",
      Real(+2.6900129137681610087717280694e+03) );
    // The netlib-documented value for cycle is -5.2263930249e+00, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "cycle",
      Real(-5.2263930248941017172447e+00) );
    // The netlib-documented value for czprob is +2.1851966989e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "czprob",
      Real(+2.185196698856577485895115595e+06) );
    // The netlib-documented value for d2q06c is +1.2278423615e+05, which
    // appears to only be accurate in its first seven digits.
    problemObjectives.emplace_back( "d2q06c",
      Real(+1.2278421081418945895739128812e+05) );
    // The netlib-documented value for d6cube is +3.1549166667e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "d6cube",
      Real(+3.1549166666666666666666666667e+02) );
    // The netlib-documented value for degen2 is -1.4351780000e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "degen2",
      Real(-1.435178000000000000000000000e+03) );
    // The netlib-documented value for degen3 is -9.8729400000e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "degen3",
      Real(-9.872940000000000000000000000e+02) );
    // The netlib-documented value for dfl001 is +1.126639607e+07, which appears
    // to only be accurate in its first ten digits.
    problemObjectives.emplace_back( "dfl001",
      Real(+1.126639604667139220238e+07) );
    // The netlib-documented value for e226 is -1.8751929066e+01, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "e226",
      Real(-1.875192906637054910260569e+01) );
    // The netlib-documented value for etamacro is -7.5571521774e+02, which
    // appears to only be accurate in its first seven digits.
    problemObjectives.emplace_back( "etamacro",
      Real(-7.5571523337491333507925837e+02) );
    // The netlib-documented value for fffff800 is +5.5567961165e+05, which
    // appears to only be accurate in its first seven digits.
    problemObjectives.emplace_back( "fffff800",
      Real(+5.556795648174963765328643789e+05) );
    // The netlib-documented value for finnis is +1.7279096547e+05, which
    // appears to only be accurate in its first ten digits.
    problemObjectives.emplace_back( "finnis",
      Real(+1.7279106559561159432297900375e+05) );
    // The netlib-documented value for fit1d is -9.1463780924e+03, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "fit1d",
      Real(-9.1463780924209269467749025e+03) );
    // The netlib-documented value for fit1p is +9.1463780924e+03, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "fit1p",
      Real(+9.1463780924209269467749025e+03) );
    // The netlib-documented value for fit2d is -6.8464293294e+04, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "fit2d",
      Real(-6.846429329383206957594351843583743643335e+04) );
    // The netlib-documented value for fit2p is +6.8464293232e+04, which
    // appears to only be accurate in its first nine digits.
    problemObjectives.emplace_back( "fit2p",
      Real(+6.84642932938320695759435184358374e+04) );
    // The netlib-documented value for forplan is -6.6421873953e+02, which
    // appears to only be accurate in its first six digits.
    problemObjectives.emplace_back( "forplan",
      Real(-6.642189612722045748123511970169171790752754902660249031235e+02) );
    // The netlib-documented value for ganges is -1.0958636356e+05, which
    // appears to only be accurate in its first six digits.
    problemObjectives.emplace_back( "ganges",
      Real(-1.09585736129278116234449e+05) );
    // The netlib-documented value for gfrd-pnc is +6.9022359995e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "gfrd-pnc",
      Real(+6.9022359995488088296e+06) );
    // The netlib-documented value for greenbea is -7.2462405908e+07, which
    // appears to only be accurate in its first two (!!!) digits.
    problemObjectives.emplace_back( "greenbea",
      Real(-7.25552481298459874575578706e+07) );
    // The netlib-documented value for greenbeb is -4.3021476065e+06, which
    // appears to only be accurate in its first four digits.
    problemObjectives.emplace_back( "greenbeb",
      Real(-4.30226026120658675392137e+06) );
    // The netlib-documented value for grow15 is -1.0687094129e+08, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "grow15",
      Real(-1.068709412935753367160404093031e+08) );
    // The netlib-documented value for grow22 is -1.6083433648e+08, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "grow22",
      Real(-1.608343364825629671845603998261e+08) );
    // The netlib-documented value for grow7 is -4.7787811815e+07, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "grow7",
      Real(-4.778781181471150261677e+07) );
    // The netlib-documented value for israel is -8.9664482186e+05, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "israel",
      Real(-8.96644821863045729662e+05) );
    // The netlib-documented value for kb2 is -1.7499001299e+03, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "kb2",
      Real(-1.749900129906206e+03) );
    // The netlib-documented value for lotfi is -2.5264706062e+01, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "lotfi",
      Real(-2.52647060618800000000e+01) );
    // The netlib-documented value for maros is -5.8063743701e+04, which appears
    // to be accurate to all eleven digits.
    problemObjectives.emplace_back( "maros",
      Real(-5.80637437011258954012085e+04) );
    // The netlib-documented value for maros-r7 is +1.4971851665e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "maros-r7",
      Real(+1.4971851664796437907337544e+06) );
    // The netlib-documented value for modszk1 is +3.2061972906e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "modszk1",
      Real(+3.20619729064315804943e+02) );
    // The netlib-documented value for nesm is +1.4076073035e+07, which appears
    // to only be accurate in its first six digits.
    problemObjectives.emplace_back( "nesm",
      Real(+1.407603648756272834e+07) );
    // The netlib-documented value for perold is -9.3807580773e+03, which
    // appears to only be accurate in its first six digits.
    problemObjectives.emplace_back( "perold",
      Real(-9.380755278235160673483327038651070221818e+03) );
    // The netlib-documented value for pilot is -5.5740430007e+02, which appears
    // to only be accurate in its first four digits.
    problemObjectives.emplace_back( "pilot",
      Real(-5.574897292840681807303e+02) );
    // The netlib-documented value for pilot.ja is -6.1131344111e+03, which
    // appears to only be accurate in its first six digits.
    problemObjectives.emplace_back( "pilot.ja",
      Real(-6.1131364655813432748849e+03) );
    // The netlib-documented value for pilot.we is -2.7201027439e+06, which
    // appears to only be accurate in its first six digits.
    problemObjectives.emplace_back( "pilot.we",
      Real(-2.72010753284496396294e+06) );
    // The netlib-documented value for pilot4 is -2.5811392641e+03, which
    // appears to only be accurate in its first eight digits.
    problemObjectives.emplace_back( "pilot4",
      Real(-2.581139258883888674583099727e+03) );
    // The netlib-documented value for pilot87 is +3.0171072827e+02, which
    // appears to only be accurate in its first six digits.
    problemObjectives.emplace_back( "pilot87",
      Real(+3.01710347333110527721660e+02) );
    // The netlib-documented value for pilotnov is -4.4972761882e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "pilotnov",
      Real(-4.49727618821887114309962117839e+03) );
    /*
    problemObjectives.emplace_back( "qap8",  Real(+2.0350000000e+02) );
    problemObjectives.emplace_back( "qap12", Real(+5.2289435056e+02) );
    problemObjectives.emplace_back( "qap15", Real(+1.0409940410e+03) );
    */
    // The netlib-documented value for recipe is -2.6661600000e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "recipe",
      Real(-2.66616000000000000000000000000e+02) );
    // The netlib-documented value for sc105 is -5.2202061212e+01, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "sc105",
      Real(-5.22020612117072480626280e+01) );
    // The netlib-documented value for sc205 is -5.2202061212e+01, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "sc205",
      Real(-5.220206121170724806262801086e+01) );
    // The netlib-documented value for sc50a is -6.4575077059e+01, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "sc50a",
      Real(-6.4575077058564509026860413915e+01) );
    // The netlib-documented value for sc50b is -7.0000000000e+01, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "sc50b",
      Real(-7.000000000000000e+01) );
    // The netlib-documented value for scagr25 is -1.4753433061e+07, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "scagr25",
      Real(-1.475343306076852316779092508e+07) );
    // The netlib-documented value for scagr7 is -2.3313892548e+06, which
    // appears to only be accurate in its first seven digits.
    problemObjectives.emplace_back( "scagr7",
      Real(-2.331389824330984000e+06) );
    // The netlib-documented value for scfxm1 is +1.8416759028e+04, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "scfxm1",
      Real(+1.84167590283489436836e+04) );
    // The netlib-documented value for scfxm2 is +3.6660261565e+04, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "scfxm2",
      Real(+3.6660261564998812956939505e+04) );
    // The netlib-documented value for scfxm3 is +5.4901254550e+04, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "scfxm3",
      Real(+5.49012545497514546238887249e+04) );
    // The netlib-documented value for scorpion is +1.8781248227e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "scorpion",
      Real(+1.87812482273810662964794117636e+03) );
    // The netlib-documented value for scrs8 is +9.0429998619e+02, which
    // appears to only be accurate in its first seven digits.
    problemObjectives.emplace_back( "scrs8",
      Real(+9.0429695380079143579923108e+02) );
    // The netlib-documented value for scsd1 is +8.6666666743e+00, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "scsd1",
      Real(+8.666666674333364729253e+00) );
    // The netlib-documented value for scsd6 is +5.0500000078e+01, which
    // appears to only be accurate in its first ten digits.
    problemObjectives.emplace_back( "scsd6",
      Real(+5.050000007714434528998549e+01) );
    // The netlib-documented value for scsd8 is +9.0499999993e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "scsd8",
      Real(+9.049999999254644009e+02) );
    // The netlib-documented value for sctap1 is +1.4122500000e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "sctap1",
      Real(+1.41225000000000000000000000000000e+03) );
    // The netlib-documented value for sctap2 is +1.7248071429e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "sctap2",
      Real(+1.724807142857142857143e+03) );
    // The netlib-documented value for sctap3 is +1.4240000000e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "sctap3",
      Real(+1.42400000000000000000e+03) );
    // The netlib-documented value for seba is +1.5711600000e+04, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "seba",
      Real(+1.57116000000000000000000000000000000000000000000e+04) );
    // The netlib-documented value for share1b is -7.6589318579e+04, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "share1b",
      Real(-7.65893185791856811280e+04) );
    // The netlib-documented value for share2b is -4.1573224074e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "share2b",
      Real(-4.1573224074141948654519910873905e+02) );
    // The netlib-documented value for shell is +1.2088253460e+09, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "shell",
      Real(+1.2088253460000000000000000000000e+09) );
    // The netlib-documented value for ship04l is +1.7933245380e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "ship04l",
      Real(+1.79332453797035576255625563e+06) );
    // The netlib-documented value for ship04s is +1.7987147004e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "ship04s",
      Real(+1.798714700445392171995408272605e+06) );
    // The netlib-documented value for ship08l is +1.9090552114e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "ship08l",
      Real(+1.9090552113891315179803820e+06) );
    // The netlib-documented value for ship08s is +1.9200982105e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "ship08s",
      Real(+1.9200982105346197101726954747e+06) );
    // The netlib-documented value for ship12l is +1.4701879193e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "ship12l",
      Real(+1.4701879193292648227021755439e+06) );
    // The netlib-documented value for ship12s is +1.4892361344e+06, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "ship12s",
      Real(+1.4892361344061291565686421605401e+06) );
    // The netlib-documented value for sierra is +1.5394362184e+07, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "sierra",
      Real(+1.5394362183631930000e+07) );
    // The netlib-documented value for stair is -2.5126695119e+02, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "stair",
      Real(-2.5126695119296330352803637106304e+02) );
    // The netlib-documented value for standata is +1.2576995000e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "standata",
      Real(+1.25769950000000000000000000e+03) );
    /* Skipping 'standgub' */
    // The netlib-documented value for standmps is +1.4060175000e+03, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "standmps",
      Real(+1.406017500000000000000000e+03) );
    // The netlib-documented value for stocfor1 is -4.1131976219e+04, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "stocfor1",
      Real(-4.11319762194364060656827607315e+04) );
    // The netlib-documented value for stocfor2 is -3.9024408538e+04, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "stocfor2",
      Real(-3.902440853788202960e+04) );
    /*
    problemObjectives.emplace_back( "stocfor3", Real(-3.9976661576e+04) );
    problemObjectives.emplace_back( "truss",    Real(+4.5881584719e+05) );
    */
    // The netlib-documented value for tuff is +2.9214776509e-01, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "tuff",
      Real(+2.9214776509361284488226e-01) );
    // The netlib-documented value for vtp.base is +1.2983146246e+05, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "vtp.base",
      Real(+1.298314624613613657396e+05) );
    // The netlib-documented value for wood1p is +1.4429024116e+00, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "wood1p",
      Real(+1.44290241157340924000109e+00) );
    // The netlib-documented value for woodw is +1.3044763331e+00, which
    // appears to be accurate to all eleven digits.
    problemObjectives.emplace_back( "woodw",
      Real(+1.304476333084229269005552085e+00) );

    const Real epsilon = El::limits::Epsilon<Real>();
    const Real demandedRelativeError =
      El::Pow( epsilon, relativeObjectiveGapTolLogEps );
    std::vector<std::pair<Real,std::string>> problemRelErrors;
    std::vector<std::pair<El::Int,std::string>> problemIterations;
    for( const auto& problemObjective : problemObjectives )
    {
        const std::string& name = problemObjective.first;
        const Real& objective = problemObjective.second;
        El::Output("Testing ",name);
        const std::string filename = directory + "/" + name + ".mps";
        auto info =
          SparseLoadAndSolve<Real>
          ( filename, metadataSummary, print, progress, outerEquil,
            infeasibilityTolLogEps,
            relativeObjectiveGapTolLogEps,
            relativeComplementarityGapTolLogEps,
            xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
            xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
            lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
            forceSameStep, compositeNewton, maxIter );
        const Real relativeError =
          El::Abs(objective - info.ipmInfo.primalObjective) /
          El::Abs(objective);
        if( relativeError > demandedRelativeError )
        {
            El::Output
            ("WARNING: Only solved ",name," to objective ",
             info.ipmInfo.primalObjective," (vs. ",objective,")");
        }
       problemRelErrors.emplace_back( relativeError, name );
       problemIterations.emplace_back( info.ipmInfo.numIterations, name );
    }
    El::Output("");

    const El::Int numWorstErrors = 10;
    std::sort( problemRelErrors.rbegin(), problemRelErrors.rend() );
    El::Output(numWorstErrors," problems with highest error:");
    for( El::Int i=0; i<El::Min(numWorstErrors,problemRelErrors.size()); ++i )
    {
        const auto& problemRelError = problemRelErrors[i];
        El::Output
        (problemRelError.second,
         " was only solved to relative accuracy of ",problemRelError.first);
    }
    El::Output("");

    const El::Int numMostIters = 10;
    std::sort( problemIterations.rbegin(), problemIterations.rend() );
    El::Output(numMostIters," problems taking the most iterations:");
    for( El::Int i=0; i<El::Min(numMostIters,problemIterations.size()); ++i )
    {
        const auto& problemIteration = problemIterations[i];
        El::Output
        (problemIteration.second," took ",problemIteration.first," iterations");
    }
}

int main( int argc, char* argv[] )
{
    El::Environment env( argc, argv );

    try
    {
        const std::string filename =
          El::Input
          ("--filename","MPS filename",
           std::string("../data/optimization/lp_data/share1b.mps"));
        const bool metadataSummary =
          El::Input("--metadataSummary","summarize MPS metadata?",true);
        const bool testDense =
          El::Input("--testDense","test with dense matrices?",false);
        const bool testDouble =
          El::Input("--testDouble","test double-precision?",true);
        const bool testArbitrary =
          El::Input("--testArbitrary","test arbitrary precision?",false);
        const bool print = El::Input("--print","print matrices?",false);
        const bool progress = El::Input("--progress","IPM progress?",true);
        const bool outerEquil =
          El::Input("--outerEquil","outer equilibration?",true);
        const double infeasibilityTolLogEps =
          El::Input
          ("--infeasibilityTolLogEps","log_eps(infeasibilityTol)",0.45);
        const double relativeObjectiveGapTolLogEps =
          El::Input
          ("--relativeObjectiveGapTolLogEps",
           "log_eps(relativeObjectiveGapTol)",0.05);
        const double relativeComplementarityGapTolLogEps =
          El::Input
          ("--relativeComplementarityGapTolLogEps",
           "log_eps(relativeComplementarityGapTol)",0.3);
        const double xRegSmallLogEps =
          El::Input("--xRegSmallLogEps","log_eps(xRegSmall)",0.8);
        const double yRegSmallLogEps =
          El::Input("--yRegSmallLogEps","log_eps(yRegSmall)",0.8);
        const double zRegSmallLogEps =
          El::Input("--zRegSmallLogEps","log_eps(zRegSmall)",0.8);
        const double xRegLargeLogEps =
          El::Input("--xRegLargeLogEps","log_eps(xRegLarge)",0.6);
        const double yRegLargeLogEps =
          El::Input("--yRegLargeLogEps","log_eps(yRegLarge)",0.6);
        const double zRegLargeLogEps =
          El::Input("--zRegLargeLogEps","log_eps(zRegLarge)",0.6);
        const double lowerTargetRatioLogCompRatio =
          El::Input
          ("--lowerTargetRatioLogCompRatio","log_compratio(lowerTargetRatio)",
           -0.25);
        const double upperTargetRatioLogCompRatio =
          El::Input
          ("--upperTargetRatioLogCompRatio","log_compratio(upperTargetRatio)",
           0.25);
        const bool forceSameStep =
          El::Input
          ("--forceSameStep","force same primal and dual step sizes?",false);
        const bool compositeNewton =
          El::Input("--compositeNewton","Mehrotra predictor-corrector?",false);
        const El::Int maxIter =
          El::Input("--maxIter","max. IPM iterations",1000);
        const bool testNetlib =
          El::Input("--testNetlib","test entire netlib LP_data?",false);
        const std::string netlibDirectory =
          El::Input
          ("--netlibDirectory","path to netlib LP_data MPS files",
           std::string("../data/optimization/lp_data"));
#ifdef EL_HAVE_MPC
        const mpfr_prec_t prec = El::Input("--mpfr_prec","MPFR precision",512);
#endif
        El::ProcessInput();
        El::PrintInputReport();

        if( testNetlib )
        {
            // TODO(poulson): Provide a command-line interface for testing the
            // netlib suite with other datatypes.
            SparseNetlibLPData<double>
            ( netlibDirectory, metadataSummary, print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton, maxIter );
            return 0;
        }

        if( testDense )
        {
            if( testDouble )
                DenseLoadAndSolve<double>
                ( filename, metadataSummary,
                  print, progress, outerEquil,
                  infeasibilityTolLogEps,
                  relativeObjectiveGapTolLogEps,
                  relativeComplementarityGapTolLogEps,
                  xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
                  xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
                  lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
                  forceSameStep, compositeNewton, maxIter );
#ifdef EL_HAVE_QD
            DenseLoadAndSolve<El::DoubleDouble>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton, maxIter );
            DenseLoadAndSolve<El::QuadDouble>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton, maxIter );
#endif
#ifdef EL_HAVE_QUAD
            DenseLoadAndSolve<El::Quad>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton, maxIter );
#endif
#ifdef EL_HAVE_MPC
            if( testArbitrary )
            {
                // TODO(poulson): Make this configurable.
                El::mpfr::SetPrecision( prec );
                DenseLoadAndSolve<El::BigFloat>
                ( filename, metadataSummary,
                  print, progress, outerEquil,
                  infeasibilityTolLogEps,
                  relativeObjectiveGapTolLogEps,
                  relativeComplementarityGapTolLogEps,
                  xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
                  xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
                  lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
                  forceSameStep, compositeNewton, maxIter );
            }
#endif
        }

        if( testDouble )
            SparseLoadAndSolve<double>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton, maxIter );
#ifdef EL_HAVE_QD
        //SparseLoadAndTestPrimalFeasibility<El::DoubleDouble>
        SparseLoadAndSolve<El::DoubleDouble>
        ( filename, metadataSummary,
          print, progress, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          forceSameStep, compositeNewton, maxIter );
        //SparseLoadAndTestPrimalFeasibility<El::QuadDouble>
        SparseLoadAndSolve<El::QuadDouble>
        ( filename, metadataSummary,
          print, progress, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          forceSameStep, compositeNewton, maxIter );
#endif
#ifdef EL_HAVE_QUAD
        //SparseLoadAndTestPrimalFeasibility<El::Quad>
        SparseLoadAndSolve<El::Quad>
        ( filename, metadataSummary,
          print, progress, outerEquil,
          infeasibilityTolLogEps,
          relativeObjectiveGapTolLogEps,
          relativeComplementarityGapTolLogEps,
          xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
          xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
          lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
          forceSameStep, compositeNewton, maxIter );
#endif
#ifdef EL_HAVE_MPC
        if( testArbitrary )
        {
            // TODO(poulson): Make this configurable
            El::mpfr::SetPrecision( prec );
            //SparseLoadAndTestPrimalFeasibility<El::BigFloat>
            SparseLoadAndSolve<El::BigFloat>
            ( filename, metadataSummary,
              print, progress, outerEquil,
              infeasibilityTolLogEps,
              relativeObjectiveGapTolLogEps,
              relativeComplementarityGapTolLogEps,
              xRegSmallLogEps, yRegSmallLogEps, zRegSmallLogEps,
              xRegLargeLogEps, yRegLargeLogEps, zRegLargeLogEps,
              lowerTargetRatioLogCompRatio, upperTargetRatioLogCompRatio,
              forceSameStep, compositeNewton, maxIter );
        }
#endif
    }
    catch( std::exception& e ) { El::ReportException(e); }

    return 0;
}
