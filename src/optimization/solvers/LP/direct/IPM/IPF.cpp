/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "./util.hpp"

namespace El {
namespace lp {
namespace direct {

// The following solves a pair of linear programs in "direct" conic form:
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
// using a simple Infeasible Path Following (IPF) scheme. This routine
// should only be used for academic purposes, as the Mehrotra alternative
// typically requires an order of magnitude fewer iterations.

template<typename Real>
void IPF
( const Matrix<Real>& APre, 
  const Matrix<Real>& bPre, const Matrix<Real>& cPre,
        Matrix<Real>& x,          Matrix<Real>& y, 
        Matrix<Real>& z,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::IPF"))    

    // Equilibrate the LP by diagonally scaling A
    auto A = APre;
    Matrix<Real> dRow, dCol;
    GeomEquil( A, dRow, dCol );
    const Int m = A.Height();
    const Int n = A.Width();
    auto b = bPre;
    auto c = cPre;
    DiagonalSolve( LEFT, NORMAL, dRow, b );
    DiagonalSolve( LEFT, NORMAL, dCol, c );

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );

    // TODO: Expose this as a parameter to IPFCtrl
    const bool standardShift = true;
    Initialize
    ( A, b, c, x, y, z, ctrl.primalInitialized, ctrl.dualInitialized,
      standardShift );

    Matrix<Real> J, d, 
                 rc, rb, rmu, 
                 dx, dy, dz;
#ifndef EL_RELEASE
    Matrix<Real> dxError, dyError, dzError, prod;
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = NumNonPositive( x );
        const Int zNumNonPos = NumNonPositive( z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y); 
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        Scale( Real(-1), rb );
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Axpy( Real(-1), z, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // Now check the pieces
        // --------------------
        if( ctrl.print )
            cout << " iter " << numIts << ":\n"
                 << "  |primal - dual| / (1 + |primal|) = "
                 << objConv << "\n"
                 << "  || r_b ||_2 / (1 + || b ||_2)   = "
                 << rbConv << "\n"
                 << "  || r_c ||_2 / (1 + || c ||_2)   = "
                 << rcConv << endl;
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Compute the duality measure and r_mu = x o z - tau e
        // ====================================================
        const Real mu = Dot(x,z) / n;
        rmu = z;
        DiagonalScale( LEFT, NORMAL, x, rmu );
        Shift( rmu, -ctrl.centering*mu );

        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // =============================
            KKT( A, x, z, J );
            KKTRHS( rc, rb, rmu, z, d );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, d );
            ExpandSolution( m, n, d, dx, dy, dz );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the reduced KKT system
            // ================================
            AugmentedKKT( A, x, z, J );
            AugmentedKKTRHS( x, rc, rb, rmu, d );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, d );
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }
        else if( ctrl.system == NORMAL_KKT )
        {
            // Construct the reduced KKT system
            // ================================
            NormalKKT( A, x, z, J );
            NormalKKTRHS( A, x, z, rc, rb, rmu, dy );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, dy );
            ExpandNormalSolution( A, c, x, z, rc, rmu, dx, dy, dz );
        }
        else
            LogicError("Invalid KKT system choice");
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Axpy( Real(-1), dz, dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        const Real rmuNrm2 = Nrm2( rmu );
        dzError = rmu;
        prod = dz;
        DiagonalScale( LEFT, NORMAL, x, prod );
        Axpy( Real(1), prod, dzError );
        prod = dx;
        DiagonalScale( LEFT, NORMAL, z, prod ); 
        Axpy( Real(1), prod, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        if( ctrl.print )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = " 
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = " 
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_mu ||_2) = " 
                 << dzErrorNrm2/(1+rmuNrm2) << endl;
#endif

        // Take a step in the computed direction
        // =====================================
        const Real alphaPrimal = MaxStepInPositiveCone( x, dx, Real(1) );
        const Real alphaDual = MaxStepInPositiveCone( z, dz, Real(1) );
        const Real alphaMax = Min(alphaPrimal,alphaDual);
        if( ctrl.print )
            cout << "alphaMax = " << alphaMax << endl;
        const Real alpha =
          IPFLineSearch
          ( A, b, c, x, y, z, dx, dy, dz, 
            Real(0.99)*alphaMax,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.lineSearchCtrl );
        if( ctrl.print )
            cout << "  alpha = " << alpha << endl;
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
    }

    // Unequilibrate the LP
    DiagonalSolve( LEFT, NORMAL, dCol, x );
    DiagonalSolve( LEFT, NORMAL, dRow, y );
    DiagonalScale( LEFT, NORMAL, dCol, z );
}

template<typename Real>
void IPF
( const AbstractDistMatrix<Real>& APre, 
  const AbstractDistMatrix<Real>& bPre, const AbstractDistMatrix<Real>& cPre,
        AbstractDistMatrix<Real>& xPre,       AbstractDistMatrix<Real>& yPre, 
        AbstractDistMatrix<Real>& zPre, 
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::IPF"))    
    const Grid& grid = APre.Grid();
    const int commRank = grid.Rank();

    // Ensure that the inputs have the appropriate read/write properties
    DistMatrix<Real> A(grid), b(grid), c(grid);
    A.Align(0,0);
    b.Align(0,0);
    c.Align(0,0);
    A = APre;
    b = bPre;
    c = cPre;
    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    // NOTE: x does not need to be a read proxy when !ctrl.primalInitialized
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,control); auto& x = *xPtr;
    // NOTE: {y,z} do not need to be read proxies when !ctrl.dualInitialized
    auto yPtr = ReadWriteProxy<Real,MC,MR>(&yPre,control); auto& y = *yPtr;
    auto zPtr = ReadWriteProxy<Real,MC,MR>(&zPre,control); auto& z = *zPtr;

    // Equilibrate the LP by diagonally scaling A
    DistMatrix<Real,MC,STAR> dRow(grid);
    DistMatrix<Real,MR,STAR> dCol(grid);
    GeomEquil( A, dRow, dCol );
    const Int m = A.Height();
    const Int n = A.Width();
    DiagonalSolve( LEFT, NORMAL, dRow, b );
    DiagonalSolve( LEFT, NORMAL, dCol, c );

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );

    // TODO: Expose this as a parameter to IPFCtrl
    const bool standardShift = true;
    Initialize
    ( A, b, c, x, y, z, ctrl.primalInitialized, ctrl.dualInitialized,
      standardShift );

    DistMatrix<Real> J(grid), d(grid), 
                     rc(grid), rb(grid), rmu(grid),
                     dx(grid), dy(grid), dz(grid);
    dx.AlignWith( x );
    dz.AlignWith( x );
    rmu.AlignWith( x );
#ifndef EL_RELEASE
    DistMatrix<Real> dxError(grid), dyError(grid), dzError(grid), prod(grid);
    dzError.AlignWith( dz );
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = NumNonPositive( x );
        const Int zNumNonPos = NumNonPositive( z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        Scale( Real(-1), rb );
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Axpy( Real(-1), z, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // Now check the pieces
        // --------------------
        if( ctrl.print && commRank == 0 )
            cout << " iter " << numIts << ":\n"
                 << "  |primal - dual| / (1 + |primal|) = "
                 << objConv << "\n"
                 << "  || r_b ||_2 / (1 + || b ||_2)   = "
                 << rbConv << "\n"
                 << "  || r_c ||_2 / (1 + || c ||_2)   = "
                 << rcConv << endl;
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Compute the duality measure and r_mu = x o z - tau e
        // ====================================================
        const Real mu = Dot(x,z) / n;
        rmu = z;
        DiagonalScale( LEFT, NORMAL, x, rmu );
        Shift( rmu, -ctrl.centering*mu );

        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // =============================
            KKT( A, x, z, J );
            KKTRHS( rc, rb, rmu, z, d );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, d );
            ExpandSolution( m, n, d, dx, dy, dz );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the reduced KKT system
            // ================================
            AugmentedKKT( A, x, z, J );
            AugmentedKKTRHS( x, rc, rb, rmu, d );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, d );
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }
        else if( ctrl.system == NORMAL_KKT )
        {
            // Construct the reduced KKT system
            // ================================
            NormalKKT( A, x, z, J );
            NormalKKTRHS( A, x, z, rc, rb, rmu, dy );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, dy );
            ExpandNormalSolution( A, c, x, z, rc, rmu, dx, dy, dz );
        }
        else
            LogicError("Invalid KKT system choice");
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Axpy( Real(-1), dz, dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        const Real rmuNrm2 = Nrm2( rmu );
        dzError = rmu;
        prod = dz;
        DiagonalScale( LEFT, NORMAL, x, prod );
        Axpy( Real(1), prod, dzError );
        prod = dx;
        DiagonalScale( LEFT, NORMAL, z, prod );
        Axpy( Real(1), prod, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        if( ctrl.print && commRank == 0 )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = " 
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = " 
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_mu ||_2) = " 
                 << dzErrorNrm2/(1+rmuNrm2) << endl;
#endif

        // Take a step in the computed direction
        // =====================================
        const Real alphaPrimal = MaxStepInPositiveCone( x, dx, Real(1) );
        const Real alphaDual = MaxStepInPositiveCone( z, dz, Real(1) );
        const Real alphaMax = Min(alphaPrimal,alphaDual);
        if( ctrl.print && commRank == 0 )
            cout << "alphaMax = " << alphaMax << endl;
        const Real alpha =
          IPFLineSearch
          ( A, b, c, x, y, z, dx, dy, dz, 
            Real(0.99)*alphaMax,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.lineSearchCtrl );
        if( ctrl.print && commRank == 0 )
            cout << "  alpha = " << alpha << endl;
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
    }

    // Unequilibrate the LP
    DiagonalSolve( LEFT, NORMAL, dCol, x );
    DiagonalSolve( LEFT, NORMAL, dRow, y );
    DiagonalScale( LEFT, NORMAL, dCol, z );
}

template<typename Real>
void IPF
( const SparseMatrix<Real>& APre, 
  const Matrix<Real>& bPre,       const Matrix<Real>& cPre,
        Matrix<Real>& x,                Matrix<Real>& y, 
        Matrix<Real>& z,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::IPF"))    
    const Real epsilon = lapack::MachineEpsilon<Real>();

    // Equilibrate the LP by diagonally scaling A
    auto A = APre;
    Matrix<Real> dRow, dCol;
    GeomEquil( A, dRow, dCol );
    const Int m = A.Height();
    const Int n = A.Width();
    auto b = bPre;
    auto c = cPre;
    DiagonalSolve( LEFT, NORMAL, dRow, b );
    DiagonalSolve( LEFT, NORMAL, dCol, c );

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );

    vector<Int> map, invMap;
    SymmNodeInfo info;
    Separator rootSep;
    // The initialization involves an augmented KKT system, and so we can
    // only reuse the factorization metadata if the this IPM is using the
    // augmented formulation
    // TODO: Expose this as a parameter to IPFCtrl
    const bool standardShift = true;
    if( ctrl.system == AUGMENTED_KKT )
    {
        Initialize
        ( A, b, c, x, y, z, map, invMap, rootSep, info,
          ctrl.primalInitialized, ctrl.dualInitialized, standardShift,
          ctrl.print );
    }
    else
    {
        vector<Int> augMap, augInvMap;
        SymmNodeInfo augInfo;
        Separator augRootSep;
        Initialize
        ( A, b, c, x, y, z, augMap, augInvMap, augRootSep, augInfo, 
          ctrl.primalInitialized, ctrl.dualInitialized, standardShift,
          ctrl.print );
    }

    SparseMatrix<Real> J;
    SymmFront<Real> JFront;
    Matrix<Real> d,
                 rc, rb, rmu, 
                 dx, dy, dz;

    Matrix<Real> regCand, reg;
    // TODO: Dynamically modify these values in the manner suggested by 
    //       Altman and Gondzio based upon the number of performed steps of
    //       iterative refinement
    if( ctrl.system == FULL_KKT )
    {
        const Real regMagPrimal = Pow(epsilon,Real(0.75));
        const Real regMagLagrange = Pow(epsilon,Real(0.5));
        const Real regMagDual = Pow(epsilon,Real(0.5));
        regCand.Resize( m+2*n, 1 );
        for( Int i=0; i<m+2*n; ++i )
        {
            if( i < n )
                regCand.Set( i, 0, regMagPrimal );
            else if( i < n+m )
                regCand.Set( i, 0, -regMagLagrange );
            else
                regCand.Set( i, 0, -regMagDual );
        }
    }
    else if( ctrl.system == AUGMENTED_KKT )
    {
        const Real regMagPrimal = Pow(epsilon,Real(0.75));
        const Real regMagLagrange = Pow(epsilon,Real(0.5));
        regCand.Resize( n+m, 1 );
        for( Int i=0; i<n+m; ++i )
        {
            if( i < n )
                regCand.Set( i, 0, regMagPrimal );
            else
                regCand.Set( i, 0, -regMagLagrange );
        }
    }
    MatrixNode<Real> regCandNodal, regNodal;
    bool increasedReg = false;

#ifndef EL_RELEASE
    Matrix<Real> dxError, dyError, dzError, prod;
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = NumNonPositive( x );
        const Int zNumNonPos = NumNonPositive( z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        Scale( Real(-1), rb );
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Axpy( Real(-1), z, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // Now check the pieces
        // --------------------
        if( ctrl.print )
            cout << " iter " << numIts << ":\n"
                 << "  |primal - dual| / (1 + |primal|) = "
                 << objConv << "\n"
                 << "  || r_b ||_2 / (1 + || b ||_2)   = "
                 << rbConv << "\n"
                 << "  || r_c ||_2 / (1 + || c ||_2)   = "
                 << rcConv << endl;
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Compute the duality measure and r_mu = x o z - tau e
        // ====================================================
        const Real mu = Dot(x,z) / n;
        rmu = z;
        DiagonalScale( LEFT, NORMAL, x, rmu );
        Shift( rmu, -ctrl.centering*mu );

        // Compute the search direction
        // ============================
        // TODO: Expose these as control parameters
        const Real minReductionFactor = 2;
        const Int maxRefineIts = 50;
        bool aPriori = true;
        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // -----------------------------
            // TODO: Add default regularization
            KKT( A, x, z, J, false );
            KKTRHS( rc, rb, rmu, z, d );
            const Real pivTol = MaxNorm(J)*epsilon;
            // Do not use any a priori regularization
            Zeros( reg, m+2*n, 1 );

            // Factor the KKT system using dynamic regularization
            // --------------------------------------------------
            if( numIts == 0 )
            {
                NestedDissection( J.LockedGraph(), map, rootSep, info );
                InvertMap( map, invMap );
            }
            JFront.Pull( J, map, info );
            regCandNodal.Pull( invMap, info, regCand );
            regNodal.Pull( invMap, info, reg );
            RegularizedQSDLDL
            ( info, JFront, pivTol, regCandNodal, regNodal, aPriori, LDL_1D );
            regNodal.Push( invMap, info, reg );

            // Compute the proposed step from the regularized KKT system
            // ---------------------------------------------------------
            const Int numLargeRefines = reg_qsd_ldl::SolveAfter
            ( J, reg, invMap, info, JFront, d, 
              REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
            if( numLargeRefines > 3 && !increasedReg )
            {
                Scale( Real(10), regCand );
                increasedReg = true;
            }
            ExpandSolution( m, n, d, dx, dy, dz );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the "augmented" KKT system
            // ------------------------------------
            // TODO: Add default regularization
            AugmentedKKT( A, x, z, J, false );
            AugmentedKKTRHS( x, rc, rb, rmu, d );
            const Real pivTol = MaxNorm(J)*epsilon;
            // Do not use any a priori regularization
            Zeros( reg, m+n, 1 );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            if( ctrl.primalInitialized && ctrl.dualInitialized && numIts == 0 )
            {
                NestedDissection( J.LockedGraph(), map, rootSep, info );
                InvertMap( map, invMap );
            }
            JFront.Pull( J, map, info );
            regCandNodal.Pull( invMap, info, regCand );
            regNodal.Pull( invMap, info, reg );
            RegularizedQSDLDL
            ( info, JFront, pivTol, regCandNodal, regNodal, aPriori, LDL_1D );
            regNodal.Push( invMap, info, reg );

            // Compute the proposed step from the regularized KKT system
            // ---------------------------------------------------------
            const Int numLargeRefines = reg_qsd_ldl::SolveAfter
            ( J, reg, invMap, info, JFront, d, 
              REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
            if( numLargeRefines > 3 && !increasedReg )
            {
                Scale( Real(10), regCand );
                increasedReg = true;
            }
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }
        else if( ctrl.system == NORMAL_KKT )
        {
            // Construct the reduced KKT system, J dy = d
            // ------------------------------------------
            // NOTE: Explicit symmetry is currently required for both METIS and
            //       for the frontal tree initialization
            NormalKKT( A, x, z, J, false );
            NormalKKTRHS( A, x, z, rc, rb, rmu, dy );

            // Factor the KKT system
            // ---------------------
            if( numIts == 0 )
            {
                NestedDissection( J.LockedGraph(), map, rootSep, info );
                InvertMap( map, invMap );
            }
            JFront.Pull( J, map, info );
            LDL( info, JFront ); 

            // Compute the proposed step
            // -------------------------
            ldl::SolveWithIterativeRefinement
            ( J, invMap, info, JFront, dy, minReductionFactor, maxRefineIts );
            ExpandNormalSolution( A, c, x, z, rc, rmu, dx, dy, dz );
        }
        else
            LogicError("Invalid KKT system choice");
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Multiply( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Axpy( Real(-1), dz, dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        const Real rmuNrm2 = Nrm2( rmu );
        dzError = rmu;
        prod = dz;
        DiagonalScale( LEFT, NORMAL, x, prod );
        Axpy( Real(1), prod, dzError );
        prod = dx;
        DiagonalScale( LEFT, NORMAL, z, prod );
        Axpy( Real(1), prod, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        // TODO: Also compute and print the residuals with regularization

        if( ctrl.print )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = " 
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = " 
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_mu ||_2) = " 
                 << dzErrorNrm2/(1+rmuNrm2) << endl;
#endif

        // Take a step in the computed direction
        // =====================================
        const Real alphaPrimal = MaxStepInPositiveCone( x, dx, Real(1) );
        const Real alphaDual = MaxStepInPositiveCone( z, dz, Real(1) );
        const Real alphaMax = Min(alphaPrimal,alphaDual);
        if( ctrl.print )
            cout << "alphaMax = " << alphaMax << endl;
        const Real alpha = 
          IPFLineSearch
          ( A, b, c, x, y, z, dx, dy, dz, 
            Real(0.99)*alphaMax,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.lineSearchCtrl );
        if( ctrl.print )
            cout << "  alpha = " << alpha << endl; 
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
    }

    // Unequilibrate the LP
    DiagonalSolve( LEFT, NORMAL, dCol, x );
    DiagonalSolve( LEFT, NORMAL, dRow, y );
    DiagonalScale( LEFT, NORMAL, dCol, z );
}

template<typename Real>
void IPF
( const DistSparseMatrix<Real>& APre, 
  const DistMultiVec<Real>& bPre,     const DistMultiVec<Real>& cPre,
        DistMultiVec<Real>& x,              DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::direct::IPF"))    
    mpi::Comm comm = APre.Comm();
    const int commRank = mpi::Rank(comm);
    const Real epsilon = lapack::MachineEpsilon<Real>();

    // Equilibrate the LP by diagonally scaling A
    auto A = APre;
    DistMultiVec<Real> dRow(comm), dCol(comm);
    GeomEquil( A, dRow, dCol );
    const Int m = A.Height();
    const Int n = A.Width();
    auto b = bPre;
    auto c = cPre;
    DiagonalSolve( LEFT, NORMAL, dRow, b );
    DiagonalSolve( LEFT, NORMAL, dCol, c );

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );

    DistMap map, invMap;
    DistSymmNodeInfo info;
    DistSeparator rootSep;
    // The initialization involves an augmented KKT system, and so we can
    // only reuse the factorization metadata if the this IPM is using the
    // augmented formulation
    // TODO: Expose this as a parameter to IPFCtrl
    const bool standardShift = true;
    if( ctrl.system == AUGMENTED_KKT )
    {
        Initialize
        ( A, b, c, x, y, z, map, invMap, rootSep, info,
          ctrl.primalInitialized, ctrl.dualInitialized, standardShift,
          ctrl.print );
    }
    else
    {
        DistMap augMap, augInvMap;
        DistSymmNodeInfo augInfo;
        DistSeparator augRootSep;
        Initialize
        ( A, b, c, x, y, z, augMap, augInvMap, augRootSep, augInfo, 
          ctrl.primalInitialized, ctrl.dualInitialized, standardShift,
          ctrl.print );
    }

    DistSparseMatrix<Real> J(comm);
    DistSymmFront<Real> JFront;
    DistMultiVec<Real> d(comm),
                       rc(comm), rb(comm), rmu(comm), 
                       dx(comm), dy(comm), dz(comm);

    DistMultiVec<Real> regCand(comm), reg(comm);
    // TODO: Dynamically modify these values in the manner suggested by 
    //       Altman and Gondzio based upon the number of performed steps of
    //       iterative refinement
    if( ctrl.system == FULL_KKT )
    {
        const Real regMagPrimal = Pow(epsilon,Real(0.75));
        const Real regMagLagrange = Pow(epsilon,Real(0.5));
        const Real regMagDual = Pow(epsilon,Real(0.5));
        regCand.Resize( m+2*n, 1 );
        for( Int iLoc=0; iLoc<regCand.LocalHeight(); ++iLoc )
        {
            const Int i = regCand.GlobalRow(iLoc);
            if( i < n )
                regCand.SetLocal( iLoc, 0, regMagPrimal );
            else if( i < n+m )
                regCand.SetLocal( iLoc, 0, -regMagLagrange );
            else
                regCand.SetLocal( iLoc, 0, -regMagDual );
        }
    }
    else if( ctrl.system == AUGMENTED_KKT )
    {
        const Real regMagPrimal = Pow(epsilon,Real(0.75));
        const Real regMagLagrange = Pow(epsilon,Real(0.5));
        regCand.Resize( n+m, 1 );
        for( Int iLoc=0; iLoc<regCand.LocalHeight(); ++iLoc )
        {
            const Int i = regCand.GlobalRow(iLoc);
            if( i < n )
                regCand.SetLocal( iLoc, 0, regMagPrimal );
            else
                regCand.SetLocal( iLoc, 0, -regMagLagrange );
        }
    }
    DistMultiVecNode<Real> regCandNodal, regNodal;
    bool increasedReg = false;

#ifndef EL_RELEASE
    DistMultiVec<Real> dxError(comm), dyError(comm), dzError(comm), prod(comm);
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = NumNonPositive( x );
        const Int zNumNonPos = NumNonPositive( z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        Scale( Real(-1), rb );
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Axpy( Real(-1), z, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // Now check the pieces
        // --------------------
        if( ctrl.print && commRank == 0 )
            cout << " iter " << numIts << ":\n"
                 << "  |primal - dual| / (1 + |primal|) = "
                 << objConv << "\n"
                 << "  || r_b ||_2 / (1 + || b ||_2)   = "
                 << rbConv << "\n"
                 << "  || r_c ||_2 / (1 + || c ||_2)   = "
                 << rcConv << endl;
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Compute the duality measure and r_mu = x o z - tau e
        // ====================================================
        const Real mu = Dot(x,z) / n;
        rmu = z;
        DiagonalScale( LEFT, NORMAL, x, rmu );
        Shift( rmu, -ctrl.centering*mu );

        // Compute the search direction
        // ============================
        // TODO: Expose these as control parameters
        const Real minReductionFactor = 2;
        const Int maxRefineIts = 50;
        bool aPriori = true;
        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // -----------------------------
            // TODO: Add default regularization
            KKT( A, x, z, J, false );
            KKTRHS( rc, rb, rmu, z, d );
            const Real pivTol = MaxNorm(J)*epsilon;
            // Do not use any a priori regularization
            Zeros( reg, m+2*n, 1 );

            // Factor the KKT system using dynamic regularization
            // --------------------------------------------------
            if( numIts == 0 )
            {
                NestedDissection( J.LockedDistGraph(), map, rootSep, info );
                InvertMap( map, invMap );
            }
            JFront.Pull( J, map, rootSep, info );
            regCandNodal.Pull( invMap, info, regCand );
            regNodal.Pull( invMap, info, reg );
            RegularizedQSDLDL
            ( info, JFront, pivTol, regCandNodal, regNodal, aPriori, LDL_1D );
            regNodal.Push( invMap, info, reg );

            // Compute the proposed step from the regularized KKT system
            // ---------------------------------------------------------
            const Int numLargeRefines = reg_qsd_ldl::SolveAfter
            ( J, reg, invMap, info, JFront, d, 
              REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
            if( numLargeRefines > 3 && !increasedReg )
            {
                Scale( Real(10), regCand );
                increasedReg = true;
            }
            ExpandSolution( m, n, d, dx, dy, dz );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the "augmented" KKT system
            // ------------------------------------
            // TODO: Add default regularization
            AugmentedKKT( A, x, z, J, false );
            AugmentedKKTRHS( x, rc, rb, rmu, d );
            const Real pivTol = MaxNorm(J)*epsilon;
            // Do not use any a priori regularization
            Zeros( reg, m+n, 1 );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            if( ctrl.primalInitialized && ctrl.dualInitialized && numIts == 0 )
            {
                NestedDissection( J.LockedDistGraph(), map, rootSep, info );
                InvertMap( map, invMap );
            }
            JFront.Pull( J, map, rootSep, info );
            regCandNodal.Pull( invMap, info, regCand );
            regNodal.Pull( invMap, info, reg );
            RegularizedQSDLDL
            ( info, JFront, pivTol, regCandNodal, regNodal, aPriori, LDL_1D );
            regNodal.Push( invMap, info, reg );

            // Compute the proposed step from the regularized KKT system
            // ---------------------------------------------------------
            const Int numLargeRefines = reg_qsd_ldl::SolveAfter
            ( J, reg, invMap, info, JFront, d, 
              REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
            if( numLargeRefines > 3 && !increasedReg )
            {
                Scale( Real(10), regCand );
                increasedReg = true;
            }
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }
        else if( ctrl.system == NORMAL_KKT )
        {
            // Construct the reduced KKT system, J dy = d
            // ------------------------------------------
            // NOTE: Explicit symmetry is currently required for both METIS and
            //       for the frontal tree initialization
            NormalKKT( A, x, z, J, false );
            NormalKKTRHS( A, x, z, rc, rb, rmu, dy );

            // Factor the KKT system
            // ---------------------
            if( numIts == 0 )
            {
                NestedDissection( J.LockedDistGraph(), map, rootSep, info );
                InvertMap( map, invMap );
            }
            JFront.Pull( J, map, rootSep, info );
            LDL( info, JFront, LDL_INTRAPIV_1D ); 

            // Compute the proposed step
            // -------------------------
            ldl::SolveWithIterativeRefinement
            ( J, invMap, info, JFront, dy, minReductionFactor, maxRefineIts );
            ExpandNormalSolution( A, c, x, z, rc, rmu, dx, dy, dz );
        }
        else
            LogicError("Invalid KKT system choice");
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Multiply( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Axpy( Real(-1), dz, dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        const Real rmuNrm2 = Nrm2( rmu );
        dzError = rmu;
        prod = dz;
        DiagonalScale( LEFT, NORMAL, x, prod );
        Axpy( Real(1), prod, dzError );
        prod = dx;
        DiagonalScale( LEFT, NORMAL, z, prod );
        Axpy( Real(1), prod, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        // TODO: Also compute and print the residuals with regularization

        if( ctrl.print && commRank == 0 )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = " 
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = " 
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_mu ||_2) = " 
                 << dzErrorNrm2/(1+rmuNrm2) << endl;
#endif

        // Take a step in the computed direction
        // =====================================
        const Real alphaPrimal = MaxStepInPositiveCone( x, dx, Real(1) );
        const Real alphaDual = MaxStepInPositiveCone( z, dz, Real(1) );
        const Real alphaMax = Min(alphaPrimal,alphaDual);
        if( ctrl.print && commRank == 0 )
            cout << "alphaMax = " << alphaMax << endl;
        const Real alpha = 
          IPFLineSearch
          ( A, b, c, x, y, z, dx, dy, dz, 
            Real(0.99)*alphaMax,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.lineSearchCtrl );
        if( ctrl.print && commRank == 0 )
            cout << "  alpha = " << alpha << endl; 
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
    }

    // Unequilibrate the LP
    DiagonalSolve( LEFT, NORMAL, dCol, x );
    DiagonalSolve( LEFT, NORMAL, dRow, y );
    DiagonalScale( LEFT, NORMAL, dCol, z );
}

#define PROTO(Real) \
  template void IPF \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
          Matrix<Real>& x,       Matrix<Real>& y, \
          Matrix<Real>& z, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
          AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
          Matrix<Real>& x,             Matrix<Real>& y, \
          Matrix<Real>& z, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x,           DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const IPFCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace lp
} // namespace El
