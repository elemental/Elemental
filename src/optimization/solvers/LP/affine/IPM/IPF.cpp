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
// which corresponds to G = -I and h = 0, using a simple Infeasible Path 
// Following (IPF) scheme. 
//
// NOTE: This routine should only be used for academic purposes, as the 
//       Mehrotra alternative typically requires an order of magnitude fewer 
//       iterations.

template<typename Real>
void IPF
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,       Matrix<Real>& y, 
        Matrix<Real>& z,       Matrix<Real>& s,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::affine::IPF"))    
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    // TODO: Equilibrate LP using GeomEquil on [A;G] here

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    // TODO: Expose this as a parameter of IPFCtrl
    const bool standardShift = true;
    Initialize
    ( A, G, b, c, h, x, y, z, s, 
      ctrl.primalInitialized, ctrl.dualInitialized, standardShift );

    Matrix<Real> J, d,
                 rmu, rc, rb, rh,
                 dx, dy, dz, ds;
#ifndef EL_RELEASE
    Matrix<Real> dxError, dyError, dzError;
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = NumNonPositive( s );
        const Int zNumNonPos = NumNonPositive( z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) - Dot(h,z);
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
        Gemv( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h;
        Scale( Real(-1), rh );
        Gemv( NORMAL, Real(1), G, x, Real(1), rh );
        Axpy( Real(1), s, rh ); 
        const Real rhNrm2 = Nrm2( rh );
        const Real rhConv = rhNrm2 / (Real(1)+hNrm2);
        // Now check the pieces
        // --------------------
        if( ctrl.print )
            cout << " iter " << numIts << ":\n"
                 << "  |primal - dual| / (1 + |primal|) = "
                 << objConv << "\n"
                 << "  || r_b ||_2 / (1 + || b ||_2)   = "
                 << rbConv << "\n"
                 << "  || r_c ||_2 / (1 + || c ||_2)   = "
                 << rcConv << "\n"
                 << "  || r_h ||_2 / (1 + || h ||_2)   = "
                 << rhConv << endl;
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && 
            rcConv  <= ctrl.tol && rhConv <= ctrl.tol )
            break;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Form the residual for the scaled equation, s o z - sigma mu e
        // =============================================================
        const Real mu = Dot(s,z) / k;
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );
        Shift( rmu, -ctrl.centering*mu );

        // Construct the full KKT system
        // =============================
        KKT( A, G, s, z, J );
        KKTRHS( rc, rb, rh, rmu, z, d );

        // Compute the proposed step from the KKT system
        // =============================================
        SymmetricSolve( LOWER, NORMAL, J, d );
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dz, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dx, Real(1), dzError );
        Axpy( Real(1), ds, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        // TODO: dmuError

        if( ctrl.print )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_mu ||_2) = "
                 << dzErrorNrm2/(1+rhNrm2) << endl;
#endif

        // Take a step in the computed direction
        // =====================================
        const Real alphaPrimal = MaxStepInPositiveCone( s, ds, Real(1) );
        const Real alphaDual = MaxStepInPositiveCone( z, dz, Real(1) );
        const Real alphaMax = Min(alphaPrimal,alphaDual);
        if( ctrl.print )
            cout << "alphaMax = " << alphaMax << endl;
        const Real alpha =
          IPFLineSearch
          ( A, G, b, c, h, x, y, z, s, dx, dy, dz, ds,
            Real(0.99)*alphaMax,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.tol*(1+hNrm2),
            ctrl.lineSearchCtrl );
        if( ctrl.print )
            cout << "  alpha = " << alpha << endl;
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
        Axpy( alpha, ds, s );
    }
    // TODO: Unequilibrate LP here
}

template<typename Real>
void IPF
( const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& GPre,
  const AbstractDistMatrix<Real>& bPre, const AbstractDistMatrix<Real>& cPre,
  const AbstractDistMatrix<Real>& hPre,
        AbstractDistMatrix<Real>& xPre,       AbstractDistMatrix<Real>& yPre, 
        AbstractDistMatrix<Real>& zPre,       AbstractDistMatrix<Real>& sPre,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::affine::IPF"))    

    ProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = 0;
    proxCtrl.rowAlign = 0;
    auto APtr = ReadProxy<Real,MC,MR>(&APre,proxCtrl);      auto& A = *APtr;
    auto bPtr = ReadProxy<Real,MC,MR>(&bPre,proxCtrl);      auto& b = *bPtr;
    auto cPtr = ReadProxy<Real,MC,MR>(&cPre,proxCtrl);      auto& c = *cPtr;
    auto hPtr = ReadProxy<Real,MC,MR>(&hPre,proxCtrl);      auto& h = *hPtr;
    auto GPtr = ReadProxy<Real,MC,MR>(&GPre,proxCtrl);      auto& G = *GPtr;
    // NOTE: {x,s} do not need to be read proxies when !ctrl.primalInitialized
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,proxCtrl); auto& x = *xPtr;
    auto sPtr = ReadWriteProxy<Real,MC,MR>(&sPre,proxCtrl); auto& s = *sPtr;
    // NOTE: {y,z} do not need to be read proxies when !ctrl.dualInitialized
    auto yPtr = ReadWriteProxy<Real,MC,MR>(&yPre,proxCtrl); auto& y = *yPtr;
    auto zPtr = ReadWriteProxy<Real,MC,MR>(&zPre,proxCtrl); auto& z = *zPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Grid& grid = A.Grid();
    const Int commRank = A.Grid().Rank();

    // TODO: Equilibrate LP using GeomEquil on [A;G] here

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    // TODO: Expose this as a parameter of IPFCtrl
    const bool standardShift = true;
    Initialize
    ( A, G, b, c, h, x, y, z, s, 
      ctrl.primalInitialized, ctrl.dualInitialized, standardShift );

    DistMatrix<Real> J(grid), d(grid), 
                     rc(grid), rb(grid), rh(grid), rmu(grid),
                     dx(grid), dy(grid), dz(grid), ds(grid);
    ds.AlignWith( s );
    dz.AlignWith( s );
    rmu.AlignWith( s );
#ifndef EL_RELEASE
    DistMatrix<Real> dxError(grid), dyError(grid), dzError(grid);
    dzError.AlignWith( s );
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = NumNonPositive( s );
        const Int zNumNonPos = NumNonPositive( z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) - Dot(h,z);
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
        Gemv( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h;
        Scale( Real(-1), rh );
        Gemv( NORMAL, Real(1), G, x, Real(1), rh );
        Axpy( Real(1), s, rh ); 
        const Real rhNrm2 = Nrm2( rh );
        const Real rhConv = rhNrm2 / (Real(1)+hNrm2);
        // Now check the pieces
        // --------------------
        if( ctrl.print && commRank == 0 )
            cout << " iter " << numIts << ":\n"
                 << "  |primal - dual| / (1 + |primal|) = "
                 << objConv << "\n"
                 << "  || r_b ||_2 / (1 + || b ||_2)   = "
                 << rbConv << "\n"
                 << "  || r_c ||_2 / (1 + || c ||_2)   = "
                 << rcConv << "\n"
                 << "  || r_h ||_2 / (1 + || h ||_2)   = "
                 << rhConv << endl;
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && 
            rcConv  <= ctrl.tol && rhConv <= ctrl.tol )
            break;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Form the residual for the scaled equation, s o z - sigma mu e
        // =============================================================
        const Real mu = Dot(s,z) / k;
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );
        Shift( rmu, -ctrl.centering*mu );

        // Construct the full KKT system
        // =============================
        KKT( A, G, s, z, J );
        KKTRHS( rc, rb, rh, rmu, z, d );

        // Compute the proposed step from the KKT system
        // =============================================
        SymmetricSolve( LOWER, NORMAL, J, d );
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dz, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dx, Real(1), dzError );
        Axpy( Real(1), ds, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        // TODO: dmuError

        if( ctrl.print && commRank == 0 )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_mu ||_2) = "
                 << dzErrorNrm2/(1+rhNrm2) << endl;
#endif

        // Take a step in the computed direction
        // =====================================
        const Real alphaPrimal = MaxStepInPositiveCone( s, ds, Real(1) );
        const Real alphaDual = MaxStepInPositiveCone( z, dz, Real(1) );
        const Real alphaMax = Min(alphaPrimal,alphaDual);
        if( ctrl.print && commRank == 0 )
            cout << "alphaMax = " << alphaMax << endl;
        const Real alpha =
          IPFLineSearch
          ( A, G, b, c, h, x, y, z, s, dx, dy, dz, ds,
            Real(0.99)*alphaMax,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.tol*(1+hNrm2),
            ctrl.lineSearchCtrl );
        if( ctrl.print && commRank == 0 )
            cout << "  alpha = " << alpha << endl;
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
        Axpy( alpha, ds, s );
    }
    // TODO: Unequilibrate LP here
}

template<typename Real>
void IPF
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,             Matrix<Real>& y, 
        Matrix<Real>& z,             Matrix<Real>& s,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::affine::IPF"))    

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Real epsilon = lapack::MachineEpsilon<Real>();

    // TODO: Equilibrate LP using GeomEquil on [A;G] here

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    vector<Int> map, invMap;
    SymmNodeInfo info;
    Separator rootSep;
    // TODO: Expose this as a parameter of IPFCtrl
    const bool standardShift = true;
    Initialize
    ( A, G, b, c, h, x, y, z, s, map, invMap, rootSep, info, 
      ctrl.primalInitialized, ctrl.dualInitialized, standardShift, ctrl.print );

    SparseMatrix<Real> J;
    SymmFront<Real> JFront;
    Matrix<Real> d,
                 rc, rb, rh, rmu,
                 dx, dy, dz, ds;

    // TODO: Dynamically modify these values in the manner suggested by 
    //       Altman and Gondzio based upon the number of performed steps of
    //       iterative refinement
    const Real regMagPrimal = Pow(epsilon,Real(0.75));
    const Real regMagLagrange = Pow(epsilon,Real(0.5));
    const Real regMagDual = Pow(epsilon,Real(0.5));
    Matrix<Real> regCand, reg;
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
    MatrixNode<Real> regCandNodal, regNodal;
    bool increasedReg = false;

#ifndef EL_RELEASE
    Matrix<Real> dxError, dyError, dzError;
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = NumNonPositive( s );
        const Int zNumNonPos = NumNonPositive( z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) - Dot(h,z);
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
        Multiply( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h;
        Scale( Real(-1), rh );
        Multiply( NORMAL, Real(1), G, x, Real(1), rh );
        Axpy( Real(1), s, rh ); 
        const Real rhNrm2 = Nrm2( rh );
        const Real rhConv = rhNrm2 / (Real(1)+hNrm2);
        // Now check the pieces
        // --------------------
        if( ctrl.print )
            cout << " iter " << numIts << ":\n"
                 << "  |primal - dual| / (1 + |primal|) = "
                 << objConv << "\n"
                 << "  || r_b ||_2 / (1 + || b ||_2)   = "
                 << rbConv << "\n"
                 << "  || r_c ||_2 / (1 + || c ||_2)   = "
                 << rcConv << "\n"
                 << "  || r_h ||_2 / (1 + || h ||_2)   = "
                 << rhConv << endl;
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && 
            rcConv  <= ctrl.tol && rhConv <= ctrl.tol )
            break;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Form the residual for the scaled equation, s o z - sigma mu e
        // =============================================================
        const Real mu = Dot(s,z) / k;
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );
        Shift( rmu, -ctrl.centering*mu );

        const Real minReductionFactor = 2;
        const Int maxRefineIts = 50;
        bool aPriori = true;
        {
            // TODO: Add default regularization
            KKT( A, G, s, z, J, false );
            KKTRHS( rc, rb, rh, rmu, z, d );
            const Real pivTol = MaxNorm(J)*epsilon;
            // Do not use any a priori regularization
            Zeros( reg, m+n+k, 1 );

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

            const Int numLargeRefines = reg_qsd_ldl::SolveAfter
            ( J, reg, invMap, info, JFront, d, 
              REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
            if( numLargeRefines > 3 && !increasedReg )
            {
                Scale( Real(10), regCand );
                increasedReg = true;
            }
            ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        }
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Multiply( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Multiply( TRANSPOSE, Real(1), G, dz, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Multiply( NORMAL, Real(1), G, dx, Real(1), dzError );
        Axpy( Real(1), ds, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        // TODO: dmuError
        // TODO: Also compute and print the residuals with regularization

        if( ctrl.print )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_mu ||_2) = "
                 << dzErrorNrm2/(1+rhNrm2) << endl;
#endif

        // Take a step in the computed direction
        // =====================================
        const Real alphaPrimal = MaxStepInPositiveCone( s, ds, Real(1) );
        const Real alphaDual = MaxStepInPositiveCone( z, dz, Real(1) );
        const Real alphaMax = Min(alphaPrimal,alphaDual);
        if( ctrl.print )
            cout << "alphaMax = " << alphaMax << endl;
        const Real alpha =
          IPFLineSearch
          ( A, G, b, c, h, x, y, z, s, dx, dy, dz, ds,
            Real(0.99)*alphaMax,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.tol*(1+hNrm2),
            ctrl.lineSearchCtrl );
        if( ctrl.print )
            cout << "  alpha = " << alpha << endl;
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
        Axpy( alpha, ds, s );
    }
    // TODO: Unequilibrate LP here
}

template<typename Real>
void IPF
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,           DistMultiVec<Real>& s,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::affine::IPF"))    

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);
    const Real epsilon = lapack::MachineEpsilon<Real>();

    // TODO: Equilibrate LP using GeomEquil on [A;G] here

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    DistMap map, invMap;
    DistSymmNodeInfo info;
    DistSeparator rootSep;
    // TODO: Expose this as a parameter of IPFCtrl
    const bool standardShift = true;
    Initialize
    ( A, G, b, c, h, x, y, z, s, map, invMap, rootSep, info, 
      ctrl.primalInitialized, ctrl.dualInitialized, standardShift, ctrl.print );

    DistSparseMatrix<Real> J(comm);
    DistSymmFront<Real> JFront;
    DistMultiVec<Real> d(comm),
                       rc(comm), rb(comm), rh(comm), rmu(comm),
                       dx(comm), dy(comm), dz(comm), ds(comm);

    // TODO: Dynamically modify these values in the manner suggested by 
    //       Altman and Gondzio based upon the number of performed steps of
    //       iterative refinement
    const Real regMagPrimal = Pow(epsilon,Real(0.75));
    const Real regMagLagrange = Pow(epsilon,Real(0.5));
    const Real regMagDual = Pow(epsilon,Real(0.5));
    DistMultiVec<Real> regCand(comm), reg(comm);
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
    DistMultiVecNode<Real> regCandNodal, regNodal;
    bool increasedReg = false;

#ifndef EL_RELEASE
    DistMultiVec<Real> dxError(comm), dyError(comm), dzError(comm);
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = NumNonPositive( s );
        const Int zNumNonPos = NumNonPositive( z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) - Dot(h,z);
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
        Multiply( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h;
        Scale( Real(-1), rh );
        Multiply( NORMAL, Real(1), G, x, Real(1), rh );
        Axpy( Real(1), s, rh ); 
        const Real rhNrm2 = Nrm2( rh );
        const Real rhConv = rhNrm2 / (Real(1)+hNrm2);
        // Now check the pieces
        // --------------------
        if( ctrl.print && commRank == 0 )
            cout << " iter " << numIts << ":\n"
                 << "  |primal - dual| / (1 + |primal|) = "
                 << objConv << "\n"
                 << "  || r_b ||_2 / (1 + || b ||_2)   = "
                 << rbConv << "\n"
                 << "  || r_c ||_2 / (1 + || c ||_2)   = "
                 << rcConv << "\n"
                 << "  || r_h ||_2 / (1 + || h ||_2)   = "
                 << rhConv << endl;
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && 
            rcConv  <= ctrl.tol && rhConv <= ctrl.tol )
            break;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Form the residual for the scaled equation, s o z - sigma mu e
        // =============================================================
        const Real mu = Dot(s,z) / k;
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );
        Shift( rmu, -ctrl.centering*mu );

        const Real minReductionFactor = 2;
        const Int maxRefineIts = 50;
        bool aPriori = true;
        {
            // TODO: Add default regularization
            KKT( A, G, s, z, J, false );
            KKTRHS( rc, rb, rh, rmu, z, d );
            const Real pivTol = MaxNorm(J)*epsilon;
            // Do not use any a priori regularization
            Zeros( reg, m+n+k, 1 );

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

            const Int numLargeRefines = reg_qsd_ldl::SolveAfter
            ( J, reg, invMap, info, JFront, d, 
              REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
            if( numLargeRefines > 3 && !increasedReg )
            {
                Scale( Real(10), regCand );
                increasedReg = true;
            }
            ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        }
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Multiply( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Multiply( TRANSPOSE, Real(1), G, dz, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Multiply( NORMAL, Real(1), G, dx, Real(1), dzError );
        Axpy( Real(1), ds, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        // TODO: dmuError
        // TODO: Also compute and print the residuals with regularization

        if( ctrl.print && commRank == 0 )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_mu ||_2) = "
                 << dzErrorNrm2/(1+rhNrm2) << endl;
#endif

        // Take a step in the computed direction
        // =====================================
        const Real alphaPrimal = MaxStepInPositiveCone( s, ds, Real(1) );
        const Real alphaDual = MaxStepInPositiveCone( z, dz, Real(1) );
        const Real alphaMax = Min(alphaPrimal,alphaDual);
        if( ctrl.print && commRank == 0 )
            cout << "alphaMax = " << alphaMax << endl;
        const Real alpha =
          IPFLineSearch
          ( A, G, b, c, h, x, y, z, s, dx, dy, dz, ds,
            Real(0.99)*alphaMax,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.tol*(1+hNrm2),
            ctrl.lineSearchCtrl );
        if( ctrl.print && commRank == 0 )
            cout << "  alpha = " << alpha << endl;
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
        Axpy( alpha, ds, s );
    }
    // TODO: Unequilibrate LP here
}

#define PROTO(Real) \
  template void IPF \
  ( const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,       Matrix<Real>& y, \
          Matrix<Real>& z,       Matrix<Real>& s, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
          AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,             Matrix<Real>& y, \
          Matrix<Real>& z,             Matrix<Real>& s, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x,           DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z,           DistMultiVec<Real>& s, \
    const IPFCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace affine
} // namespace lp
} // namespace El
