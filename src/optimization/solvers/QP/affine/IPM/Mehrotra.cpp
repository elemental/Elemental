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
namespace qp {
namespace affine {

// The following solves a pair of quadratic programs in "affine" conic form:
//
//   min (1/2) x^T Q x + c^T x
//   s.t. A x = b, G x + s = h, s >= 0,
//
//   max (1/2) (A^T y + G^T z + c)^T pinv(Q) (A^T y + G^T z + c) - b^T y - h^T z
//   s.t. A^T y + G^T z + c in range(Q), z >= 0,
//
// as opposed to the more specific "direct" conic form:
//
//   min (1/2) x^T Q x + c^T x
//   s.t. A x = b, x >= 0,
//
//   max (1/2) (A^T y - z + c)^T pinv(Q) (A^T y - z + c) - b^T y
//   s.t. A^T y - z + c in range(Q), z >= 0,  
//
// which corresponds to G = -I and h = 0, using a Mehrotra Predictor-Corrector 
// scheme.
//

template<typename Real>
void Mehrotra
( const Matrix<Real>& Q,
  const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,       Matrix<Real>& y, 
        Matrix<Real>& z,       Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::Mehrotra"))    
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    // TODO: Equilibrate the QP here by calling GeomEquil on [A;G]

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    // TODO: Expose this as a parameter to MehrotraCtrl
    const bool standardShift = true;
    Initialize
    ( Q, A, G, b, c, h, x, y, z, s, 
      ctrl.primalInitialized, ctrl.dualInitialized, standardShift );

    Matrix<Real> J, d,
                 rmu,   rc,    rb,    rh,
                 dxAff, dyAff, dzAff, dsAff,
                 dx,    dy,    dz,    ds;
    Matrix<Real> dSub;
    Matrix<Int> p;
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
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        Zeros( d, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj =  xTQx/2 + Dot(c,x);
        const Real dualObj = -xTQx/2 - Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b; Scale( Real(-1), rb );
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Hemv( LOWER,     Real(1), Q, x, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h; Scale( Real(-1), rh );
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

        // r_mu := s o z
        // =============
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );

        // Compute the affine search direction
        // ===================================
        // Construct the full KKT system
        // -----------------------------
        KKT( Q, A, G, s, z, J );
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        LDL( J, dSub, p, false );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Hemv( LOWER,     Real(1), Q, dxAff, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
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

        // Compute the max affine [0,1]-step which keeps s and z in the cone
        // =================================================================
        const Real alphaAffPri = MaxStepInPositiveCone( s, dsAff, Real(1) );
        const Real alphaAffDual = MaxStepInPositiveCone( z, dzAff, Real(1) );
        if( ctrl.print )
            cout << "  alphaAffPri = " << alphaAffPri
                 << ", alphaAffDual = " << alphaAffDual << endl;

        // Compute what the new duality measure would become
        // =================================================
        const Real mu = Dot(s,z) / k;
        // NOTE: dz and ds are used as temporaries
        ds = s;
        dz = z;
        Axpy( alphaAffPri,  dsAff, ds );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(ds,dz) / k;

        // Compute a centrality parameter using Mehrotra's formula
        // =======================================================
        // TODO: Allow the user to override this function
        const Real sigma = Pow(muAff/mu,Real(3));
        if( ctrl.print )
            cout << "  muAff = " << muAff
                 << ", mu = " << mu
                 << ", sigma = " << sigma << endl;

        // Solve for the centering-corrector
        // =================================
        Zeros( rc, n, 1 ); 
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        // r_mu := dsAff o dzAff - sigma*mu
        // --------------------------------
        rmu = dzAff;
        DiagonalScale( LEFT, NORMAL, dsAff, rmu );
        Shift( rmu, -sigma*mu );
        // Construct the new full KKT RHS
        // ------------------------------
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        // TODO: Residual checks for center-corrector

        // Add in the affine search direction
        // ==================================
        Axpy( Real(1), dxAff, dx );
        Axpy( Real(1), dyAff, dy );
        Axpy( Real(1), dzAff, dz );
        Axpy( Real(1), dsAff, ds );

        // Compute max [0,1/maxStepRatio] step which keeps s and z in the cone
        // ===================================================================
        Real alphaPri = MaxStepInPositiveCone( s, ds, 1/ctrl.maxStepRatio );
        Real alphaDual = MaxStepInPositiveCone( z, dz, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.print )
            cout << "  alphaPri = " << alphaPri
                 << ", alphaDual = " << alphaDual << endl;

        // Update the current estimates
        // ============================
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
    }

    // TODO: Unequilibrate the QP here
}

template<typename Real>
void Mehrotra
( const AbstractDistMatrix<Real>& QPre, 
  const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& GPre,
  const AbstractDistMatrix<Real>& bPre, const AbstractDistMatrix<Real>& cPre,
  const AbstractDistMatrix<Real>& hPre,
        AbstractDistMatrix<Real>& xPre,       AbstractDistMatrix<Real>& yPre, 
        AbstractDistMatrix<Real>& zPre,       AbstractDistMatrix<Real>& sPre,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::Mehrotra"))    

    ProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = 0;
    proxCtrl.rowAlign = 0;
    auto QPtr = ReadProxy<Real,MC,MR>(&QPre,proxCtrl);      auto& Q = *QPtr;
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

    // TODO: Equilibrate the QP here by calling GeomEquil on [A;G]

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    // TODO: Expose this as a parameter to MehrotraCtrl and extend to other
    //       QP and LP IPMs
    bool forceSameStep = false;

    // TODO: Expose this as a parameter to MehrotraCtrl
    const bool standardShift = true;
    Initialize
    ( Q, A, G, b, c, h, x, y, z, s, 
      ctrl.primalInitialized, ctrl.dualInitialized, standardShift );

    DistMatrix<Real> J(grid),     d(grid), 
                     rc(grid),    rb(grid),    rh(grid),    rmu(grid),
                     dxAff(grid), dyAff(grid), dzAff(grid), dsAff(grid),
                     dx(grid),    dy(grid),    dz(grid),    ds(grid);
    dsAff.AlignWith( s );
    dzAff.AlignWith( s );
    ds.AlignWith( s );
    dz.AlignWith( s );
    rmu.AlignWith( s );
    DistMatrix<Real> dSub(grid);
    DistMatrix<Int> p(grid);
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
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        Zeros( d, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj =  xTQx/2 + Dot(c,x);
        const Real dualObj = -xTQx/2 - Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b; Scale( Real(-1), rb );
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Hemv( LOWER,     Real(1), Q, x, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h; Scale( Real(-1), rh );
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

        // r_mu := s o z
        // =============
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );

        // Compute the affine search direction
        // ===================================
        // Construct the full KKT system
        // -----------------------------
        KKT( Q, A, G, s, z, J );
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        LDL( J, dSub, p, false );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Hemv( LOWER,     Real(1), Q, dxAff, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
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
 
        // Compute the max affine [0,1]-step which keeps s and z in the cone
        // =================================================================
        const Real alphaAffPri = MaxStepInPositiveCone( s, dsAff, Real(1) );
        const Real alphaAffDual = MaxStepInPositiveCone( z, dzAff, Real(1) );
        if( ctrl.print && commRank == 0 )
            cout << "  alphaAffPri = " << alphaAffPri
                 << ", alphaAffDual = " << alphaAffDual << endl;

        // Compute what the new duality measure would become
        // =================================================
        const Real mu = Dot(s,z) / k;
        // NOTE: dz and ds are used as temporaries
        ds = s;
        dz = z;
        if( forceSameStep )
        {
            const Real alphaAff = Min(alphaAffPri,alphaAffDual);
            Axpy( alphaAff, dsAff, ds );
            Axpy( alphaAff, dzAff, dz );
        }
        else
        {
            Axpy( alphaAffPri,  dsAff, ds );
            Axpy( alphaAffDual, dzAff, dz );
        }
        const Real muAff = Dot(ds,dz) / k;

        // Compute a centrality parameter using Mehrotra's formula
        // =======================================================
        // TODO: Allow the user to override this function
        const Real sigma = Pow(muAff/mu,Real(3));
        if( ctrl.print && commRank == 0 )
            cout << "  muAff = " << muAff
                 << ", mu = " << mu
                 << ", sigma = " << sigma << endl;

        // Solve for the centering-corrector
        // =================================
        Zeros( rc, n, 1 );
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        // r_mu := dsAff o dzAff - sigma*mu
        // --------------------------------
        rmu = dzAff;
        DiagonalScale( LEFT, NORMAL, dsAff, rmu ); 
        Shift( rmu, -sigma*mu );
        // Construct the new full KKT RHS
        // ------------------------------
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        // TODO: Residual checks for center-corrector

        // Add in the affine search direction
        // ==================================
        Axpy( Real(1), dxAff, dx );
        Axpy( Real(1), dyAff, dy );
        Axpy( Real(1), dzAff, dz );
        Axpy( Real(1), dsAff, ds );

        // Compute max [0,1/maxStepRatio] step which keeps s and z in the cone
        // ===================================================================
        Real alphaPri = MaxStepInPositiveCone( s, ds, 1/ctrl.maxStepRatio );
        Real alphaDual = MaxStepInPositiveCone( z, dz, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.print && commRank == 0 )
            cout << "  alphaPri = " << alphaPri
                 << ", alphaDual = " << alphaDual << endl;

        // Update the current estimates
        // ============================
        if( forceSameStep )
        {
            Real alpha = Min(alphaPri,alphaDual);
            Axpy( alpha, dx, x );
            Axpy( alpha, ds, s );
            Axpy( alpha, dy, y );
            Axpy( alpha, dz, z );
        }
        else
        {
            Axpy( alphaPri,  dx, x );
            Axpy( alphaPri,  ds, s );
            Axpy( alphaDual, dy, y );
            Axpy( alphaDual, dz, z );
        }
    }

    // TODO: Unequilibrate the QP here
}

template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& Q,
  const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,             Matrix<Real>& y, 
        Matrix<Real>& z,             Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::Mehrotra"))    
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Real epsilon = lapack::MachineEpsilon<Real>();

    // TODO: Equilibrate the QP here by calling GeomEquil on [A;G]

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    vector<Int> map, invMap;
    SymmNodeInfo info;
    Separator rootSep;
    // TODO: Expose this as a parameter to MehrotraCtrl
    const bool standardShift = true;
    Initialize
    ( Q, A, G, b, c, h, x, y, z, s, map, invMap, rootSep, info, 
      ctrl.primalInitialized, ctrl.dualInitialized, standardShift, ctrl.print );

    SparseMatrix<Real> J;
    SymmFront<Real> JFront;
    Matrix<Real> d,
                 rc,    rb,    rh,    rmu,
                 dxAff, dyAff, dzAff, dsAff,
                 dx,    dy,    dz,    ds;

    Matrix<Real> regCand, reg;
    // TODO: Dynamically modify these values in the manner suggested by 
    //       Altman and Gondzio based upon the number of performed steps of
    //       iterative refinement
    const Real regMagPrimal = Pow(epsilon,Real(0.75));
    const Real regMagLagrange = Pow(epsilon,Real(0.5));
    const Real regMagDual = Pow(epsilon,Real(0.5));
    regCand.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
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
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        Zeros( d, n, 1 );
        // NOTE: The following assumes that Q is explicitly symmetric
        Multiply( NORMAL, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj =  xTQx/2 + Dot(c,x);
        const Real dualObj = -xTQx/2 - Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b; Scale( Real(-1), rb );
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( NORMAL,    Real(1), Q, x, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h; Scale( Real(-1), rh );
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

        // r_mu := s o z
        // =============
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );

        // Compute the affine search direction
        // ===================================
        const Real minReductionFactor = 2;
        const Int maxRefineIts = 50;
        bool aPriori = true;
        Int numLargeAffineRefines = 0;
        {
            // Construct the full KKT system
            // -----------------------------
            // TODO: Add default regularization
            KKT( Q, A, G, s, z, J, false );
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

            numLargeAffineRefines = reg_qsd_ldl::SolveAfter
            ( J, reg, invMap, info, JFront, d,
              REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
            ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );
        }
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Multiply( NORMAL,    Real(1), Q, dxAff, Real(1), dyError );
        Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Multiply( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Multiply( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
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

        // Compute the max affine [0,1]-step which keeps s and z in the cone
        // =================================================================
        const Real alphaAffPri = MaxStepInPositiveCone( s, dsAff, Real(1) );
        const Real alphaAffDual = MaxStepInPositiveCone( z, dzAff, Real(1) );
        if( ctrl.print )
            cout << "  alphaAffPri = " << alphaAffPri
                 << ", alphaAffDual = " << alphaAffDual << endl;

        // Compute what the new duality measure would become
        // =================================================
        const Real mu = Dot(s,z) / k;
        // NOTE: dz and ds are used as temporaries
        ds = s;
        dz = z;
        Axpy( alphaAffPri,  dsAff, ds );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(ds,dz) / k;

        // Compute a centrality parameter using Mehrotra's formula
        // =======================================================
        // TODO: Allow the user to override this function
        const Real sigma = Pow(muAff/mu,Real(3));
        if( ctrl.print )
            cout << "  muAff = " << muAff
                 << ", mu = " << mu
                 << ", sigma = " << sigma << endl;

        // Solve for the centering-corrector
        // =================================
        Zeros( rc, n, 1 );
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        // r_mu := dsAff o dzAff - sigma*mu
        // --------------------------------
        rmu = dzAff;
        DiagonalScale( LEFT, NORMAL, dsAff, rmu );
        Shift( rmu, -sigma*mu );
        // Construct the new full KKT RHS
        // ------------------------------
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        const Int numLargeCorrectorRefines = reg_qsd_ldl::SolveAfter
        ( J, reg, invMap, info, JFront, d,
          REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        if( Max(numLargeAffineRefines,numLargeCorrectorRefines) > 3 && 
            !increasedReg )
        {
            Scale( Real(10), regCand );
            increasedReg = true;
        }

        // Add in the affine search direction
        // ==================================
        Axpy( Real(1), dxAff, dx );
        Axpy( Real(1), dyAff, dy );
        Axpy( Real(1), dzAff, dz );
        Axpy( Real(1), dsAff, ds );

        // Compute max [0,1/maxStepRatio] step which keeps s and z in the cone
        // ===================================================================
        Real alphaPri = MaxStepInPositiveCone( s, ds, 1/ctrl.maxStepRatio );
        Real alphaDual = MaxStepInPositiveCone( z, dz, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.print )
            cout << "  alphaPri = " << alphaPri
                 << ", alphaDual = " << alphaDual << endl;

        // Update the current estimates
        // ============================
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
    }

    // TODO: Unequilibrate the QP here
}

template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& Q,
  const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,           DistMultiVec<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::affine::Mehrotra"))    
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);
    const Real epsilon = lapack::MachineEpsilon<Real>();

    // TODO: Equilibrate the QP here by calling GeomEquil on [A;G]

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    DistMap map, invMap;
    DistSymmNodeInfo info;
    DistSeparator rootSep;
    // TODO: Expose this as a parameter to MehrotraCtrl
    const bool standardShift = true;
    Initialize
    ( Q, A, G, b, c, h, x, y, z, s, map, invMap, rootSep, info, 
      ctrl.primalInitialized, ctrl.dualInitialized, standardShift, ctrl.print );

    DistSparseMatrix<Real> J(comm);
    DistSymmFront<Real> JFront;
    DistMultiVec<Real> d(comm),
                       rc(comm),    rb(comm),    rh(comm),    rmu(comm),
                       dxAff(comm), dyAff(comm), dzAff(comm), dsAff(comm),
                       dx(comm),    dy(comm),    dz(comm),    ds(comm);

    DistMultiVec<Real> regCand(comm), reg(comm);
    // TODO: Dynamically modify these values in the manner suggested by 
    //       Altman and Gondzio based upon the number of performed steps of
    //       iterative refinement
    const Real regMagPrimal = Pow(epsilon,Real(0.75));
    const Real regMagLagrange = Pow(epsilon,Real(0.5));
    const Real regMagDual = Pow(epsilon,Real(0.5));
    regCand.Resize( n+m+k, 1 );
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
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        Zeros( d, n, 1 );
        // NOTE: The following assumes that Q is explicitly symmetric
        Multiply( NORMAL, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj =  xTQx/2 + Dot(c,x);
        const Real dualObj = -xTQx/2 - Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b; Scale( Real(-1), rb );
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( NORMAL,    Real(1), Q, x, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h; Scale( Real(-1), rh );
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

        // r_mu := s o z
        // =============
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );

        // Compute the affine search direction
        // ===================================
        const Real minReductionFactor = 2;
        const Int maxRefineIts = 50;
        bool aPriori = true;
        Int numLargeAffineRefines = 0;
        {
            // Construct the full KKT system
            // -----------------------------
            // TODO: Add default regularization
            KKT( Q, A, G, s, z, J, false );
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

            numLargeAffineRefines = reg_qsd_ldl::SolveAfter
            ( J, reg, invMap, info, JFront, d,
              REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
            ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );
        }
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dxError = rb;
        Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Multiply( NORMAL,    Real(1), Q, dxAff, Real(1), dyError );
        Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Multiply( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Multiply( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
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

        // Compute the max affine [0,1]-step which keeps s and z in the cone
        // =================================================================
        const Real alphaAffPri = MaxStepInPositiveCone( s, dsAff, Real(1) );
        const Real alphaAffDual = MaxStepInPositiveCone( z, dzAff, Real(1) );
        if( ctrl.print && commRank == 0 )
            cout << "  alphaAffPri = " << alphaAffPri
                 << ", alphaAffDual = " << alphaAffDual << endl;

        // Compute what the new duality measure would become
        // =================================================
        const Real mu = Dot(s,z) / k;
        // NOTE: dz and ds are used as temporaries
        ds = s;
        dz = z;
        Axpy( alphaAffPri,  dsAff, ds );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(ds,dz) / k;

        // Compute a centrality parameter using Mehrotra's formula
        // =======================================================
        // TODO: Allow the user to override this function
        const Real sigma = Pow(muAff/mu,Real(3));
        if( ctrl.print && commRank == 0 )
            cout << "  muAff = " << muAff
                 << ", mu = " << mu
                 << ", sigma = " << sigma << endl;

        // Solve for the centering-corrector
        // =================================
        Zeros( rc, n, 1 );
        Zeros( rb, m, 1 );
        Zeros( rh, k, 1 );
        // r_mu := dsAff o dzAff - sigma*mu
        // --------------------------------
        rmu = dzAff;
        DiagonalScale( LEFT, NORMAL, dsAff, rmu );
        Shift( rmu, -sigma*mu );
        // Construct the new full KKT RHS
        // ------------------------------
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        const Int numLargeCorrectorRefines = reg_qsd_ldl::SolveAfter
        ( J, reg, invMap, info, JFront, d,
          REG_REFINE_FGMRES, minReductionFactor, maxRefineIts, ctrl.print );
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        if( Max(numLargeAffineRefines,numLargeCorrectorRefines) > 3 && 
            !increasedReg )
        {
            Scale( Real(10), regCand );
            increasedReg = true;
        }

        // Add in the affine search direction
        // ==================================
        Axpy( Real(1), dxAff, dx );
        Axpy( Real(1), dyAff, dy );
        Axpy( Real(1), dzAff, dz );
        Axpy( Real(1), dsAff, ds );

        // Compute max [0,1/maxStepRatio] step which keeps s and z in the cone
        // ===================================================================
        Real alphaPri = MaxStepInPositiveCone( s, ds, 1/ctrl.maxStepRatio );
        Real alphaDual = MaxStepInPositiveCone( z, dz, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.print && commRank == 0 )
            cout << "  alphaPri = " << alphaPri
                 << ", alphaDual = " << alphaDual << endl;

        // Update the current estimates
        // ============================
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
    }

    // TODO: Unequilibrate the QP here
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,       Matrix<Real>& y, \
          Matrix<Real>& z,       Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AbstractDistMatrix<Real>& Q, \
    const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
          AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,             Matrix<Real>& y, \
          Matrix<Real>& z,             Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x,           DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z,           DistMultiVec<Real>& s, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace dual
} // namespace qp
} // namespace El
