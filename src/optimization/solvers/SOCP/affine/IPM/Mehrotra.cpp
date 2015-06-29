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
namespace socp {
namespace affine {

// The following solves a pair of SOC programs in "affine" conic form:
//
//   min c^T x
//   s.t. A x = b, G x + s = h, s in K,
//
//   max -b^T y - h^T z
//   s.t. A^T y + G^T z + c = 0, z in K,
//
// as opposed to the more specific "direct" conic form:
//
//   min c^T x
//   s.t. A x = b, x in K,
//
//   max -b^T y
//   s.t. A^T y - z + c = 0, z in K,
//
// which corresponds to G = -I and h = 0, using a Mehrotra Predictor-Corrector 
// scheme.
//

template<typename Real>
void Mehrotra
( const Matrix<Real>& APre, 
  const Matrix<Real>& GPre,
  const Matrix<Real>& bPre, 
  const Matrix<Real>& cPre,
  const Matrix<Real>& hPre,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& labels,
        Matrix<Real>& x, 
        Matrix<Real>& y, 
        Matrix<Real>& z, 
        Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("socp::affine::Mehrotra"))    
    const bool forceSameStep = true;
    const bool vandenbergheSigma = true;

    // Equilibrate the SOCP by diagonally scaling [A;G]
    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto c = cPre;
    auto h = hPre;
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = MaxNorm(labels)+1;
    Matrix<Real> dRowA, dRowG, dCol;
    // TODO: Outer equilibration support
    Ones( dRowA, m, 1 );
    Ones( dRowG, k, 1 );
    Ones( dCol,  n, 1 );

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    // TODO: Expose this as a parameter to MehrotraCtrl
    const bool standardShift = true;
    Initialize
    ( A, G, b, c, h, orders, firstInds, labels, x, y, z, s,
      ctrl.primalInit, ctrl.dualInit, standardShift );

    Real relError = 1;
    Matrix<Real> J, d, 
                 w, wRoot, wRootInv,
                 l, lInv,
                 rmu,   rc,    rb,    rh,
                 dxAff, dyAff, dzAff, dsAff,
                 dx,    dy,    dz,    ds,
                 dzAffScaled, dsAffScaled;
    Matrix<Real> dSub;
    Matrix<Int> p;
#ifndef EL_RELEASE
    Matrix<Real> dxError, dyError, dzError, dmuError;
#endif
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Real minDist = Epsilon<Real>();
        ForceIntoSOC( s, orders, firstInds, minDist );
        ForceIntoSOC( z, orders, firstInds, minDist );
        SOCNesterovTodd( s, z, w, orders, firstInds ); 

        // Check for convergence
        // =====================
        // TODO: Adjust convergence criteria
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
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
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
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Reached maximum number of iterations, ",ctrl.maxIts,
             ", with rel. error ",relError," which does not meet the minimum ",
             "tolerance of ",ctrl.minTol);

        // Compute the affine search direction
        // ===================================
        const Real wMaxNorm = MaxNorm(w);
        const Real wMaxNormLimit = Real(100)/relError;
        if( wMaxNorm > wMaxNormLimit )
        {
            ForcePairIntoSOC( s, z, w, orders, firstInds, wMaxNormLimit );
            SOCNesterovTodd( s, z, w, orders, firstInds );
        }
        SOCSquareRoot( w, wRoot, orders, firstInds );
        SOCInverse( wRoot, wRootInv, orders, firstInds );
        SOCApplyQuadratic( wRoot, z, l, orders, firstInds );
        SOCInverse( l, lInv, orders, firstInds );
        const Real mu = Dot(s,z) / degree;

        // r_mu := l
        // ---------
        rmu = l;

        // Construct the KKT system
        // ------------------------
        KKT( A, G, w, orders, firstInds, labels, J );
        KKTRHS( rc, rb, rh, rmu, wRoot, orders, firstInds, labels, d );

        // Solve for the direction
        // -----------------------
        try 
        { 
            LDL( J, dSub, p, false );
            ldl::SolveAfter( J, dSub, p, d, false ); 
        }
        catch(...)
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Solve failed with rel. error ",relError,
                 " which does not meet the minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution
        ( m, n, d, rmu, wRoot, orders, firstInds, labels, 
          dxAff, dyAff, dzAff, dsAff );
        SOCApplyQuadratic( wRoot, dzAff, dzAffScaled, orders, firstInds );
        SOCApplyQuadratic( wRootInv, dsAff, dsAffScaled, orders, firstInds );

#ifndef EL_RELEASE
        // Sanity checks
        // -------------
        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        dmuError = dzAffScaled;
        Axpy( Real(1), dsAffScaled, dmuError );
        Axpy( Real(1), l, dmuError );
        const Real rmuNrm2 = Nrm2( rmu );
        const Real dmuErrorNrm2 = Nrm2( dmuError );
        // TODO: Also compute and print the residuals with regularization

        if( ctrl.print )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_h ||_2) = "
                 << dzErrorNrm2/(1+rhNrm2) << "\n"
                 << "  || dmuError ||_2 / (1 + || r_mu ||_2) = "
                 << dmuErrorNrm2/(1+rmuNrm2) << endl;
#endif

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = 
          MaxStepInSOC( s, dsAff, orders, firstInds, Real(1) );
        Real alphaAffDual = 
          MaxStepInSOC( z, dzAff, orders, firstInds, Real(1) );
        if( forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            cout << "  alphaAffPri = " << alphaAffPri
                 << ", alphaAffDual = " << alphaAffDual << endl;
        Real sigma;
        if( vandenbergheSigma )
        {
            sigma = Pow(Real(1)-Min(alphaAffPri,alphaAffDual),Real(3));
        }
        else
        {
            // NOTE: dz and ds are used as temporaries
            ds = s;
            dz = z;
            Axpy( alphaAffPri,  dsAff, ds );
            Axpy( alphaAffDual, dzAff, dz );
            const Real muAff = Dot(ds,dz) / degree;
            // TODO: Allow the user to override this function
            sigma = Pow(muAff/mu,Real(3));
            sigma = Min(sigma,Real(1));
            if( ctrl.print )
                cout << "  muAff = " << muAff << ", mu = " << mu << endl;
        }
        if( ctrl.print )
            cout << "  sigma = " << sigma << endl;

        // Solve for the combined direction
        // ================================
        Scale( Real(1)-sigma, rc );
        Scale( Real(1)-sigma, rb );
        Scale( Real(1)-sigma, rh );
        // r_mu := l + inv(l) o ((inv(W)^T dsAff) o (W dzAff) - sigma*mu)
        // --------------------------------------------------------------
        SOCApply( dsAffScaled, dzAffScaled, rmu, orders, firstInds );
        SOCShift( rmu, -sigma*mu, orders, firstInds );
        SOCApply( lInv, rmu, orders, firstInds );
        Axpy( Real(1), l, rmu );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS( rc, rb, rh, rmu, wRoot, orders, firstInds, labels, d );
        try { ldl::SolveAfter( J, dSub, p, d, false ); }
        catch(...)
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Solve failed with rel. error ",relError,
                 " which does not meet the minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution
        ( m, n, d, rmu, wRoot, orders, firstInds, labels, dx, dy, dz, ds );
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = 
          MaxStepInSOC( s, ds, orders, firstInds, 1/ctrl.maxStepRatio );
        Real alphaDual = 
          MaxStepInSOC( z, dz, orders, firstInds, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            cout << "  alphaPri = " << alphaPri
                 << ", alphaDual = " << alphaDual << endl;
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Zero step length computed before reaching minimum tolerance "
                 "of ",ctrl.minTol);
        }
    }
}

template<typename Real>
void Mehrotra
( const AbstractDistMatrix<Real>& APre, 
  const AbstractDistMatrix<Real>& GPre,
  const AbstractDistMatrix<Real>& bPre, 
  const AbstractDistMatrix<Real>& cPre,
  const AbstractDistMatrix<Real>& hPre,
  const AbstractDistMatrix<Int>& ordersPre,
  const AbstractDistMatrix<Int>& firstIndsPre,
  const AbstractDistMatrix<Int>& labelsPre,
        AbstractDistMatrix<Real>& xPre, 
        AbstractDistMatrix<Real>& yPre, 
        AbstractDistMatrix<Real>& zPre, 
        AbstractDistMatrix<Real>& sPre,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("socp::affine::Mehrotra"))    
    const Grid& grid = APre.Grid();
    const int commRank = grid.Rank();
    const bool onlyLower = true;

    // TODO: Expose thesa as tuning parameters
    const Int cutoffPar = 1000;
    const bool forceSameStep = true;
    const bool vandenbergheSigma = true;

    // Ensure that the inputs have the appropriate read/write properties
    DistMatrix<Real> A(grid), G(grid), b(grid), c(grid), h(grid);
    A.Align(0,0);
    G.Align(0,0);
    b.Align(0,0);
    c.Align(0,0);
    A = APre;
    G = GPre;
    b = bPre;
    c = cPre;
    h = hPre;
    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    // NOTE: {x,s} do not need to be a read proxy when !ctrl.primalInit
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,control); auto& x = *xPtr;
    auto sPtr = ReadWriteProxy<Real,MC,MR>(&sPre,control); auto& s = *sPtr;
    // NOTE: {y,z} do not need to be read proxies when !ctrl.dualInit
    auto yPtr = ReadWriteProxy<Real,MC,MR>(&yPre,control); auto& y = *yPtr;
    auto zPtr = ReadWriteProxy<Real,MC,MR>(&zPre,control); auto& z = *zPtr;

    auto ordersPtr = ReadProxy<Int,VC,STAR>(&ordersPre);
    auto firstIndsPtr = ReadProxy<Int,VC,STAR>(&firstIndsPre);
    auto labelsPtr = ReadProxy<Int,VC,STAR>(&labelsPre);
    auto& orders = *ordersPtr;
    auto& firstInds = *firstIndsPtr;
    auto& labels = *labelsPtr;

    // Equilibrate the SOCP by diagonally scaling [A;G]
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = MaxNorm(labels)+1;
    DistMatrix<Real,MC,STAR> dRowA(grid),
                             dRowG(grid);
    DistMatrix<Real,MR,STAR> dCol(grid);
    // TODO: Outer equilibration support
    Ones( dRowA, m, 1 );
    Ones( dRowG, k, 1 );
    Ones( dCol,  n, 1 );

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    // TODO: Expose this as a parameter to MehrotraCtrl
    const bool standardShift = true;
    Initialize
    ( A, G, b, c, h, orders, firstInds, labels, x, y, z, s,
      ctrl.primalInit, ctrl.dualInit, standardShift, cutoffPar );

    Real relError = 1;
    DistMatrix<Real> J(grid),     d(grid),     
                     w(grid),     wRoot(grid), wRootInv(grid),
                     l(grid),     lInv(grid),
                     rc(grid),    rb(grid),    rh(grid),    rmu(grid),
                     dxAff(grid), dyAff(grid), dzAff(grid), dsAff(grid),
                     dx(grid),    dy(grid),    dz(grid),    ds(grid),
                     dsAffScaled(grid), dzAffScaled(grid);
    w.AlignWith( s );
    wRoot.AlignWith( s );
    wRootInv.AlignWith( s );
    l.AlignWith( s );
    lInv.AlignWith( s );
    dsAff.AlignWith( s );
    dzAff.AlignWith( s );
    ds.AlignWith( s );
    dz.AlignWith( s );
    rmu.AlignWith( s );
    DistMatrix<Real> dSub(grid);
    DistMatrix<Int> p(grid);
#ifndef EL_RELEASE
    DistMatrix<Real> 
      dxError(grid), dyError(grid), dzError(grid), dmuError(grid);
    dzError.AlignWith( s );
#endif
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Real minDist = Epsilon<Real>();
        ForceIntoSOC( s, orders, firstInds, minDist, cutoffPar );
        ForceIntoSOC( z, orders, firstInds, minDist, cutoffPar );
        SOCNesterovTodd( s, z, w, orders, firstInds, cutoffPar );

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
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
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
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Reached maximum number of iterations, ",ctrl.maxIts,
             ", with rel. error ",relError," which does not meet the minimum ",
             "tolerance of ",ctrl.minTol);

        // Compute the affine search direction
        // ===================================
        const Real wMaxNorm = MaxNorm(w);
        const Real wMaxNormLimit = Real(100)/relError;
        if( wMaxNorm > wMaxNormLimit )
        {
            ForcePairIntoSOC
            ( s, z, w, orders, firstInds, wMaxNormLimit, cutoffPar );
            SOCNesterovTodd( s, z, w, orders, firstInds, cutoffPar );
        }
        SOCSquareRoot( w, wRoot, orders, firstInds, cutoffPar );
        SOCInverse( wRoot, wRootInv, orders, firstInds, cutoffPar );
        SOCApplyQuadratic( wRoot, z, l, orders, firstInds, cutoffPar );
        SOCInverse( l, lInv, orders, firstInds, cutoffPar );
        const Real mu = Dot(s,z) / degree;

        // r_mu := l
        // ---------
        rmu = l;

        // Construct the KKT system
        // ------------------------
        KKT( A, G, w, orders, firstInds, labels, J, onlyLower, cutoffPar );
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, orders, firstInds, labels, d, cutoffPar );

        // Solve for the direction
        // -----------------------
        try 
        { 
            LDL( J, dSub, p, false );
            ldl::SolveAfter( J, dSub, p, d, false ); 
        }
        catch(...)
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Solve failed with rel. error ",relError,
                 " which does not meet the minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution
        ( m, n, d, rmu, wRoot, orders, firstInds, labels, 
          dxAff, dyAff, dzAff, dsAff, cutoffPar );
        SOCApplyQuadratic
        ( wRoot, dzAff, dzAffScaled, orders, firstInds, cutoffPar );
        SOCApplyQuadratic
        ( wRootInv, dsAff, dsAffScaled, orders, firstInds, cutoffPar );

#ifndef EL_RELEASE
        // Sanity checks
        // -------------
        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        dmuError = dzAffScaled;
        Axpy( Real(1), dsAffScaled, dmuError );
        Axpy( Real(1), l, dmuError );
        const Real rmuNrm2 = Nrm2( rmu );
        const Real dmuErrorNrm2 = Nrm2( dmuError );
        // TODO: Also compute and print the residuals with regularization

        if( ctrl.print && commRank == 0 )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_h ||_2) = "
                 << dzErrorNrm2/(1+rhNrm2) << "\n"
                 << "  || dmuError ||_2 / (1 + || r_mu ||_2) = "
                 << dmuErrorNrm2/(1+rmuNrm2) << endl;
#endif

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = 
          MaxStepInSOC( s, dsAff, orders, firstInds, Real(1), cutoffPar );
        Real alphaAffDual = 
          MaxStepInSOC( z, dzAff, orders, firstInds, Real(1), cutoffPar );
        if( forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            cout << "  alphaAffPri = " << alphaAffPri
                 << ", alphaAffDual = " << alphaAffDual << endl;
        Real sigma;
        if( vandenbergheSigma )
        {
            sigma = Pow(Real(1)-Min(alphaAffPri,alphaAffDual),Real(3));
        }
        else
        {
            // NOTE: dz and ds are used as temporaries
            ds = s;
            dz = z;
            Axpy( alphaAffPri,  dsAff, ds );
            Axpy( alphaAffDual, dzAff, dz );
            const Real muAff = Dot(ds,dz) / degree;
            // TODO: Allow the user to override this function
            sigma = Pow(muAff/mu,Real(3));
            sigma = Min(sigma,Real(1));
            if( ctrl.print && commRank == 0 )
                cout << "  muAff = " << muAff << ", mu = " << mu << endl;
        }
        if( ctrl.print && commRank == 0 )
            cout << "  sigma = " << sigma << endl;

        // Solve for the combined direction
        // ================================
        Scale( Real(1)-sigma, rc );
        Scale( Real(1)-sigma, rb );
        Scale( Real(1)-sigma, rh );
        // r_mu := l + inv(l) o ((inv(W)^T dsAff) o (W dzAff) - sigma*mu)
        // --------------------------------------------------------------
        SOCApply( dsAffScaled, dzAffScaled, rmu, orders, firstInds, cutoffPar );
        SOCShift( rmu, -sigma*mu, orders, firstInds );
        SOCApply( lInv, rmu, orders, firstInds, cutoffPar );
        Axpy( Real(1), l, rmu );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, orders, firstInds, labels, d, cutoffPar );
        try { ldl::SolveAfter( J, dSub, p, d, false ); }
        catch(...)
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Solve failed with rel. error ",relError,
                 " which does not meet the minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution
        ( m, n, d, rmu, wRoot, orders, firstInds, labels, dx, dy, dz, ds,
          cutoffPar );
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = 
          MaxStepInSOC
          ( s, ds, orders, firstInds, 1/ctrl.maxStepRatio, cutoffPar );
        Real alphaDual = 
          MaxStepInSOC
          ( z, dz, orders, firstInds, 1/ctrl.maxStepRatio, cutoffPar );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            cout << "  alphaPri = " << alphaPri
                 << ", alphaDual = " << alphaDual << endl;
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Zero step length computed before reaching minimum tolerance "
                 "of ",ctrl.minTol);
        }
    }
}

template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& APre, 
  const SparseMatrix<Real>& GPre,
  const Matrix<Real>& bPre, 
  const Matrix<Real>& cPre,
  const Matrix<Real>& hPre,
  const Matrix<Int>& orders,
  const Matrix<Int>& firstInds,
  const Matrix<Int>& labels,
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("socp::affine::Mehrotra"))    

    const bool forceSameStep = true;
    const bool vandenbergheSigma = true;

    // Equilibrate the SOCP by diagonally scaling [A;G]
    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto c = cPre;
    auto h = hPre;
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = MaxNorm(labels)+1;
    Matrix<Real> dRowA, dRowG, dCol;
    // TODO: Add outer equilibration support
    Ones( dRowA, m, 1 );
    Ones( dRowG, k, 1 );
    Ones( dCol,  n, 1 );

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    vector<Int> map, invMap;
    ldl::NodeInfo info;
    ldl::Separator rootSep;
    // TODO: Expose this as a parameter to MehrotraCtrl
    const bool standardShift = true;
    Initialize
    ( A, G, b, c, h, orders, firstInds, labels, x, y, z, s,
      map, invMap, rootSep, info, 
      ctrl.primalInit, ctrl.dualInit, standardShift, ctrl.qsdCtrl );

    SparseMatrix<Real> J, JOrig;
    ldl::Front<Real> JFront;
    Matrix<Real> d, 
                 w,     wRoot, wRootInv,
                 l,     lInv,
                 rc,    rb,    rh,    rmu,
                 dxAff, dyAff, dzAff, dsAff,
                 dx,    dy,    dz,    ds,
                 dzAffScaled, dsAffScaled;

    Matrix<Real> regPerm, regTmp;
    regPerm.Resize( n+m+k, 1 );
    regTmp.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )
        {
            // TODO: Generalize regularization control structure
            regPerm.Set( i, 0, ctrl.targetTol/10 );
            regTmp.Set( i, 0, ctrl.qsdCtrl.regPrimal );
        }
        else
        {
            regPerm.Set( i, 0, -ctrl.targetTol/10 );
            regTmp.Set( i, 0, -ctrl.qsdCtrl.regDual );
        }
    }

    Real relError = 1;
    Matrix<Real> dInner;
#ifndef EL_RELEASE
    Matrix<Real> dxError, dyError, dzError, dmuError;
#endif
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Real minDist = Epsilon<Real>();
        ForceIntoSOC( s, orders, firstInds, minDist );
        ForceIntoSOC( z, orders, firstInds, minDist );
        SOCNesterovTodd( s, z, w, orders, firstInds );

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
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
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
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Reached maximum number of iterations, ",ctrl.maxIts,
             ", with rel. error ",relError," which does not meet the minimum ",
             "tolerance of ",ctrl.minTol);
        const Real wMaxLimit = Pow(Epsilon<Real>(),Real(0.4));
        const Real wMaxNorm = MaxNorm(w);
        if( wMaxNorm >= wMaxLimit && relError <= ctrl.minTol )
            break;

        // Compute the affine search direction
        // ===================================
        const Real wMaxNormLimit = Real(100)/relError;
        if( wMaxNorm > wMaxNormLimit )
        {
            ForcePairIntoSOC( s, z, w, orders, firstInds, wMaxNormLimit );
            SOCNesterovTodd( s, z, w, orders, firstInds );
        }
        SOCSquareRoot( w, wRoot, orders, firstInds );
        SOCInverse( wRoot, wRootInv, orders, firstInds );
        SOCApplyQuadratic( wRoot, z, l, orders, firstInds );
        SOCInverse( l, lInv, orders, firstInds );
        const Real mu = Dot(s,z) / degree;

        // r_mu := l
        // ---------
        rmu = l;

        // Form the KKT system
        // -------------------
        KKT( A, G, w, orders, firstInds, labels, JOrig, false );
        UpdateRealPartOfDiagonal( JOrig, Real(1), regPerm );
        KKTRHS( rc, rb, rh, rmu, wRoot, orders, firstInds, labels, d );

        // Solve for the direction
        // -----------------------
        try 
        {
            J = JOrig;
            const bool innerGeomEquil = false;
            SymmetricEquil
            ( J, dInner,
              innerGeomEquil, ctrl.innerEquil, 
              ctrl.scaleTwoNorm, ctrl.basisSize, ctrl.print );
            UpdateRealPartOfDiagonal( J, Real(1), regTmp );
            if( ctrl.primalInit && ctrl.dualInit && numIts == 0 )
            {
                NestedDissection( J.LockedGraph(), map, rootSep, info );
                InvertMap( map, invMap );
            }
            JFront.Pull( J, map, info );

            LDL( info, JFront, LDL_2D );
            reg_qsd_ldl::SolveAfter
            ( JOrig, regTmp, dInner, invMap, info, JFront, d, 
              ctrl.qsdCtrl );
        } 
        catch(...)
        {
            if( relError < ctrl.minTol )
                break;
            else
                RuntimeError
                ("Solve failed with rel. error ",relError,
                 " which does not meet the minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution
        ( m, n, d, rmu, wRoot, orders, firstInds, labels, 
          dxAff, dyAff, dzAff, dsAff );
        SOCApplyQuadratic( wRoot, dzAff, dzAffScaled, orders, firstInds );
        SOCApplyQuadratic( wRootInv, dsAff, dsAffScaled, orders, firstInds );

#ifndef EL_RELEASE
        // Sanity checks
        // -------------
        dxError = rb;
        Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Multiply( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Multiply( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        dmuError = dzAffScaled;
        Axpy( Real(1), dsAffScaled, dmuError );
        Axpy( Real(1), l, dmuError );
        const Real rmuNrm2 = Nrm2( rmu );
        const Real dmuErrorNrm2 = Nrm2( dmuError );
        // TODO: Also compute and print the residuals with regularization

        if( ctrl.print )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_h ||_2) = "
                 << dzErrorNrm2/(1+rhNrm2) << "\n"
                 << "  || dmuError ||_2 / (1 + || r_mu ||_2) = "
                 << dmuErrorNrm2/(1+rmuNrm2) << endl;
#endif

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = 
          MaxStepInSOC( s, dsAff, orders, firstInds, Real(1) );
        Real alphaAffDual = 
          MaxStepInSOC( z, dzAff, orders, firstInds, Real(1) );
        if( forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            cout << "  alphaAffPri = " << alphaAffPri
                 << ", alphaAffDual = " << alphaAffDual << endl;
        Real sigma;
        if( vandenbergheSigma )
        {
            sigma = Pow(Real(1)-Min(alphaAffPri,alphaAffDual),Real(3));
        }
        else
        {
            // NOTE: dz and ds are used as temporaries
            ds = s;
            dz = z;
            Axpy( alphaAffPri,  dsAff, ds );
            Axpy( alphaAffDual, dzAff, dz );
            const Real muAff = Dot(ds,dz) / degree;
            // TODO: Allow the user to override this function
            sigma = Pow(muAff/mu,Real(3));
            sigma = Min(sigma,Real(1));
            if( ctrl.print )
                cout << "  muAff = " << muAff << ", mu = " << mu << endl;
        }
        if( ctrl.print )
            cout << "  sigma = " << sigma << endl;

        // Solve for the combined direction
        // ================================
        Scale( Real(1)-sigma, rc );
        Scale( Real(1)-sigma, rb );
        Scale( Real(1)-sigma, rh );
        // r_mu := l + inv(l) o ((inv(W)^T dsAff) o (W dzAff) - sigma*mu)
        // --------------------------------------------------------------
        SOCApplyQuadratic( wRootInv, dsAff, orders, firstInds );
        SOCApplyQuadratic( wRoot,    dzAff, orders, firstInds );
        SOCApply( dsAff, dzAff, rmu, orders, firstInds );
        SOCShift( rmu, -sigma*mu, orders, firstInds );
        SOCApply( lInv, rmu, orders, firstInds );
        Axpy( Real(1), l, rmu );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS( rc, rb, rh, rmu, wRoot, orders, firstInds, labels, d );
        try 
        {
            reg_qsd_ldl::SolveAfter
            ( JOrig, regTmp, dInner, invMap, info, JFront, d, 
              ctrl.qsdCtrl );
        } 
        catch(...)
        {
            if( relError < ctrl.minTol )
                break;
            else
                RuntimeError
                ("Solve failed with rel. error ",relError,
                 " which does not meet the minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution
        ( m, n, d, rmu, wRoot, orders, firstInds, labels, dx, dy, dz, ds );

        // Update the current estimates
        // ============================
        Real alphaPri = 
          MaxStepInSOC( s, ds, orders, firstInds, 1/ctrl.maxStepRatio );
        Real alphaDual = 
          MaxStepInSOC( z, dz, orders, firstInds, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            cout << "  alphaPri = " << alphaPri
                 << ", alphaDual = " << alphaDual << endl;
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( relError < ctrl.minTol )
                break;
            else
                RuntimeError
                ("Zero step length computed before reaching minimum tolerance "
                 "of ",ctrl.minTol);
        }
    }
}

template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& APre,
  const DistSparseMatrix<Real>& GPre,
  const DistMultiVec<Real>& bPre,
  const DistMultiVec<Real>& cPre,
  const DistMultiVec<Real>& hPre,
  const DistMultiVec<Int>& orders,
  const DistMultiVec<Int>& firstInds,
  const DistMultiVec<Int>& labels,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("socp::affine::Mehrotra"))    
    mpi::Comm comm = APre.Comm();
    const int commRank = mpi::Rank(comm);
    const bool onlyLower = false;
    Timer timer;

    // TODO: Expose as tuning parameters
    const Int cutoffPar = 1000;
    const bool forceSameStep = true;
    const bool vandenbergheSigma = true;

    // Equilibrate the SOCP by diagonally scaling [A;G]
    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto h = hPre;
    auto c = cPre;

    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = MaxNorm(labels)+1;
    DistMultiVec<Real> dRowA(comm), dRowG(comm), dCol(comm);
    // TODO: Outer equilibration support
    Ones( dRowA, m, 1 );
    Ones( dRowG, k, 1 );
    Ones( dCol,  n, 1 );

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    DistMap map, invMap;
    ldl::DistNodeInfo info;
    ldl::DistSeparator rootSep;
    // TODO: Expose this as a parameter to MehrotraCtrl
    const bool standardShift = true;
    if( commRank == 0 && ctrl.time )
        timer.Start();
    Initialize
    ( A, G, b, c, h, orders, firstInds, labels, x, y, z, s,
      map, invMap, rootSep, info, 
      ctrl.primalInit, ctrl.dualInit, standardShift, cutoffPar, 
      ctrl.qsdCtrl );
    if( commRank == 0 && ctrl.time )
        cout << "  Init: " << timer.Stop() << " secs" << endl;

    DistSparseMultMeta metaOrig, meta;
    DistSparseMatrix<Real> J(comm), JOrig(comm);
    ldl::DistFront<Real> JFront;
    DistMultiVec<Real> d(comm),
                       w(comm),     wRoot(comm), wRootInv(comm),
                       l(comm),     lInv(comm),
                       rc(comm),    rb(comm),    rh(comm),    rmu(comm),
                       dxAff(comm), dyAff(comm), dzAff(comm), dsAff(comm),
                       dx(comm),    dy(comm),    dz(comm),    ds(comm),
                       dzAffScaled(comm), dsAffScaled(comm);

    DistMultiVec<Real> regPerm(comm), regTmp(comm);
    regTmp.Resize( n+m+k, 1 );
    regPerm.Resize( n+m+k, 1 );
    for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
    {
        const Int i = regTmp.GlobalRow(iLoc);
        if( i < n )
        {
            regPerm.SetLocal( iLoc, 0, ctrl.targetTol/10 );
            regTmp.SetLocal( iLoc, 0, ctrl.qsdCtrl.regPrimal );
        }
        else
        {
            regPerm.SetLocal( iLoc, 0, -ctrl.targetTol/10 );
            regTmp.SetLocal( iLoc, 0, -ctrl.qsdCtrl.regDual );
        }
    }

    Real relError = 1;
    DistMultiVec<Real> dInner(comm);
#ifndef EL_RELEASE
    DistMultiVec<Real> dxError(comm), dyError(comm), 
                       dzError(comm), dmuError(comm);
#endif
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        // TODO: Let this be tunable
        const Real minDist = Epsilon<Real>();
        ForceIntoSOC( s, orders, firstInds, minDist, cutoffPar );
        ForceIntoSOC( z, orders, firstInds, minDist, cutoffPar );
        SOCNesterovTodd( s, z, w, orders, firstInds, cutoffPar );

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
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
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
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Reached maximum number of iterations, ",ctrl.maxIts,
             ", with rel. error ",relError," which does not meet the minimum ",
             "tolerance of ",ctrl.minTol);
        const Real wMaxLimit = Pow(Epsilon<Real>(),Real(0.4));
        const Real wMaxNorm = MaxNorm(w);
        if( wMaxNorm >= wMaxLimit && relError <= ctrl.minTol )
            break;

        // Compute the affine search direction
        // ===================================
        const Real wMaxNormLimit = Real(100)/relError;
        if( wMaxNorm > wMaxNormLimit )
        {
            ForcePairIntoSOC
            ( s, z, w, orders, firstInds, wMaxNormLimit, cutoffPar );
            SOCNesterovTodd( s, z, w, orders, firstInds, cutoffPar );
        }
        SOCSquareRoot( w, wRoot, orders, firstInds, cutoffPar );
        SOCInverse( wRoot, wRootInv, orders, firstInds, cutoffPar );
        SOCApplyQuadratic( wRoot, z, l, orders, firstInds, cutoffPar );
        SOCInverse( l, lInv, orders, firstInds, cutoffPar );
        const Real mu = Dot(s,z) / degree;

        // r_mu := l
        // ---------
        rmu = l;

        // Construct the KKT system
        // ------------------------
        KKT( A, G, w, orders, firstInds, labels, JOrig, onlyLower, cutoffPar );
        UpdateRealPartOfDiagonal( JOrig, Real(1), regPerm );
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, orders, firstInds, labels, d, cutoffPar );

        // Solve for the direction
        // -----------------------
        try
        {
            // Cache the metadata for the finalized JOrig
            if( numIts == 0 )
                metaOrig = JOrig.InitializeMultMeta();
            else
                JOrig.multMeta = metaOrig;
            J = JOrig;
            if( commRank == 0 && ctrl.time )
                timer.Start();
            const bool innerGeomEquil = false;
            SymmetricEquil
            ( J, dInner, 
              innerGeomEquil, ctrl.innerEquil, 
              ctrl.scaleTwoNorm, ctrl.basisSize, ctrl.print, ctrl.time );
            if( commRank == 0 && ctrl.time )
                cout << "  Equilibration: " << timer.Stop() << " secs" << endl;
            UpdateRealPartOfDiagonal( J, Real(1), regTmp );
            // Cache the metadata for the finalized J
            if( numIts == 0 )
            {
                meta = J.InitializeMultMeta();
                if( ctrl.primalInit && ctrl.dualInit )
                {
                    if( commRank == 0 && ctrl.time )
                        timer.Start();
                    NestedDissection( J.LockedDistGraph(), map, rootSep, info );
                    if( commRank == 0 && ctrl.time )
                        cout << "  ND: " << timer.Stop() << " secs" << endl;
                    InvertMap( map, invMap );
                }
            }
            else
                J.multMeta = meta;
            JFront.Pull( J, map, rootSep, info );

            if( commRank == 0 && ctrl.time )
                timer.Start();
            LDL( info, JFront, LDL_2D );
            if( commRank == 0 && ctrl.time )
                cout << "  LDL: " << timer.Stop() << " secs" << endl;

            if( commRank == 0 && ctrl.time )
                timer.Start();
            reg_qsd_ldl::SolveAfter
            ( JOrig, regTmp, dInner, invMap, info, JFront, d, 
              ctrl.qsdCtrl );
            if( commRank == 0 && ctrl.time )
                cout << "  Affine: " << timer.Stop() << " secs" << endl;
        }
        catch(...)
        {
            if( relError < ctrl.minTol )
                break;
            else
                RuntimeError
                ("Solve failed with rel. error ",relError,
                 " which does not meet the minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution
        ( m, n, d, rmu, wRoot, orders, firstInds, labels, 
          dxAff, dyAff, dzAff, dsAff, cutoffPar );
        SOCApplyQuadratic
        ( wRoot, dzAff, dzAffScaled, orders, firstInds, cutoffPar );
        SOCApplyQuadratic
        ( wRootInv, dsAff, dsAffScaled, orders, firstInds, cutoffPar );

#ifndef EL_RELEASE
        // Sanity checks
        // -------------
        dxError = rb;
        Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dyError = rc;
        Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Multiply( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dzError = rh;
        Multiply( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        dmuError = dzAffScaled;
        Axpy( Real(1), dsAffScaled, dmuError );
        Axpy( Real(1), l, dmuError );
        const Real rmuNrm2 = Nrm2( rmu );
        const Real dmuErrorNrm2 = Nrm2( dmuError );
        // TODO: Also compute and print the residuals with regularization

        if( ctrl.print && commRank == 0 )
            cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                 << dxErrorNrm2/(1+rbNrm2) << "\n"
                 << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                 << dyErrorNrm2/(1+rcNrm2) << "\n"
                 << "  || dzError ||_2 / (1 + || r_h ||_2) = "
                 << dzErrorNrm2/(1+rhNrm2) << "\n"
                 << "  || dmuError ||_2 / (1 + || r_mu ||_2) = " 
                 << dmuErrorNrm2/(1+rmuNrm2) << endl;
#endif

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = 
          MaxStepInSOC( s, dsAff, orders, firstInds, Real(1), cutoffPar );
        Real alphaAffDual = 
          MaxStepInSOC( z, dzAff, orders, firstInds, Real(1), cutoffPar );
        if( forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            cout << "  alphaAffPri = " << alphaAffPri
                 << ", alphaAffDual = " << alphaAffDual << endl;
        Real sigma;
        if( vandenbergheSigma )
        {
            sigma = Pow(Real(1)-Min(alphaAffPri,alphaAffDual),Real(3));
        }
        else
        {
            // NOTE: dz and ds are used as temporaries
            ds = s;
            dz = z;
            Axpy( alphaAffPri,  dsAff, ds );
            Axpy( alphaAffDual, dzAff, dz );
            const Real muAff = Dot(ds,dz) / degree;
            // TODO: Allow the user to override this function
            sigma = Pow(muAff/mu,Real(3));
            sigma = Min(sigma,Real(1));
            if( ctrl.print && commRank == 0 )
                cout << "  muAff = " << muAff << ", mu = " << mu << endl;
        }
        if( ctrl.print && commRank == 0 )
            cout << "  sigma = " << sigma << endl;

        // Solve for the combined direction
        // ================================
        Scale( Real(1)-sigma, rc );
        Scale( Real(1)-sigma, rb );
        Scale( Real(1)-sigma, rh );
        // r_mu := l + inv(l) o ((inv(W)^T dsAff) o (W dzAff) - sigma*mu)
        // --------------------------------------------------------------
        SOCApplyQuadratic( wRootInv, dsAff, orders, firstInds );
        SOCApplyQuadratic( wRoot,    dzAff, orders, firstInds );
        SOCApply( dsAff, dzAff, rmu, orders, firstInds );
        SOCShift( rmu, -sigma*mu, orders, firstInds );
        SOCApply( lInv, rmu, orders, firstInds );
        Axpy( Real(1), l, rmu );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, orders, firstInds, labels, d, cutoffPar );
        try
        {
            if( commRank == 0 && ctrl.time )
                timer.Start();
            reg_qsd_ldl::SolveAfter
            ( JOrig, regTmp, dInner, invMap, info, JFront, d, 
              ctrl.qsdCtrl );
            if( commRank == 0 && ctrl.time )
                cout << "  Corrector: " << timer.Stop() << " secs" << endl;
        }
        catch(...)
        {
            if( relError < ctrl.minTol )
                break;
            else
                RuntimeError
                ("Solve failed with rel. error ",relError,
                 " which does not meet the minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution
        ( m, n, d, rmu, wRoot, orders, firstInds, labels, dx, dy, dz, ds, 
          cutoffPar );

        // Update the current estimates
        // ============================
        Real alphaPri = 
          MaxStepInSOC
          ( s, ds, orders, firstInds, 1/ctrl.maxStepRatio, cutoffPar );
        Real alphaDual = 
          MaxStepInSOC
          ( z, dz, orders, firstInds, 1/ctrl.maxStepRatio, cutoffPar );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            cout << "  alphaPri = " << alphaPri << "\n"
                 << "  alphaDual = " << alphaDual << endl;
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( relError < ctrl.minTol )
                break;
            else
                RuntimeError
                ("Zero step length computed before reaching minimum tolerance "
                 "of ",ctrl.minTol);
        }
    }
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& labels, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, \
    const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
    const AbstractDistMatrix<Int>& orders, \
    const AbstractDistMatrix<Int>& firstInds, \
    const AbstractDistMatrix<Int>& labels, \
          AbstractDistMatrix<Real>& x, \
          AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z, \
          AbstractDistMatrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
    const Matrix<Int>& labels, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
    const DistMultiVec<Int>& orders, \
    const DistMultiVec<Int>& firstInds, \
    const DistMultiVec<Int>& labels, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMultiVec<Real>& s, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace affine
} // namespace socp
} // namespace El
