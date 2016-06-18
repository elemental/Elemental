/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include <El.hpp>
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
        Matrix<Real>& x, 
        Matrix<Real>& y, 
        Matrix<Real>& z, 
        Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Real eps = limits::Epsilon<Real>();

    // TODO: Move these into the control structure
    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const bool standardShift = true;

    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto c = cPre;
    auto h = hPre;
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = soc::Degree( firstInds );
    Matrix<Real> dRowA, dRowG, dCol;
    if( ctrl.outerEquil )
    {
        cone::GeomEquil
        ( A, G, dRowA, dRowG, dCol, orders, firstInds, ctrl.print );
        DiagonalSolve( LEFT, NORMAL, dRowA, b );
        DiagonalSolve( LEFT, NORMAL, dRowG, h );
        DiagonalSolve( LEFT, NORMAL, dCol,  c );
        if( ctrl.primalInit )
        {
            DiagonalScale( LEFT, NORMAL, dCol,  x );
            DiagonalSolve( LEFT, NORMAL, dRowG, s );
        }
        if( ctrl.dualInit )
        {
            DiagonalScale( LEFT, NORMAL, dRowA, y );
            DiagonalScale( LEFT, NORMAL, dRowG, z );
        }
    }
    else
    {
        Ones( dRowA, m, 1 );
        Ones( dRowG, k, 1 );
        Ones( dCol,  n, 1 );
    }

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );
    if( ctrl.print )
    {
        const Real ANrm1 = OneNorm( A );
        const Real GNrm1 = OneNorm( G );
        Output("|| A ||_1 = ",ANrm1);
        Output("|| G ||_1 = ",GNrm1);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
        Output("|| h ||_2 = ",hNrm2);
    }

    Initialize
    ( A, G, b, c, h, orders, firstInds, x, y, z, s,
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
    Permutation p;
    Matrix<Real> dxError, dyError, dzError, dmuError;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Real minDist = eps;
        soc::PushInto( s, orders, firstInds, minDist );
        soc::PushInto( z, orders, firstInds, minDist );
        soc::NesterovTodd( s, z, w, orders, firstInds ); 

        // Check for convergence
        // =====================
        // TODO: Adjust convergence criteria
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        rb *= -1;
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h;
        rh *= -1;
        Gemv( NORMAL, Real(1), G, x, Real(1), rh );
        rh += s;
        const Real rhNrm2 = Nrm2( rh );
        const Real rhConv = rhNrm2 / (1+hNrm2);

        // Now check the pieces
        // --------------------
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( x );
            const Real yNrm2 = Nrm2( y );
            const Real zNrm2 = Nrm2( z );
            const Real sNrm2 = Nrm2( s );
            Output
            ("iter ",numIts,":\n",Indent(),
             "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
             "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
             "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
             "  ||  s  ||_2 = ",sNrm2,"\n",Indent(),
             "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
             "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
             "  || r_h ||_2 = ",rhNrm2,"\n",Indent(),
             "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
             "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
             "  || r_h ||_2 / (1 + || h ||_2) = ",rhConv,"\n",Indent(),
             "  primal = ",primObj,"\n",Indent(),
             "  dual   = ",dualObj,"\n",Indent(),
             "  |primal - dual| / (1 + |primal|) = ",objConv);
        }
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Reached maximum number of iterations, ",ctrl.maxIts,
             ", with rel. error ",relError," which does not meet the minimum ",
             "tolerance of ",ctrl.minTol);

        // Compute the affine search direction
        // ===================================
        Real wMaxNorm = MaxNorm(w);
        const Real wMaxNormLimit = 
          Max(ctrl.wSafeMaxNorm,10/Min(Real(1),relError));
        if( wMaxNorm > wMaxNormLimit )
        {
            soc::PushPairInto( s, z, w, orders, firstInds, wMaxNormLimit );
            soc::NesterovTodd( s, z, w, orders, firstInds );
            wMaxNorm = MaxNorm(w);
        }
        soc::SquareRoot( w, wRoot, orders, firstInds );
        soc::Inverse( wRoot, wRootInv, orders, firstInds );
        soc::ApplyQuadratic( wRoot, z, l, orders, firstInds );
        soc::Inverse( l, lInv, orders, firstInds );
        const Real mu = Dot(s,z) / degree;

        // r_mu := l
        // ---------
        rmu = l;

        // Construct the KKT system
        // ------------------------
        KKT( A, G, w, orders, firstInds, J );
        KKTRHS( rc, rb, rh, rmu, wRoot, orders, firstInds, d );

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
        ( m, n, d, rmu, wRoot, orders, firstInds,
          dxAff, dyAff, dzAff, dsAff );
        soc::ApplyQuadratic( wRoot, dzAff, dzAffScaled, orders, firstInds );
        soc::ApplyQuadratic( wRootInv, dsAff, dsAffScaled, orders, firstInds );

        if( ctrl.checkResiduals && ctrl.print )
        {
            dxError = rb;
            Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
            const Real dxErrorNrm2 = Nrm2( dxError );

            dyError = rc;
            Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
            Gemv( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
            const Real dyErrorNrm2 = Nrm2( dyError );

            dzError = rh;
            Gemv( NORMAL, Real(1), G, dxAff, Real(1), dzError );
            dzError += dsAff;
            const Real dzErrorNrm2 = Nrm2( dzError );

            dmuError = dzAffScaled;
            dmuError += dsAffScaled;
            dmuError += l;
            const Real rmuNrm2 = Nrm2( rmu );
            const Real dmuErrorNrm2 = Nrm2( dmuError );

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+rhNrm2),"\n",Indent(),
             "|| dmuError ||_2 / (1 + || r_mu ||_2) = ",
             dmuErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = soc::MaxStep(s,dsAff,orders,firstInds,Real(1));
        Real alphaAffDual = soc::MaxStep(z,dzAff,orders,firstInds,Real(1));
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and ds are used as temporaries
        ds = s;
        dz = z;
        Axpy( alphaAffPri,  dsAff, ds );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(ds,dz) / degree;
        if( ctrl.print )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma = centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        rc *= 1-sigma;
        rb *= 1-sigma;
        rh *= 1-sigma;
        if( ctrl.mehrotra )
        {
            // r_mu := l + inv(l) o ((inv(W)^T dsAff) o (W dzAff) - sigma*mu)
            // --------------------------------------------------------------
            soc::Apply( dsAffScaled, dzAffScaled, rmu, orders, firstInds );
            soc::Shift( rmu, -sigma*mu, orders, firstInds );
            soc::Apply( lInv, rmu, orders, firstInds );
            rmu += l;
        }
        else
        {
            // r_mu -= sigma*mu*inv(l)
            // -----------------------
            Axpy( -sigma*mu, lInv, rmu );
        }

        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS( rc, rb, rh, rmu, wRoot, orders, firstInds, d );
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
        ( m, n, d, rmu, wRoot, orders, firstInds, dx, dy, dz, ds );

        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = 
          soc::MaxStep( s, ds, orders, firstInds, 1/ctrl.maxStepRatio );
        Real alphaDual = 
          soc::MaxStep( z, dz, orders, firstInds, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
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
    SetIndent( indent );
    if( ctrl.outerEquil )
    {
        DiagonalSolve( LEFT, NORMAL, dCol,  x );
        DiagonalSolve( LEFT, NORMAL, dRowA, y );
        DiagonalSolve( LEFT, NORMAL, dRowG, z );
        DiagonalScale( LEFT, NORMAL, dRowG, s );
    }
}

template<typename Real>
void Mehrotra
( const ElementalMatrix<Real>& APre, 
  const ElementalMatrix<Real>& GPre,
  const ElementalMatrix<Real>& bPre, 
  const ElementalMatrix<Real>& cPre,
  const ElementalMatrix<Real>& hPre,
  const ElementalMatrix<Int>& ordersPre,
  const ElementalMatrix<Int>& firstIndsPre,
        ElementalMatrix<Real>& xPre, 
        ElementalMatrix<Real>& yPre, 
        ElementalMatrix<Real>& zPre, 
        ElementalMatrix<Real>& sPre,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Real eps = limits::Epsilon<Real>();
    const bool onlyLower = true;

    // TODO: Move these into the control structure
    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const Int cutoffPar = 1000;
    const bool standardShift = true;

    const Grid& grid = APre.Grid();
    const int commRank = grid.Rank();

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

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    DistMatrixReadWriteProxy<Real,Real,MC,MR>
    // NOTE: {x,s} do not need to be a read proxy when !ctrl.primalInit
      xProx( xPre, control ),
      sProx( sPre, control ),
    // NOTE: {y,z} do not need to be read proxies when !ctrl.dualInit
      yProx( yPre, control ),
      zProx( zPre, control );
    auto& x = xProx.Get();
    auto& s = sProx.Get();
    auto& y = yProx.Get();
    auto& z = zProx.Get();

    DistMatrixReadProxy<Int,Int,VC,STAR>
      ordersProx( ordersPre ),
      firstIndsProx( firstIndsPre );
    auto& orders = ordersProx.GetLocked();
    auto& firstInds = firstIndsProx.GetLocked();

    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = soc::Degree( firstInds );
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
    if( ctrl.print )
    {
        const Real ANrm1 = OneNorm( A );
        const Real GNrm1 = OneNorm( G );
        if( commRank == 0 )
        {
            Output("|| A ||_1 = ",ANrm1);
            Output("|| G ||_1 = ",GNrm1);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
            Output("|| h ||_2 = ",hNrm2);
        }
    }

    Initialize
    ( A, G, b, c, h, orders, firstInds, x, y, z, s,
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
    DistPermutation p(grid);
    DistMatrix<Real> 
      dxError(grid), dyError(grid), dzError(grid), dmuError(grid);
    dzError.AlignWith( s );
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Real minDist = eps;
        soc::PushInto( s, orders, firstInds, minDist, cutoffPar );
        soc::PushInto( z, orders, firstInds, minDist, cutoffPar );
        soc::NesterovTodd( s, z, w, orders, firstInds, cutoffPar );

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        rb *= -1;
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h;
        rh *= -1;
        Gemv( NORMAL, Real(1), G, x, Real(1), rh );
        rh += s;
        const Real rhNrm2 = Nrm2( rh );
        const Real rhConv = rhNrm2 / (1+hNrm2);

        // Now check the pieces
        // --------------------
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( x );
            const Real yNrm2 = Nrm2( y );
            const Real zNrm2 = Nrm2( z );
            const Real sNrm2 = Nrm2( s );
            if( commRank == 0 )
                Output
                ("iter ",numIts,":\n",Indent(),
                 "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
                 "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
                 "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
                 "  ||  s  ||_2 = ",sNrm2,"\n",Indent(),
                 "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
                 "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
                 "  || r_h ||_2 = ",rhNrm2,"\n",Indent(),
                 "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
                 "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
                 "  || r_h ||_2 / (1 + || h ||_2) = ",rhConv,"\n",Indent(),
                 "  primal = ",primObj,"\n",Indent(),
                 "  dual   = ",dualObj,"\n",Indent(),
                 "  |primal - dual| / (1 + |primal|) = ",objConv);
        }
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Reached maximum number of iterations, ",ctrl.maxIts,
             ", with rel. error ",relError," which does not meet the minimum ",
             "tolerance of ",ctrl.minTol);

        // Compute the affine search direction
        // ===================================
        Real wMaxNorm = MaxNorm(w);
        const Real wMaxNormLimit =
          Max(ctrl.wSafeMaxNorm,10/Min(Real(1),relError));
        if( wMaxNorm > wMaxNormLimit )
        {
            soc::PushPairInto
            ( s, z, w, orders, firstInds, wMaxNormLimit, cutoffPar );
            soc::NesterovTodd( s, z, w, orders, firstInds, cutoffPar );
            wMaxNorm = MaxNorm(w);
        }
        soc::SquareRoot( w, wRoot, orders, firstInds, cutoffPar );
        soc::Inverse( wRoot, wRootInv, orders, firstInds, cutoffPar );
        soc::ApplyQuadratic( wRoot, z, l, orders, firstInds, cutoffPar );
        soc::Inverse( l, lInv, orders, firstInds, cutoffPar );
        const Real mu = Dot(s,z) / degree;

        // r_mu := l
        // ---------
        rmu = l;

        // Construct the KKT system
        // ------------------------
        KKT( A, G, w, orders, firstInds, J, onlyLower, cutoffPar );
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, orders, firstInds, d, cutoffPar );

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
        ( m, n, d, rmu, wRoot, orders, firstInds,
          dxAff, dyAff, dzAff, dsAff, cutoffPar );
        soc::ApplyQuadratic
        ( wRoot, dzAff, dzAffScaled, orders, firstInds, cutoffPar );
        soc::ApplyQuadratic
        ( wRootInv, dsAff, dsAffScaled, orders, firstInds, cutoffPar );

        if( ctrl.checkResiduals && ctrl.print )
        {
            dxError = rb;
            Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
            const Real dxErrorNrm2 = Nrm2( dxError );

            dyError = rc;
            Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
            Gemv( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
            const Real dyErrorNrm2 = Nrm2( dyError );

            dzError = rh;
            Gemv( NORMAL, Real(1), G, dxAff, Real(1), dzError );
            dzError += dsAff;
            const Real dzErrorNrm2 = Nrm2( dzError );

            dmuError = dzAffScaled;
            dmuError += dsAffScaled;
            dmuError += l;
            const Real rmuNrm2 = Nrm2( rmu );
            const Real dmuErrorNrm2 = Nrm2( dmuError );

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rhNrm2),"\n",Indent(),
                 "|| dmuError ||_2 / (1 + || r_mu ||_2) = ",
                 dmuErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = 
          soc::MaxStep( s, dsAff, orders, firstInds, Real(1), cutoffPar );
        Real alphaAffDual = 
          soc::MaxStep( z, dzAff, orders, firstInds, Real(1), cutoffPar );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and ds are used as temporaries
        ds = s;
        dz = z;
        Axpy( alphaAffPri,  dsAff, ds );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(ds,dz) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma = centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        rc *= 1-sigma;
        rb *= 1-sigma;
        rh *= 1-sigma;
        if( ctrl.mehrotra )
        {
            // r_mu := l + inv(l) o ((inv(W)^T dsAff) o (W dzAff) - sigma*mu)
            // --------------------------------------------------------------
            soc::Apply
            ( dsAffScaled, dzAffScaled, rmu, orders, firstInds, cutoffPar );
            soc::Shift( rmu, -sigma*mu, orders, firstInds );
            soc::Apply( lInv, rmu, orders, firstInds, cutoffPar );
            rmu += l;
        }
        else
        {
            // r_mu -= sigma*mu*inv(l)
            // -----------------------
            Axpy( -sigma*mu, lInv, rmu );
        }

        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, orders, firstInds, d, cutoffPar );
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
        ( m, n, d, rmu, wRoot, orders, firstInds, dx, dy, dz, ds,
          cutoffPar );
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = 
          soc::MaxStep
          ( s, ds, orders, firstInds, 1/ctrl.maxStepRatio, cutoffPar );
        Real alphaDual = 
          soc::MaxStep
          ( z, dz, orders, firstInds, 1/ctrl.maxStepRatio, cutoffPar );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
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
    SetIndent( indent );
    if( ctrl.outerEquil )
    {
        DiagonalSolve( LEFT, NORMAL, dCol,  x );
        DiagonalSolve( LEFT, NORMAL, dRowA, y );
        DiagonalSolve( LEFT, NORMAL, dRowG, z );
        DiagonalScale( LEFT, NORMAL, dRowG, s );
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
        Matrix<Real>& x,
        Matrix<Real>& y,
        Matrix<Real>& z,
        Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Real eps = limits::Epsilon<Real>();
    const bool onlyLower = false;

    // TODO: Move these into the control structure
    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const bool cutoffSparse = 64;
    const bool standardShift = true;
    const Real gamma = Pow(eps,Real(0.35));
    const Real delta = Pow(eps,Real(0.35));
    const Real beta =  Pow(eps,Real(0.35));
    const Real gammaTmp = Pow(eps,Real(0.25));
    const Real deltaTmp = Pow(eps,Real(0.25));
    const Real betaTmp  = Pow(eps,Real(0.25));

    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto c = cPre;
    auto h = hPre;
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = soc::Degree( firstInds );
    Matrix<Real> dRowA, dRowG, dCol;
    if( ctrl.outerEquil )
    {
        cone::RuizEquil
        ( A, G, dRowA, dRowG, dCol, orders, firstInds, ctrl.print );

        DiagonalSolve( LEFT, NORMAL, dRowA, b );
        DiagonalSolve( LEFT, NORMAL, dRowG, h );
        DiagonalSolve( LEFT, NORMAL, dCol,  c );
        if( ctrl.primalInit )
        {
            DiagonalScale( LEFT, NORMAL, dCol,  x );
            DiagonalSolve( LEFT, NORMAL, dRowG, s );
        }
        if( ctrl.dualInit )
        {
            DiagonalScale( LEFT, NORMAL, dRowA, y );
            DiagonalScale( LEFT, NORMAL, dRowG, z );
        }
    }
    else
    {
        Ones( dRowA, m, 1 );
        Ones( dRowG, k, 1 );
        Ones( dCol,  n, 1 );
    }

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );
    const Real twoNormEstA = TwoNormEstimate( A, ctrl.basisSize );
    const Real twoNormEstG = TwoNormEstimate( G, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstA + twoNormEstG + 1;
    if( ctrl.print )
    {
        Output("|| A ||_2 estimate: ",twoNormEstA);
        Output("|| G ||_2 estimate: ",twoNormEstG);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
        Output("|| h ||_2 = ",hNrm2);
    }

    Initialize
    ( A, G, b, c, h, orders, firstInds, x, y, z, s,
      ctrl.primalInit, ctrl.dualInit, standardShift, ctrl.solveCtrl );

    // Form the offsets for the sparse embedding of the barrier's Hessian
    // ==================================================================
    Matrix<Int> 
      origToSparseOrders,    sparseToOrigOrders,
      origToSparseFirstInds, sparseToOrigFirstInds,
      sparseOrders,          sparseFirstInds;
    soc::EmbeddingMaps
    ( orders, firstInds,
      sparseOrders, sparseFirstInds,
      origToSparseOrders, origToSparseFirstInds,
      sparseToOrigOrders, sparseToOrigFirstInds, cutoffSparse );
    const Int kSparse = sparseOrders.Height();

    SparseMatrix<Real> J, JOrig;
    ldl::Front<Real> JFront;
    Matrix<Real> d, 
                 w,     wRoot, wRootInv,
                 l,     lInv,
                 rc,    rb,    rh,    rmu,
                 dxAff, dyAff, dzAff, dsAff,
                 dx,    dy,    dz,    ds,
                 dzAffScaled, dsAffScaled;

    // TODO: Expose regularization rules to user
    Matrix<Real> regTmp;
    regTmp.Resize( n+m+kSparse, 1 );
    for( Int i=0; i<n+m+kSparse; ++i )
    {
        if( i < n )
        {
            regTmp(i) = gammaTmp*gammaTmp;
        }
        else if( i < n+m )
        {
            regTmp(i) = -deltaTmp*deltaTmp;
        }
        else
        {
            const Int iCone = i-(n+m);
            const Int firstInd = sparseFirstInds(iCone);
            const Int order = sparseToOrigOrders(iCone);
            const Int sparseOrder = sparseOrders(iCone);
            const bool embedded = ( order != sparseOrder );
          
            // TODO: Use different diagonal modification for the auxiliary
            //       variables? These diagonal entries are always +-1.
            if( embedded && iCone == firstInd+sparseOrder-1 )
                regTmp(i) = betaTmp*betaTmp;
            else 
                regTmp(i) = -betaTmp*betaTmp;
        }
    }
    regTmp *= origTwoNormEst;

    // Form the static portion of the KKT system
    // =========================================
    SparseMatrix<Real> JStatic;
    StaticKKT
    ( A, G, gamma, delta, beta, 
      orders, firstInds, origToSparseOrders, origToSparseFirstInds, 
      kSparse, JStatic, onlyLower );

    vector<Int> map, invMap;
    ldl::NodeInfo info;
    ldl::Separator rootSep;
    NestedDissection( JStatic.LockedGraph(), map, rootSep, info );
    InvertMap( map, invMap );
 
    Real relError = 1;
    Matrix<Real> dInner;
    Matrix<Real> dxError, dyError, dzError, dmuError;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Real minDist = eps;
        soc::PushInto( s, orders, firstInds, minDist );
        soc::PushInto( z, orders, firstInds, minDist );
        soc::NesterovTodd( s, z, w, orders, firstInds );

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        rb *= -1;
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h;
        rh *= -1;
        Multiply( NORMAL, Real(1), G, x, Real(1), rh );
        rh += s;
        const Real rhNrm2 = Nrm2( rh );
        const Real rhConv = rhNrm2 / (1+hNrm2);

        // Now check the pieces
        // --------------------
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( x );
            const Real yNrm2 = Nrm2( y );
            const Real zNrm2 = Nrm2( z );
            const Real sNrm2 = Nrm2( s );
            Output
            ("iter ",numIts,":\n",Indent(),
             "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
             "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
             "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
             "  ||  s  ||_2 = ",sNrm2,"\n",Indent(),
             "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
             "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
             "  || r_h ||_2 = ",rhNrm2,"\n",Indent(),
             "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
             "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
             "  || r_h ||_2 / (1 + || h ||_2) = ",rhConv,"\n",Indent(),
             "  primal = ",primObj,"\n",Indent(),
             "  dual   = ",dualObj,"\n",Indent(),
             "  |primal - dual| / (1 + |primal|) = ",objConv);
        }
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Reached maximum number of iterations, ",ctrl.maxIts,
             ", with rel. error ",relError," which does not meet the minimum ",
             "tolerance of ",ctrl.minTol);
        Real wMaxNorm = MaxNorm(w);
        if( wMaxNorm >= ctrl.wMaxLimit && relError <= ctrl.minTol )
            break;

        // Compute the affine search direction
        // ===================================
        const Real wMaxNormLimit =
          Max(ctrl.wSafeMaxNorm,10/Min(Real(1),relError));
        if( wMaxNorm > wMaxNormLimit )
        {
            soc::PushPairInto( s, z, w, orders, firstInds, wMaxNormLimit );
            soc::NesterovTodd( s, z, w, orders, firstInds );
            wMaxNorm = MaxNorm(w);
        }
        soc::SquareRoot( w, wRoot, orders, firstInds );
        soc::Inverse( wRoot, wRootInv, orders, firstInds );
        soc::ApplyQuadratic( wRoot, z, l, orders, firstInds );
        soc::Inverse( l, lInv, orders, firstInds );
        const Real mu = Dot(s,z) / degree;

        // r_mu := l
        // ---------
        rmu = l;

        // Form the KKT system
        // -------------------
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        FinishKKT
        ( m, n, w, 
          orders, firstInds, 
          origToSparseOrders, origToSparseFirstInds, 
          kSparse, JOrig, onlyLower );
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, 
          orders, firstInds, origToSparseFirstInds, kSparse, d );

        // Solve for the direction
        // -----------------------
        try 
        {
            J = JOrig;
            J.FreezeSparsity();
            UpdateDiagonal( J, Real(1), regTmp );
            if( wMaxNorm >= ctrl.ruizEquilTol )
                SymmetricRuizEquil( J, dInner, ctrl.ruizMaxIter, ctrl.print );
            else if( wMaxNorm >= ctrl.diagEquilTol )
                SymmetricDiagonalEquil( J, dInner, ctrl.print );
            else
                Ones( dInner, n+m+kSparse, 1 );

            JFront.Pull( J, map, info );
            LDL( info, JFront, LDL_2D );
            if( ctrl.resolveReg )
                reg_ldl::SolveAfter
                ( JOrig, regTmp, dInner, invMap, info, JFront, d,
                  ctrl.solveCtrl );
            else
                reg_ldl::RegularizedSolveAfter
                ( JOrig, regTmp, dInner, invMap, info, JFront, d,
                  ctrl.solveCtrl.relTol, ctrl.solveCtrl.maxRefineIts, 
                  ctrl.solveCtrl.progress );
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
        ( m, n, d, rmu, wRoot, 
          orders, firstInds, 
          sparseOrders, sparseFirstInds,
          sparseToOrigOrders, sparseToOrigFirstInds,
          dxAff, dyAff, dzAff, dsAff );
        soc::ApplyQuadratic( wRoot, dzAff, dzAffScaled, orders, firstInds );
        soc::ApplyQuadratic( wRootInv, dsAff, dsAffScaled, orders, firstInds );

        if( ctrl.checkResiduals && ctrl.print )
        {
            dxError = rb;
            Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
            const Real dxErrorNrm2 = Nrm2( dxError );

            dyError = rc;
            Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
            Multiply( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
            const Real dyErrorNrm2 = Nrm2( dyError );

            dzError = rh;
            Multiply( NORMAL, Real(1), G, dxAff, Real(1), dzError );
            dzError += dsAff;
            const Real dzErrorNrm2 = Nrm2( dzError );

            dmuError = dzAffScaled;
            dmuError += dsAffScaled;
            dmuError += l;
            const Real rmuNrm2 = Nrm2( rmu );
            const Real dmuErrorNrm2 = Nrm2( dmuError );
            // TODO: Also compute and print the residuals with regularization

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+rhNrm2),"\n",Indent(),
             "|| dmuError ||_2 / (1 + || r_mu ||_2) = ",
             dmuErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = 
          soc::MaxStep( s, dsAff, orders, firstInds, Real(1) );
        Real alphaAffDual = 
          soc::MaxStep( z, dzAff, orders, firstInds, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and ds are used as temporaries
        ds = s;
        dz = z;
        Axpy( alphaAffPri,  dsAff, ds );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(ds,dz) / degree;
        if( ctrl.print )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma = centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        rc *= 1-sigma;
        rb *= 1-sigma;
        rh *= 1-sigma;
        if( ctrl.mehrotra )
        {
            // r_mu := l + inv(l) o ((inv(W)^T dsAff) o (W dzAff) - sigma*mu)
            // --------------------------------------------------------------
            soc::Apply( dsAffScaled, dzAffScaled, rmu, orders, firstInds );
            soc::Shift( rmu, -sigma*mu, orders, firstInds );
            soc::Apply( lInv, rmu, orders, firstInds );
            rmu += l;
        }
        else
        {
            // r_mu -= sigma*mu*inv(l)
            // -----------------------
            Axpy( -sigma*mu, lInv, rmu );
        }

        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, 
          orders, firstInds, origToSparseFirstInds, kSparse, d );
        try 
        {
            if( ctrl.resolveReg )
                reg_ldl::SolveAfter
                ( JOrig, regTmp, dInner, invMap, info, JFront, d,
                  ctrl.solveCtrl );
            else
                reg_ldl::RegularizedSolveAfter
                ( JOrig, regTmp, dInner, invMap, info, JFront, d,
                  ctrl.solveCtrl.relTol, ctrl.solveCtrl.maxRefineIts, 
                  ctrl.solveCtrl.progress );
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
        ( m, n, d, rmu, wRoot, 
          orders, firstInds, 
          sparseOrders, sparseFirstInds,
          sparseToOrigOrders, sparseToOrigFirstInds,
          dx, dy, dz, ds );

        // Update the current estimates
        // ============================
        Real alphaPri = 
          soc::MaxStep( s, ds, orders, firstInds, 1/ctrl.maxStepRatio );
        Real alphaDual = 
          soc::MaxStep( z, dz, orders, firstInds, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
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
    SetIndent( indent );
    if( ctrl.outerEquil )
    {
        DiagonalSolve( LEFT, NORMAL, dCol,  x );
        DiagonalSolve( LEFT, NORMAL, dRowA, y );
        DiagonalSolve( LEFT, NORMAL, dRowG, z );
        DiagonalScale( LEFT, NORMAL, dRowG, s );
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
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y,
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_CSE
    const Real eps = limits::Epsilon<Real>();
    const bool onlyLower = false;

    // TODO: Move these into the control structur
    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const Int cutoffSparse = 64;
    const Int cutoffPar = 1000;
    const bool standardShift = false;
    const Real gamma = Pow(eps,Real(0.35));
    const Real delta = Pow(eps,Real(0.35));
    const Real beta =  Pow(eps,Real(0.35));
    const Real gammaTmp = Pow(eps,Real(0.25));
    const Real deltaTmp = Pow(eps,Real(0.25));
    const Real betaTmp  = Pow(eps,Real(0.25));

    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto h = hPre;
    auto c = cPre;

    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = soc::Degree( firstInds );
    mpi::Comm comm = APre.Comm();
    const int commRank = mpi::Rank(comm);
    Timer timer, iterTimer;

    DistMultiVec<Real> dRowA(comm), dRowG(comm), dCol(comm);
    if( ctrl.outerEquil )
    {
        if( commRank == 0 && ctrl.time )
            timer.Start();
        cone::RuizEquil
        ( A, G, dRowA, dRowG, dCol, orders, firstInds, cutoffPar, ctrl.print );
        if( commRank == 0 && ctrl.time )
            Output("cone::RuizEquil: ",timer.Stop()," secs");

        DiagonalSolve( LEFT, NORMAL, dRowA, b );
        DiagonalSolve( LEFT, NORMAL, dRowG, h );
        DiagonalSolve( LEFT, NORMAL, dCol,  c );
        if( ctrl.primalInit )
        {
            DiagonalScale( LEFT, NORMAL, dCol,  x );
            DiagonalSolve( LEFT, NORMAL, dRowG, s );
        }
        if( ctrl.dualInit )
        {
            DiagonalScale( LEFT, NORMAL, dRowA, y );
            DiagonalScale( LEFT, NORMAL, dRowG, z );
        }
    }
    else
    {
        Ones( dRowA, m, 1 );
        Ones( dRowG, k, 1 );
        Ones( dCol,  n, 1 );
    }

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );
    const Real twoNormEstA = TwoNormEstimate( A, ctrl.basisSize );
    const Real twoNormEstG = TwoNormEstimate( G, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstA + twoNormEstG + 1;
    if( ctrl.print )
    {
        const double imbalanceA = A.Imbalance();
        const double imbalanceG = G.Imbalance();
        if( commRank == 0 )
        {
            Output("|| A ||_2 estimate: ",twoNormEstA);
            Output("|| G ||_2 estimate: ",twoNormEstG);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
            Output("|| h ||_2 = ",hNrm2);
            Output("Imbalance factor of A: ",imbalanceA);
            Output("Imbalance factor of G: ",imbalanceG);
        }
    }

    if( commRank == 0 && ctrl.time )
        timer.Start();
    Initialize
    ( A, G, b, c, h, orders, firstInds, x, y, z, s,
      ctrl.primalInit, ctrl.dualInit, standardShift, cutoffPar, 
      ctrl.solveCtrl );
    if( commRank == 0 && ctrl.time )
        Output("Init: ",timer.Stop()," secs");

    // Form the offsets for the sparse embedding of the barrier's Hessian
    // ================================================================== 
    DistMultiVec<Int> 
      origToSparseOrders(comm),    sparseToOrigOrders(comm),
      origToSparseFirstInds(comm), sparseToOrigFirstInds(comm), 
      sparseOrders(comm), sparseFirstInds(comm);
    soc::EmbeddingMaps
    ( orders, firstInds,
      sparseOrders, sparseFirstInds,
      origToSparseOrders, origToSparseFirstInds,
      sparseToOrigOrders, sparseToOrigFirstInds,
      cutoffSparse );
    const Int kSparse = sparseFirstInds.Height();

    auto& sparseOrdersLoc = sparseOrders.LockedMatrix();
    auto& sparseFirstIndsLoc = sparseFirstInds.LockedMatrix();
    auto& sparseToOrigOrdersLoc = sparseToOrigOrders.LockedMatrix();

    DistSparseMultMeta metaOrig;
    DistSparseMatrix<Real> J(comm), JOrig(comm);
    ldl::DistFront<Real> JFront;
    DistMultiVec<Real> d(comm),
                       w(comm),     wRoot(comm), wRootInv(comm),
                       l(comm),     lInv(comm),
                       rc(comm),    rb(comm),    rh(comm),    rmu(comm),
                       dxAff(comm), dyAff(comm), dzAff(comm), dsAff(comm),
                       dx(comm),    dy(comm),    dz(comm),    ds(comm),
                       dzAffScaled(comm), dsAffScaled(comm);

    // Form the regularization vectors
    // ===============================
    DistMultiVec<Real> regTmp(comm);
    Zeros( regTmp, n+m+kSparse, 1 );
    // Set the analytical part
    // -----------------------
    for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
    {
        const Int i = regTmp.GlobalRow(iLoc);
        if( i < n )        regTmp.SetLocal( iLoc, 0,  gammaTmp*gammaTmp );
        else if( i < n+m ) regTmp.SetLocal( iLoc, 0, -deltaTmp*deltaTmp );
        else break;
    }
    // Perform the portion that requires remote updates
    // ------------------------------------------------
    {
        const Int sparseLocalHeight = sparseFirstInds.LocalHeight();
        regTmp.Reserve( sparseLocalHeight );
        for( Int iLoc=0; iLoc<sparseLocalHeight; ++iLoc )
        {
            const Int iCone = sparseFirstInds.GlobalRow(iLoc);
            const Int order = sparseToOrigOrdersLoc(iLoc);
            const Int sparseOrder = sparseOrdersLoc(iLoc);
            const bool embedded = ( order != sparseOrder ); 
            const Int firstInd = sparseFirstIndsLoc(iLoc);
            if( embedded && iCone == firstInd+sparseOrder-1 )
                regTmp.QueueUpdate( n+m+iCone, 0, betaTmp*betaTmp );
            else
                regTmp.QueueUpdate( n+m+iCone, 0, -betaTmp*betaTmp );
        }
        regTmp.ProcessQueues();
    }
    regTmp *= origTwoNormEst;

    // Form the static portion of the KKT system
    // =========================================
    DistSparseMatrix<Real> JStatic(comm);
    StaticKKT
    ( A, G, gamma, delta, beta, 
      orders, firstInds, origToSparseOrders, origToSparseFirstInds, 
      kSparse, JStatic, onlyLower );
    if( ctrl.print )
    {
        const Real imbalanceJ = JStatic.Imbalance();
        if( commRank == 0 )
            Output("Imbalance factor of J: ",imbalanceJ);
    }

    auto meta = JStatic.InitializeMultMeta();
    if( commRank == 0 && ctrl.time )
        timer.Start();
    DistMap map, invMap;
    ldl::DistNodeInfo info;
    ldl::DistSeparator rootSep;
    NestedDissection( JStatic.LockedDistGraph(), map, rootSep, info );
    if( commRank == 0 && ctrl.time )
        Output("ND: ",timer.Stop()," secs");
    InvertMap( map, invMap );

    vector<Int> mappedSources, mappedTargets, colOffs;
    JStatic.MappedSources( map, mappedSources );
    JStatic.MappedTargets( map, mappedTargets, colOffs );

    Real relError = 1;
    DistMultiVec<Real> dInner(comm);
    DistMultiVec<Real> dxError(comm), dyError(comm), 
                       dzError(comm), dmuError(comm);
    ldl::DistMultiVecNodeMeta dmvMeta;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        if( ctrl.time && commRank == 0 )
            iterTimer.Start();
        // Ensure that s and z are in the cone
        // ===================================
        // TODO: Let this be a function of the relative error, etc.
        const Real minDist = eps;
        soc::PushInto( s, orders, firstInds, minDist, cutoffPar );
        soc::PushInto( z, orders, firstInds, minDist, cutoffPar );
        soc::NesterovTodd( s, z, w, orders, firstInds, cutoffPar );

        // Check for convergence
        // =====================
        // |c^T x - (-b^T y - h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        rb *= -1;
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h;
        rh *= -1;
        Multiply( NORMAL, Real(1), G, x, Real(1), rh );
        rh += s;
        const Real rhNrm2 = Nrm2( rh );
        const Real rhConv = rhNrm2 / (1+hNrm2);

        // Now check the pieces
        // --------------------
        relError = Max(Max(Max(objConv,rbConv),rcConv),rhConv);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( x );
            const Real yNrm2 = Nrm2( y );
            const Real zNrm2 = Nrm2( z );
            const Real sNrm2 = Nrm2( s );
            if( commRank == 0 )
                Output
                ("iter ",numIts,":\n",Indent(),
                 "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
                 "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
                 "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
                 "  ||  s  ||_2 = ",sNrm2,"\n",Indent(),
                 "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
                 "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
                 "  || r_h ||_2 = ",rhNrm2,"\n",Indent(),
                 "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
                 "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
                 "  || r_h ||_2 / (1 + || h ||_2) = ",rhConv,"\n",Indent(),
                 "  primal = ",primObj,"\n",Indent(),
                 "  dual   = ",dualObj,"\n",Indent(),
                 "  |primal - dual| / (1 + |primal|) = ",objConv);
        }
        if( relError <= ctrl.targetTol )
            break;
        if( numIts == ctrl.maxIts && relError > ctrl.minTol )
            RuntimeError
            ("Reached maximum number of iterations, ",ctrl.maxIts,
             ", with rel. error ",relError," which does not meet the minimum ",
             "tolerance of ",ctrl.minTol);
        Real wMaxNorm = MaxNorm(w);
        if( wMaxNorm >= ctrl.wMaxLimit && relError <= ctrl.minTol )
            break;

        // Compute the affine search direction
        // ===================================
        const Real wMaxNormLimit =
          Max(ctrl.wSafeMaxNorm,10/Min(Real(1),relError));
        if( ctrl.print && commRank == 0 )
            Output
            ("|| w ||_max = ",wMaxNorm,", limit is ",wMaxNormLimit);
        if( wMaxNorm > wMaxNormLimit )
        {
            if( ctrl.print && commRank == 0 )
                Output
                ("|| w ||_max = ",wMaxNorm," was larger than ",wMaxNormLimit);
            soc::PushPairInto
            ( s, z, w, orders, firstInds, wMaxNormLimit, cutoffPar );
            soc::NesterovTodd( s, z, w, orders, firstInds, cutoffPar );
            wMaxNorm = MaxNorm(w);
            if( ctrl.print && commRank == 0 )
                Output("New || w ||_max = ",wMaxNorm);
        }
        soc::SquareRoot( w, wRoot, orders, firstInds, cutoffPar );
        soc::Inverse( wRoot, wRootInv, orders, firstInds, cutoffPar );
        soc::ApplyQuadratic( wRoot, z, l, orders, firstInds, cutoffPar );
        soc::Inverse( l, lInv, orders, firstInds, cutoffPar );
        const Real mu = Dot(s,z) / degree;

        // r_mu := l
        // ---------
        rmu = l;

        // Construct the KKT system
        // ------------------------
        if( ctrl.time && commRank == 0 )
            timer.Start();
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        FinishKKT
        ( m, n, w, 
          orders, firstInds, 
          origToSparseOrders, origToSparseFirstInds,
          kSparse, JOrig, onlyLower, cutoffPar );
        if( ctrl.time && commRank == 0 )
            Output("KKT construction: ",timer.Stop()," secs");
        if( ctrl.time && commRank == 0 )
            timer.Start();
        KKTRHS
        ( rc, rb, rh, rmu, wRoot,
          orders, firstInds, origToSparseFirstInds, kSparse,
          d, cutoffPar );
        if( ctrl.time && commRank == 0 )
            Output("KKTRHS construction: ",timer.Stop()," secs");

        // Solve for the direction
        // -----------------------
        try
        {
            // Cache the metadata for the finalized JOrig
            JOrig.multMeta = meta;
            J = JOrig;
            J.FreezeSparsity();
            UpdateDiagonal( J, Real(1), regTmp );

            if( commRank == 0 && ctrl.time )
                timer.Start();
            if( wMaxNorm >= ctrl.ruizEquilTol )
                SymmetricRuizEquil( J, dInner, ctrl.ruizMaxIter, ctrl.print );
            else if( wMaxNorm >= ctrl.diagEquilTol )
                SymmetricDiagonalEquil( J, dInner, ctrl.print );
            else
                Ones( dInner, n+m+kSparse, 1 );
            if( commRank == 0 && ctrl.time )
                Output("Equilibration: ",timer.Stop()," secs");

            // Cache the metadata for the finalized J
            J.multMeta = meta;
            if( ctrl.time && commRank == 0 )
                timer.Start();
            JFront.Pull
            ( J, map, rootSep, info, mappedSources, mappedTargets, colOffs );
            if( ctrl.time && commRank == 0 )
                Output("Front pull: ",timer.Stop()," secs");

            if( commRank == 0 && ctrl.time )
                timer.Start();
            LDL( info, JFront, LDL_2D );
            if( commRank == 0 && ctrl.time )
                Output("LDL: ",timer.Stop()," secs");

            if( commRank == 0 && ctrl.time )
                timer.Start();
            if( ctrl.resolveReg )
                reg_ldl::SolveAfter
                ( JOrig, regTmp, dInner, invMap, info, JFront, d, dmvMeta, 
                  ctrl.solveCtrl );
            else
                reg_ldl::RegularizedSolveAfter
                ( JOrig, regTmp, dInner, invMap, info, JFront, d, dmvMeta,
                  ctrl.solveCtrl.relTol, ctrl.solveCtrl.maxRefineIts, 
                  ctrl.solveCtrl.progress );
            if( commRank == 0 && ctrl.time )
                Output("Affine: ",timer.Stop()," secs");
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
        ( m, n, d, rmu, wRoot, 
          orders, firstInds, 
          sparseOrders, sparseFirstInds,
          sparseToOrigOrders, sparseToOrigFirstInds,
          dxAff, dyAff, dzAff, dsAff, cutoffPar );
        soc::ApplyQuadratic
        ( wRoot, dzAff, dzAffScaled, orders, firstInds, cutoffPar );
        soc::ApplyQuadratic
        ( wRootInv, dsAff, dsAffScaled, orders, firstInds, cutoffPar );

        if( ctrl.checkResiduals && ctrl.print )
        {
            if( ctrl.time && commRank == 0 )
                timer.Start();
            dxError = rb;
            Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
            const Real dxErrorNrm2 = Nrm2( dxError );

            dyError = rc;
            Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
            Multiply( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
            const Real dyErrorNrm2 = Nrm2( dyError );

            dzError = rh;
            Multiply( NORMAL, Real(1), G, dxAff, Real(1), dzError );
            dzError += dsAff;
            const Real dzErrorNrm2 = Nrm2( dzError );

            dmuError = dzAffScaled;
            dmuError += dsAffScaled;
            dmuError += l;
            const Real rmuNrm2 = Nrm2( rmu );
            const Real dmuErrorNrm2 = Nrm2( dmuError );
            // TODO: Also compute and print the residuals with regularization
        
            if( commRank == 0 )
                Output
                ("||  r_mu    ||_2 = ",rmuNrm2,"\n",Indent(),
                 "|| dxError  ||_2 = ",dxErrorNrm2,"\n",Indent(),
                 "|| dyError  ||_2 = ",dyErrorNrm2,"\n",Indent(),
                 "|| dzError  ||_2 = ",dzErrorNrm2,"\n",Indent(),
                 "|| dmuError ||_2 = ",dmuErrorNrm2,"\n",Indent(),
                 "|| dxError ||_2 / (1 + || r_b ||_2) = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rhNrm2),"\n",Indent(),
                 "|| dmuError ||_2 / (1 + || r_mu ||_2) = ",
                 dmuErrorNrm2/(1+rmuNrm2));
            if( ctrl.time && commRank == 0 )
                Output("residual check: ",timer.Stop()," secs");
        }

        // Compute a centrality parameter
        // ==============================
        if( ctrl.time && commRank == 0 )
            timer.Start();
        Real alphaAffPri = 
          soc::MaxStep( s, dsAff, orders, firstInds, Real(1), cutoffPar );
        Real alphaAffDual = 
          soc::MaxStep( z, dzAff, orders, firstInds, Real(1), cutoffPar );
        if( ctrl.time && commRank == 0 )
            Output("Affine line search: ",timer.Stop()," secs");
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and ds are used as temporaries
        ds = s;
        dz = z;
        Axpy( alphaAffPri,  dsAff, ds );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(ds,dz) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma = centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        rc *= 1-sigma;
        rb *= 1-sigma;
        rh *= 1-sigma;
        if( ctrl.time && commRank == 0 )
            timer.Start();
        if( ctrl.mehrotra )
        {
            // r_mu := l + inv(l) o ((inv(W)^T dsAff) o (W dzAff) - sigma*mu)
            // --------------------------------------------------------------
            soc::Apply( dsAffScaled, dzAffScaled, rmu, orders, firstInds );
            soc::Shift( rmu, -sigma*mu, orders, firstInds );
            soc::Apply( lInv, rmu, orders, firstInds );
            rmu += l;
        }
        else
        {
            // r_mu -= sigma*mu*inv(l)
            // -----------------------
            Axpy( -sigma*mu, lInv, rmu );
        }
        if( ctrl.time && commRank == 0 )
            Output("r_mu formation: ",timer.Stop()," secs");

        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS
        ( rc, rb, rh, rmu, wRoot, 
          orders, firstInds, origToSparseFirstInds, kSparse,
          d, cutoffPar );
        try
        {
            if( commRank == 0 && ctrl.time )
                timer.Start();
            if( ctrl.resolveReg )
                reg_ldl::SolveAfter
                ( JOrig, regTmp, dInner, invMap, info, JFront, d, dmvMeta,
                  ctrl.solveCtrl );
            else
                reg_ldl::RegularizedSolveAfter
                ( JOrig, regTmp, dInner, invMap, info, JFront, d, dmvMeta,
                  ctrl.solveCtrl.relTol, ctrl.solveCtrl.maxRefineIts, 
                  ctrl.solveCtrl.progress );
            if( commRank == 0 && ctrl.time )
                Output("Corrector solver: ",timer.Stop()," secs");
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
        if( ctrl.time && commRank == 0 )
            timer.Start();
        ExpandSolution
        ( m, n, d, rmu, wRoot, 
          orders, firstInds, 
          sparseOrders, sparseFirstInds,
          sparseToOrigOrders, sparseToOrigFirstInds,
          dx, dy, dz, ds, cutoffPar );
        if( ctrl.time && commRank == 0 )
            Output("ExpandSolution: ",timer.Stop()," secs");

        // Update the current estimates
        // ============================
        if( ctrl.time && commRank == 0 )
            timer.Start();
        Real alphaPri = 
          soc::MaxStep
          ( s, ds, orders, firstInds, 1/ctrl.maxStepRatio, cutoffPar );
        Real alphaDual = 
          soc::MaxStep
          ( z, dz, orders, firstInds, 1/ctrl.maxStepRatio, cutoffPar );
        if( ctrl.time && commRank == 0 )
            Output("Combined line search: ",timer.Stop()," secs");
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
        if( ctrl.time && commRank == 0 )
            Output("Iteration: ",iterTimer.Stop()," secs");
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
    SetIndent( indent );
    if( ctrl.outerEquil )
    {
        DiagonalSolve( LEFT, NORMAL, dCol,  x );
        DiagonalSolve( LEFT, NORMAL, dRowA, y );
        DiagonalSolve( LEFT, NORMAL, dRowG, z );
        DiagonalScale( LEFT, NORMAL, dRowG, s );
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
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& G, \
    const ElementalMatrix<Real>& b, \
    const ElementalMatrix<Real>& c, \
    const ElementalMatrix<Real>& h, \
    const ElementalMatrix<Int>& orders, \
    const ElementalMatrix<Int>& firstInds, \
          ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
    const Matrix<Int>& orders, \
    const Matrix<Int>& firstInds, \
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
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
          DistMultiVec<Real>& s, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace affine
} // namespace socp
} // namespace El
