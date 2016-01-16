/*
   Copyright (c) 2009-2016, Jack Poulson
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
( const Matrix<Real>& QPre,
  const Matrix<Real>& APre,
  const Matrix<Real>& GPre,
  const Matrix<Real>& bPre,
  const Matrix<Real>& cPre,
  const Matrix<Real>& hPre,
        Matrix<Real>& x,
        Matrix<Real>& y, 
        Matrix<Real>& z,
        Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("qp::affine::Mehrotra"))    

    // TODO: Move these into the control structure
    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const bool standardShift = true;

    // Equilibrate the QP by diagonally scaling [A;G]
    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto c = cPre;
    auto h = hPre;
    auto Q = QPre;
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = k;
    Matrix<Real> dRowA, dRowG, dCol;
    if( ctrl.outerEquil )
    {
        StackedRuizEquil( A, G, dRowA, dRowG, dCol, ctrl.print );
        DiagonalSolve( LEFT, NORMAL, dRowA, b );
        DiagonalSolve( LEFT, NORMAL, dRowG, h );
        DiagonalSolve( LEFT, NORMAL, dCol,  c );
        // TODO: Replace with SymmetricDiagonalSolve
        {
            DiagonalSolve( LEFT, NORMAL, dCol,  Q );
            DiagonalSolve( RIGHT, NORMAL, dCol, Q );
        }
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
        const Real QNrm1 = HermitianOneNorm( LOWER, Q );
        const Real ANrm1 = OneNorm( A );
        const Real GNrm1 = OneNorm( G );
        Output("|| Q ||_1 = ",QNrm1); 
        Output("|| A ||_1 = ",ANrm1);
        Output("|| G ||_1 = ",GNrm1);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
        Output("|| h ||_2 = ",hNrm2);
    }

    Initialize
    ( Q, A, G, b, c, h, x, y, z, s, 
      ctrl.primalInit, ctrl.dualInit, standardShift );

    Real relError = 1;
    Matrix<Real> J, d,
                 rmu,   rc,    rb,    rh,
                 dxAff, dyAff, dzAff, dsAff,
                 dx,    dy,    dz,    ds;
    Matrix<Real> dSub;
    Permutation p;
    Matrix<Real> dxError, dyError, dzError;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = pos_orth::NumOutside( s );
        const Int zNumNonPos = pos_orth::NumOutside( z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure
        // ===========================
        const Real mu = Dot(s,z) / k;

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        Zeros( d, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj =  xTQx/2 + Dot(c,x);
        const Real dualObj = -xTQx/2 - Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b; rb *= -1;
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Hemv( LOWER,     Real(1), Q, x, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h; rh *= -1;
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
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := s o z
        // -------------
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );

        // Construct the full KKT system
        // -----------------------------
        KKT( Q, A, G, s, z, J );
        KKTRHS( rc, rb, rh, rmu, z, d );

        // Compute the proposed step from the KKT system
        // ---------------------------------------------
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
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );

        if( ctrl.checkResiduals && ctrl.print )
        {
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
            dzError += dsAff;
            const Real dzErrorNrm2 = Nrm2( dzError );

            // TODO: dmuError

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+rhNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = pos_orth::MaxStep( s, dsAff, Real(1) );
        Real alphaAffDual = pos_orth::MaxStep( z, dzAff, Real(1) );
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
        Shift( rmu, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dsAff o dzAff
            // ---------------------
            // NOTE: dz is used as a temporary
            dz = dzAff;
            DiagonalScale( LEFT, NORMAL, dsAff, dz );
            rmu += dz;
        }

        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        KKTRHS( rc, rb, rh, rmu, z, d );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = pos_orth::MaxStep( s, ds, 1/ctrl.maxStepRatio );
        Real alphaDual = pos_orth::MaxStep( z, dz, 1/ctrl.maxStepRatio );
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
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
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
( const ElementalMatrix<Real>& QPre, 
  const ElementalMatrix<Real>& APre,
  const ElementalMatrix<Real>& GPre,
  const ElementalMatrix<Real>& bPre,
  const ElementalMatrix<Real>& cPre,
  const ElementalMatrix<Real>& hPre,
        ElementalMatrix<Real>& xPre,
        ElementalMatrix<Real>& yPre, 
        ElementalMatrix<Real>& zPre,
        ElementalMatrix<Real>& sPre,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("qp::affine::Mehrotra"))    

    // TODO: Move these into the control structure
    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const bool standardShift = true;

    const Grid& grid = APre.Grid();
    const int commRank = grid.Rank();
    Timer timer;

    // Ensure that the inputs have the appropriate read/write properties
    DistMatrix<Real> Q(grid), A(grid), G(grid), b(grid), c(grid), h(grid);
    Q.Align(0,0);
    A.Align(0,0);
    G.Align(0,0);
    b.Align(0,0);
    c.Align(0,0);
    Q = QPre;
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

    // NOTE: {x,s} do not need to be a read proxy when !ctrl.primalInit
    DistMatrixReadWriteProxy<Real,Real,MC,MR>
      xProx( xPre, control ),
      sProx( sPre, control ),
    // NOTE: {y,z} do not need to be read proxies when !ctrl.dualInit
      yProx( yPre, control ),
      zProx( zPre, control );
    auto& x = xProx.Get();
    auto& s = sProx.Get();
    auto& y = yProx.Get();
    auto& z = zProx.Get();

    // Equilibrate the QP by diagonally scaling [A;G]
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = k;
    DistMatrix<Real,MC,STAR> dRowA(grid),
                             dRowG(grid);
    DistMatrix<Real,MR,STAR> dCol(grid);
    if( ctrl.outerEquil )
    {
        if( ctrl.time && commRank == 0 )
            timer.Start();
        StackedRuizEquil( A, G, dRowA, dRowG, dCol, ctrl.print );
        if( ctrl.time && commRank == 0 )
            Output("RuizEquil: ",timer.Stop()," secs");
        DiagonalSolve( LEFT, NORMAL, dRowA, b );
        DiagonalSolve( LEFT, NORMAL, dRowG, h );
        DiagonalSolve( LEFT, NORMAL, dCol,  c );
        // TODO: Replace with SymmetricDiagonalSolve
        {
            DiagonalSolve( LEFT, NORMAL, dCol,  Q );
            DiagonalSolve( RIGHT, NORMAL, dCol, Q );
        }
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
        const Real QNrm1 = HermitianOneNorm( LOWER, Q );
        const Real ANrm1 = OneNorm( A );
        const Real GNrm1 = OneNorm( G );
        if( commRank == 0 )
        {
            Output("|| Q ||_1 = ",QNrm1);
            Output("|| A ||_1 = ",ANrm1);
            Output("|| G ||_1 = ",GNrm1);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
            Output("|| h ||_2 = ",hNrm2);
        }
    }

    if( ctrl.time && commRank == 0 )
        timer.Start();
    Initialize
    ( Q, A, G, b, c, h, x, y, z, s, 
      ctrl.primalInit, ctrl.dualInit, standardShift );
    if( ctrl.time && commRank == 0 )
        Output("Init time: ",timer.Stop()," secs");

    Real relError = 1;
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
    DistPermutation p(grid);
    DistMatrix<Real> dxError(grid), dyError(grid), dzError(grid);
    dzError.AlignWith( s );
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = pos_orth::NumOutside( s );
        const Int zNumNonPos = pos_orth::NumOutside( z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure
        // ===========================
        const Real mu = Dot(s,z) / k;

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        Zeros( d, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj =  xTQx/2 + Dot(c,x);
        const Real dualObj = -xTQx/2 - Dot(b,y) - Dot(h,z);
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b; rb *= -1;
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Hemv( LOWER,     Real(1), Q, x, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h; rh *= -1;
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
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := s o z
        // -------------
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );

        // Construct the KKT system
        // ------------------------
        KKT( Q, A, G, s, z, J );
        KKTRHS( rc, rb, rh, rmu, z, d );

        // Solve for the direction
        // -----------------------
        try
        {
            if( ctrl.time && commRank == 0 )
                timer.Start();
            LDL( J, dSub, p, false );
            if( ctrl.time && commRank == 0 )
            {
                Output("LDL: ",timer.Stop()," secs");
                timer.Start();
            }
            ldl::SolveAfter( J, dSub, p, d, false );
            if( ctrl.time && commRank == 0 )
                Output("Affine solve: ",timer.Stop()," secs"); 
        }
        catch(...)
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );

        if( ctrl.checkResiduals && ctrl.print )
        {
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
            dzError += dsAff;
            const Real dzErrorNrm2 = Nrm2( dzError );

            // TODO: dmuError

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rhNrm2));
        }
 
        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = pos_orth::MaxStep( s, dsAff, Real(1) );
        Real alphaAffDual = pos_orth::MaxStep( z, dzAff, Real(1) );
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
        Shift( rmu, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dsAff o dzAff
            // ---------------------
            // NOTE: dz is used as a temporary
            dz = dzAff;
            DiagonalScale( LEFT, NORMAL, dsAff, dz );
            rmu += dz;
        }

        // Form the new KKT RHS
        // --------------------
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Solve for the new direction
        // ---------------------------
        try 
        { 
            if( ctrl.time && commRank == 0 )
                timer.Start();
            ldl::SolveAfter( J, dSub, p, d, false ); 
            if( ctrl.time && commRank == 0 )
                Output("Combined solve: ",timer.Stop()," secs");
        }
        catch(...)
        {
            if( relError <= ctrl.minTol )    
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = pos_orth::MaxStep( s, ds, 1/ctrl.maxStepRatio );
        Real alphaDual = pos_orth::MaxStep( z, dz, 1/ctrl.maxStepRatio );
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
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
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
( const SparseMatrix<Real>& QPre,
  const SparseMatrix<Real>& APre,
  const SparseMatrix<Real>& GPre,
  const Matrix<Real>& bPre,
  const Matrix<Real>& cPre,
  const Matrix<Real>& hPre,
        Matrix<Real>& x,
        Matrix<Real>& y, 
        Matrix<Real>& z,
        Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("qp::affine::Mehrotra"))    
    const Real eps = limits::Epsilon<Real>();

    // TODO: Move these into the control structure
    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const bool standardShift = true;
    const Real gamma = Pow(eps,Real(0.35));
    const Real delta = Pow(eps,Real(0.35));
    const Real beta  = Pow(eps,Real(0.35));
    const Real gammaTmp = Pow(eps,Real(0.25));
    const Real deltaTmp = Pow(eps,Real(0.25));
    const Real betaTmp  = Pow(eps,Real(0.25));

    // Equilibrate the QP by diagonally scaling [A;G]
    auto Q = QPre;
    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto c = cPre;
    auto h = hPre;
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = k;
    Matrix<Real> dRowA, dRowG, dCol;
    if( ctrl.outerEquil )
    {
        StackedRuizEquil( A, G, dRowA, dRowG, dCol, ctrl.print );
        DiagonalSolve( LEFT, NORMAL, dRowA, b );
        DiagonalSolve( LEFT, NORMAL, dRowG, h );
        DiagonalSolve( LEFT, NORMAL, dCol,  c );
        // TODO: Replace with SymmetricDiagonalSolve
        {
            DiagonalSolve( LEFT, NORMAL, dCol, Q );
            DiagonalSolve( RIGHT, NORMAL, dCol, Q );
        }
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
    const Real twoNormEstQ = HermitianTwoNormEstimate( Q, ctrl.basisSize );
    const Real twoNormEstA = TwoNormEstimate( A, ctrl.basisSize );
    const Real twoNormEstG = TwoNormEstimate( G, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstA + twoNormEstG + twoNormEstQ + 1;
    if( ctrl.print )
    {
        Output("|| Q ||_2 estimate: ",twoNormEstQ);
        Output("|| A ||_2 estimate: ",twoNormEstA);
        Output("|| G ||_2 estimate: ",twoNormEstG);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
        Output("|| h ||_2 = ",hNrm2);
    }

    // TODO: Expose regularization rules to user
    Matrix<Real> regTmp;
    regTmp.Resize( n+m+k, 1 );
    for( Int i=0; i<n+m+k; ++i )
    {
        if( i < n )        regTmp.Set( i, 0,  gammaTmp*gammaTmp );
        else if( i < n+m ) regTmp.Set( i, 0, -deltaTmp*deltaTmp );
        else               regTmp.Set( i, 0, -betaTmp*betaTmp );
    }
    regTmp *= origTwoNormEst;

    // Initialize the static portion of the KKT system
    // ===============================================
    SparseMatrix<Real> JStatic;
    StaticKKT( Q, A, G, gamma, delta, beta, JStatic, false );
    vector<Int> map, invMap;
    ldl::NodeInfo info;
    ldl::Separator rootSep;
    NestedDissection( JStatic.LockedGraph(), map, rootSep, info );
    InvertMap( map, invMap );

    Initialize
    ( JStatic, regTmp, b, c, h, x, y, z, s, map, invMap, rootSep, info, 
      ctrl.primalInit, ctrl.dualInit, standardShift, ctrl.solveCtrl );

    SparseMatrix<Real> J, JOrig;
    ldl::Front<Real> JFront;
    Matrix<Real> d,
                 w,
                 rc,    rb,    rh,    rmu,
                 dxAff, dyAff, dzAff, dsAff,
                 dx,    dy,    dz,    ds;

    Real relError = 1;
    Matrix<Real> dInner;
    Matrix<Real> dxError, dyError, dzError;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = pos_orth::NumOutside( s );
        const Int zNumNonPos = pos_orth::NumOutside( z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure and scaling point
        // =============================================
        const Real mu = Dot(s,z) / k;
        pos_orth::NesterovTodd( s, z, w );
        const Real wMaxNorm = MaxNorm( w );

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
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b; rb *= -1;
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        Axpy( -delta*delta, y, rb );
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( NORMAL,    Real(1), Q, x, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        Axpy( gamma*gamma, x, rc );
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h; rh *= -1;
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
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without ",
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := s o z
        // -------------
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );

        // Construct the KKT system
        // ------------------------
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        FinishKKT( m, n, s, z, JOrig );
        KKTRHS( rc, rb, rh, rmu, z, d );

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
                Ones( dInner, n+m+k, 1 );

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
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );

        if( ctrl.checkResiduals && ctrl.print )
        {
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
            dzError += dsAff;
            const Real dzErrorNrm2 = Nrm2( dzError );

            // TODO: dmuError
            // TODO: Also compute and print the residuals with regularization

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2) = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_h ||_2) = ",
             dzErrorNrm2/(1+rhNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = pos_orth::MaxStep( s, dsAff, Real(1) );
        Real alphaAffDual = pos_orth::MaxStep( z, dzAff, Real(1) );
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
        Shift( rmu, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dsAff o dzAff
            // ---------------------
            // NOTE: dz is used as a temporary
            dz = dzAff;
            DiagonalScale( LEFT, NORMAL, dsAff, dz );
            rmu += dz;
        }

        // Set up the new KKT RHS
        // ----------------------
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Solve for the new direction
        // ---------------------------
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
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );

        // Update the current estimates
        // ============================
        Real alphaPri = pos_orth::MaxStep( s, ds, 1/ctrl.maxStepRatio );
        Real alphaDual = pos_orth::MaxStep( z, dz, 1/ctrl.maxStepRatio );
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
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
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
( const DistSparseMatrix<Real>& QPre,
  const DistSparseMatrix<Real>& APre,
  const DistSparseMatrix<Real>& GPre,
  const DistMultiVec<Real>& bPre,
  const DistMultiVec<Real>& cPre,
  const DistMultiVec<Real>& hPre,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
        DistMultiVec<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CSE cse("qp::affine::Mehrotra"))    
    const Real eps = limits::Epsilon<Real>();

    // TODO: Move these into the control structure
    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const bool standardShift = true;
    const Real gamma = Pow(eps,Real(0.35));
    const Real delta = Pow(eps,Real(0.35));
    const Real beta  = Pow(eps,Real(0.35));
    const Real gammaTmp = Pow(eps,Real(0.25));
    const Real deltaTmp = Pow(eps,Real(0.25));
    const Real betaTmp  = Pow(eps,Real(0.25));
    //const Real selInvTol = Pow(eps,Real(-0.25));
    const Real selInvTol = 0;

    mpi::Comm comm = APre.Comm();
    const int commRank = mpi::Rank(comm);
    Timer timer, iterTimer;

    // Equilibrate the QP by diagonally scaling [A;G]
    auto Q = QPre;
    auto A = APre;
    auto G = GPre;
    auto b = bPre;
    auto h = hPre;
    auto c = cPre;
    const Int m = A.Height();
    const Int k = G.Height();
    const Int n = A.Width();
    const Int degree = k;
    DistMultiVec<Real> dRowA(comm), dRowG(comm), dCol(comm);
    if( ctrl.outerEquil )
    {
        if( commRank == 0 && ctrl.time )
            timer.Start();
        StackedRuizEquil( A, G, dRowA, dRowG, dCol, ctrl.print );
        if( commRank == 0 && ctrl.time )
            Output("RuizEquil: ",timer.Stop()," secs");

        DiagonalSolve( LEFT, NORMAL, dRowA, b );
        DiagonalSolve( LEFT, NORMAL, dRowG, h );
        DiagonalSolve( LEFT, NORMAL, dCol,  c );
        // TODO: Replace with SymmetricDiagonalSolve
        {
            DiagonalSolve( LEFT, NORMAL, dCol, Q );
            DiagonalSolve( RIGHT, NORMAL, dCol, Q );
        }
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
    const Real twoNormEstQ = HermitianTwoNormEstimate( Q, ctrl.basisSize );
    const Real twoNormEstA = TwoNormEstimate( A, ctrl.basisSize );
    const Real twoNormEstG = TwoNormEstimate( G, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstA + twoNormEstG + twoNormEstQ + 1;
    if( ctrl.print )
    {
        const double imbalanceQ = Q.Imbalance();
        const double imbalanceA = A.Imbalance();
        const double imbalanceG = G.Imbalance();
        if( commRank == 0 )
        {
            Output("|| Q ||_2 estimate: ",twoNormEstQ);
            Output("|| A ||_2 estimate: ",twoNormEstA);
            Output("|| G ||_2 estimate: ",twoNormEstG);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
            Output("|| h ||_2 = ",hNrm2);
            Output("Imbalance factor of Q: ",imbalanceQ);
            Output("Imbalance factor of A: ",imbalanceA);
            Output("Imbalance factor of G: ",imbalanceG);
        }
    }

    DistMultiVec<Real> regTmp(comm);
    regTmp.Resize( n+m+k, 1 );
    for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
    {
        const Int i = regTmp.GlobalRow(iLoc);
        if( i < n )        regTmp.SetLocal( iLoc, 0,  gammaTmp*gammaTmp );
        else if( i < n+m ) regTmp.SetLocal( iLoc, 0, -deltaTmp*deltaTmp );
        else               regTmp.SetLocal( iLoc, 0, -betaTmp*betaTmp );
    }
    regTmp *= origTwoNormEst;

    // Compute the static portion of the KKT system
    // ============================================
    DistSparseMatrix<Real> JStatic(comm);
    StaticKKT( Q, A, G, gamma, delta, beta, JStatic, false );
    JStatic.InitializeMultMeta();
    if( ctrl.print )
    {
        const double imbalanceJ = JStatic.Imbalance();
        if( commRank == 0 )
            Output("Imbalance factor of J: ",imbalanceJ);
    }

    DistMap map, invMap;
    ldl::DistNodeInfo info;
    ldl::DistSeparator rootSep;
    if( commRank == 0 && ctrl.time )
        timer.Start();
    NestedDissection( JStatic.LockedDistGraph(), map, rootSep, info );
    if( commRank == 0 && ctrl.time )
        Output("ND: ",timer.Stop()," secs");
    InvertMap( map, invMap );

    vector<Int> mappedSources, mappedTargets, colOffs;
    JStatic.MappedSources( map, mappedSources );
    JStatic.MappedTargets( map, mappedTargets, colOffs );

    if( commRank == 0 && ctrl.time )
        timer.Start();
    Initialize
    ( JStatic, regTmp, b, c, h, x, y, z, s, 
      map, invMap, rootSep, info, mappedSources, mappedTargets, colOffs, 
      ctrl.primalInit, ctrl.dualInit, standardShift, ctrl.solveCtrl );
    if( commRank == 0 && ctrl.time )
        Output("Init: ",timer.Stop()," secs");

    DistSparseMatrix<Real> J(comm), JOrig(comm);
    ldl::DistFront<Real> JFront;
    DistMultiVec<Real> d(comm),
                       w(comm),
                       rc(comm),    rb(comm),    rh(comm),    rmu(comm),
                       dxAff(comm), dyAff(comm), dzAff(comm), dsAff(comm),
                       dx(comm),    dy(comm),    dz(comm),    ds(comm);

    Real relError = 1;
    DistMultiVec<Real> dInner(comm);
    DistMultiVec<Real> dxError(comm), dyError(comm), dzError(comm);
    ldl::DistMultiVecNodeMeta dmvMeta;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        if( ctrl.time && commRank == 0 )
            iterTimer.Start();

        // Ensure that s and z are in the cone
        // ===================================
        const Int sNumNonPos = pos_orth::NumOutside( s );
        const Int zNumNonPos = pos_orth::NumOutside( z );
        if( sNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (sNumNonPos," entries of s were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure and scaling point
        // =============================================
        const Real mu = Dot(s,z) / k;
        pos_orth::NesterovTodd( s, z, w );
        const Real wMaxNorm = MaxNorm( w );

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
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b; rb *= -1;
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        Axpy( -delta*delta, y, rb );
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( NORMAL,    Real(1), Q, x, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), G, z, Real(1), rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        Axpy( gamma*gamma, x, rc );
        // || r_h ||_2 / (1 + || h ||_2) <= tol
        // ------------------------------------
        rh = h; rh *= -1;
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
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded without "
             "achieving minTol=",ctrl.minTol);

        // Compute the affine search direction
        // ===================================

        // r_mu := s o z
        // -------------
        rmu = z;
        DiagonalScale( LEFT, NORMAL, s, rmu );

        // Construct the KKT system
        // ------------------------
        JOrig = JStatic;
        JOrig.FreezeSparsity();
        JOrig.multMeta = JStatic.multMeta;
        FinishKKT( m, n, s, z, JOrig );
        KKTRHS( rc, rb, rh, rmu, z, d );

        // Solve for the direction
        // -----------------------
        try
        {
            J = JOrig;
            J.FreezeSparsity();
            J.multMeta = JStatic.multMeta;
            UpdateDiagonal( J, Real(1), regTmp );

            if( commRank == 0 && ctrl.time )
                timer.Start();
            if( wMaxNorm >= ctrl.ruizEquilTol )
                SymmetricRuizEquil( J, dInner, ctrl.ruizMaxIter, ctrl.print );
            else if( wMaxNorm >= ctrl.diagEquilTol )
                SymmetricDiagonalEquil( J, dInner, ctrl.print );
            else 
                Ones( dInner, n+m+k, 1 );
            if( commRank == 0 && ctrl.time )
                Output("Equilibration: ",timer.Stop()," secs");

            JFront.Pull
            ( J, map, rootSep, info, mappedSources, mappedTargets, colOffs );

            if( commRank == 0 && ctrl.time )
                timer.Start();
            if( wMaxNorm >= selInvTol )
                LDL( info, JFront, LDL_2D );
            else
                LDL( info, JFront, LDL_SELINV_2D );
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
                Output("Affine solve: ",timer.Stop()," secs");
        }
        catch(...)
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );

        if( ctrl.checkResiduals && ctrl.print )
        {
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
            dzError += dsAff;
            const Real dzErrorNrm2 = Nrm2( dzError );

            // TODO: dmuError
            // TODO: Also compute and print the residuals with regularization

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2) = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2) = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_h ||_2) = ",
                 dzErrorNrm2/(1+rhNrm2));
        }

        // Compute a centrality parameter using Mehrotra's formula
        // =======================================================
        Real alphaAffPri = pos_orth::MaxStep( s, dsAff, Real(1) );
        Real alphaAffDual = pos_orth::MaxStep( z, dzAff, Real(1) );
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
        Shift( rmu, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dsAff o dzAff
            // ---------------------
            // NOTE: dz is being used as a temporary
            dz = dzAff;
            DiagonalScale( LEFT, NORMAL, dsAff, dz );
            rmu += dz;
        }

        // Set up the new RHS
        // ------------------
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the new direction
        // -------------------------
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
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
        }
        ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );

        // Update the current estimates
        // ============================
        Real alphaPri = pos_orth::MaxStep( s, ds, 1/ctrl.maxStepRatio );
        Real alphaDual = pos_orth::MaxStep( z, dz, 1/ctrl.maxStepRatio );
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
            Output("iteration: ",iterTimer.Stop()," secs");
        if( alphaPri == Real(0) && alphaDual == Real(0) )
        {
            if( relError <= ctrl.minTol )
                break;
            else
                RuntimeError
                ("Could not achieve minimum tolerance of ",ctrl.minTol);
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
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, \
    const Matrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& G, \
    const ElementalMatrix<Real>& b, \
    const ElementalMatrix<Real>& c, \
    const ElementalMatrix<Real>& h, \
          ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
          ElementalMatrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, \
    const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
          Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, \
    const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
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
#include "El/macros/Instantiate.h"

} // namespace affine
} // namespace qp
} // namespace El
