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
namespace qp {
namespace direct {

// The following solves a pair of quadratic programs in "direct" conic form:
//
//   min (1/2) x^T Q x + c^T x
//   s.t. A x = b, x >= 0,
//
//   max (1/2) (A^T y - z + c)^T pinv(Q) (A^T y - z + c) - b^T y
//   s.t. A^T y - z + c in range(Q), z >= 0,
//
// as opposed to the more general "affine" conic form:
//
//   min (1/2) x^T Q x + c^T x
//   s.t. A x = b, G x + s = h, s >= 0,
//
//   max (1/2) (A^T y + G^T z + c)^T pinv(Q) (A^T y + G^T z + c) - b^T y - h^T z
//   s.t. A^T y + G^T z + c in range(Q), z >= 0
//
// using a Mehrotra Predictor-Corrector scheme.
//

template<typename Real>
void Mehrotra
( const Matrix<Real>& QPre,
  const Matrix<Real>& APre,
  const Matrix<Real>& bPre,
  const Matrix<Real>& cPre,
        Matrix<Real>& x,
        Matrix<Real>& y, 
        Matrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_CSE

    const bool stepLengthSigma = true;
    const bool standardShift = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;

    // Equilibrate the QP by diagonally scaling A
    auto Q = QPre;
    auto A = APre;
    auto b = bPre;
    auto c = cPre;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int degree = n;
    Matrix<Real> dRow, dCol;
    if( ctrl.outerEquil )
    {
        RuizEquil( A, dRow, dCol, ctrl.print );

        DiagonalSolve( LEFT, NORMAL, dRow, b );
        DiagonalSolve( LEFT, NORMAL, dCol, c );
        // TODO: Replace with SymmetricDiagonalSolve
        {
            DiagonalSolve( LEFT, NORMAL, dCol, Q );
            DiagonalSolve( RIGHT, NORMAL, dCol, Q );
        }
        if( ctrl.primalInit )
            DiagonalScale( LEFT, NORMAL, dCol, x );
        if( ctrl.dualInit )
        {
            DiagonalScale( LEFT, NORMAL, dRow, y );
            DiagonalSolve( LEFT, NORMAL, dCol, z );
        }
    }
    else
    {
        Ones( dRow, m, 1 );
        Ones( dCol, n, 1 );
    }

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    if( ctrl.print )
    {
        const Real QNrm1 = HermitianOneNorm( LOWER, Q );
        const Real ANrm1 = OneNorm( A );
        Output("|| Q ||_1 = ",QNrm1);
        Output("|| A ||_1 = ",ANrm1);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
    }

    Initialize
    ( Q, A, b, c, x, y, z, ctrl.primalInit, ctrl.dualInit, standardShift ); 

    Real relError = 1;
    Matrix<Real> J, d, 
                 rb,    rc,    rmu,
                 dxAff, dyAff, dzAff,
                 dx,    dy,    dz;
    Matrix<Real> dSub;
    Permutation p;
    Matrix<Real> dxError, dyError, dzError, prod;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = pos_orth::NumOutside( x );
        const Int zNumNonPos = pos_orth::NumOutside( z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure
        // ===========================
        const Real mu = Dot(x,z) / n;

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        Zeros( d, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj =  xTQx/2 + Dot(c,x);
        const Real dualObj = -xTQx/2 - Dot(b,y); 
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
        Hemv( LOWER,     Real(1), Q, x, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        rc -= z;
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // Now check the pieces
        // --------------------
        relError = Max(Max(objConv,rbConv),rcConv);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( x );
            const Real yNrm2 = Nrm2( y );
            const Real zNrm2 = Nrm2( z );
            Output
            ("iter ",numIts,":\n",Indent(),
             "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
             "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
             "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
             "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
             "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
             "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
             "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
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

        // r_mu := x o z
        // -------------
        rmu = z;
        DiagonalScale( LEFT, NORMAL, x, rmu );

        if( ctrl.system == FULL_KKT )
        {
            // Construct the KKT system
            // ------------------------
            KKT( Q, A, x, z, J );
            KKTRHS( rc, rb, rmu, z, d );

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
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandSolution( m, n, d, dxAff, dyAff, dzAff );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the KKT system
            // ------------------------
            AugmentedKKT( Q, A, x, z, J );
            AugmentedKKTRHS( x, rc, rb, rmu, d );

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
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandAugmentedSolution( x, z, rmu, d, dxAff, dyAff, dzAff );
        }
        else
            LogicError("Invalid KKT system choice");

        if( ctrl.checkResiduals && ctrl.print )
        {
            dxError = rb;
            Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
            Real dxErrorNrm2 = Nrm2( dxError );

            dyError = rc;
            Hemv( LOWER,     Real(1), Q, dxAff, Real(1), dyError );
            Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
            dyError -= dzAff;
            Real dyErrorNrm2 = Nrm2( dyError );

            Real rmuNrm2 = Nrm2( rmu );
            dzError = rmu;
            prod = dzAff;
            DiagonalScale( LEFT, NORMAL, x, prod );
            dzError += prod;
            prod = dxAff;
            DiagonalScale( LEFT, NORMAL, z, prod );
            dzError += prod;
            Real dzErrorNrm2 = Nrm2( dzError );

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2)  = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2)  = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_mu ||_2) = ",
             dzErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = pos_orth::MaxStep( x, dxAff, Real(1) );
        Real alphaAffDual = pos_orth::MaxStep( z, dzAff, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and dx are used as temporaries
        dx = x;
        dz = z;
        Axpy( alphaAffPri,  dxAff, dx );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(dx,dz) / degree;
        if( ctrl.print )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma = centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        rb *= 1-sigma;
        rc *= 1-sigma;
        Shift( rmu, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dxAff o dzAff
            // ---------------------
            // NOTE: dz is used as a temporary
            dz = dzAff;
            DiagonalScale( LEFT, NORMAL, dxAff, dz );
            rmu += dz;
        }

        if( ctrl.system == FULL_KKT )
        {
            // Construct the new KKT RHS
            // -------------------------
            KKTRHS( rc, rb, rmu, z, d );

            // Solve for the direction
            // -----------------------
            try { ldl::SolveAfter( J, dSub, p, d, false ); }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandSolution( m, n, d, dx, dy, dz );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the new KKT RHS
            // -------------------------
            AugmentedKKTRHS( x, rc, rb, rmu, d );

            // Solve for the direction
            // -----------------------
            try { ldl::SolveAfter( J, dSub, p, d, false ); }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }
        else
            LogicError("Invalid KKT system choice");
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = pos_orth::MaxStep( x, dx, 1/ctrl.maxStepRatio );
        Real alphaDual = pos_orth::MaxStep( z, dz, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  dx, x );
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
        DiagonalSolve( LEFT, NORMAL, dCol, x );
        DiagonalSolve( LEFT, NORMAL, dRow, y );
        DiagonalScale( LEFT, NORMAL, dCol, z );
    }
}

template<typename Real>
void Mehrotra
( const ElementalMatrix<Real>& QPre,
  const ElementalMatrix<Real>& APre, 
  const ElementalMatrix<Real>& bPre,
  const ElementalMatrix<Real>& cPre,
        ElementalMatrix<Real>& xPre,
        ElementalMatrix<Real>& yPre,
        ElementalMatrix<Real>& zPre,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_CSE

    const bool stepLengthSigma = true;
    const bool standardShift = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;

    const Grid& grid = APre.Grid();
    const int commRank = grid.Rank();

    // Ensure that the inputs have the appropriate read/write properties
    DistMatrix<Real> Q(grid), A(grid), b(grid), c(grid);
    Q.Align(0,0);
    A.Align(0,0);
    b.Align(0,0);
    c.Align(0,0);
    Q = QPre;
    A = APre;
    b = bPre;
    c = cPre;

    ElementalProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;

    // NOTE: x does not need to be a read proxy when !ctrl.primalInit
    DistMatrixReadWriteProxy<Real,Real,MC,MR>
      xProx( xPre, control ),
    // NOTE: {y,z} do not need to be read proxies when !ctrl.dualInit
      yProx( yPre, control ),
      zProx( zPre, control );
    auto& x = xProx.Get();
    auto& y = yProx.Get();
    auto& z = zProx.Get();

    // Equilibrate the QP by diagonally scaling A
    const Int m = A.Height();
    const Int n = A.Width();
    const Int degree = n;
    DistMatrix<Real,MC,STAR> dRow(grid);
    DistMatrix<Real,MR,STAR> dCol(grid);
    if( ctrl.outerEquil )
    {
        RuizEquil( A, dRow, dCol, ctrl.print );

        DiagonalSolve( LEFT, NORMAL, dRow, b );
        DiagonalSolve( LEFT, NORMAL, dCol, c );
        // TODO: Replace with SymmetricDiagonalSolve
        {
            DiagonalSolve( LEFT, NORMAL, dCol, Q );
            DiagonalSolve( RIGHT, NORMAL, dCol, Q );
        }
        if( ctrl.primalInit )
            DiagonalScale( LEFT, NORMAL, dCol, x );
        if( ctrl.dualInit )
        {
            DiagonalScale( LEFT, NORMAL, dRow, y );
            DiagonalSolve( LEFT, NORMAL, dCol, z );
        }
    }
    else
    {
        Ones( dRow, m, 1 );
        Ones( dCol, n, 1 );
    }

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    if( ctrl.print )
    {
        const Real QNrm1 = HermitianOneNorm( LOWER, Q );
        const Real ANrm1 = OneNorm( A );
        if( commRank == 0 )
        {
            Output("|| Q ||_1 = ",QNrm1);
            Output("|| A ||_1 = ",ANrm1);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
        }
    }

    Initialize
    ( Q, A, b, c, x, y, z, ctrl.primalInit, ctrl.dualInit, standardShift ); 

    Real relError = 1;
    DistMatrix<Real> 
        J(grid), d(grid), 
        rc(grid),    rb(grid),    rmu(grid), 
        dxAff(grid), dyAff(grid), dzAff(grid),
        dx(grid),    dy(grid),    dz(grid);
    dx.AlignWith( x );
    dz.AlignWith( x );
    dxAff.AlignWith( x );
    dzAff.AlignWith( x );
    rmu.AlignWith( x );
    DistMatrix<Real> dSub(grid);
    DistPermutation p(grid);
    DistMatrix<Real> dxError(grid), dyError(grid), dzError(grid), prod(grid);
    dzError.AlignWith( dz );
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = pos_orth::NumOutside( x );
        const Int zNumNonPos = pos_orth::NumOutside( z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure
        // ===========================
        const Real mu = Dot(x,z) / n;

        // Check for convergence
        // =====================
        // |primal - dual| / (1 + |primal|) <= tol ?
        // -----------------------------------------
        Zeros( d, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj =  xTQx/2 + Dot(c,x);
        const Real dualObj = -xTQx/2 - Dot(b,y); 
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
        Hemv( LOWER,     Real(1), Q, x, Real(1), rc );
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        rc -= z;
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // Now check the pieces
        // --------------------
        relError = Max(Max(objConv,rbConv),rcConv);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( x );
            const Real yNrm2 = Nrm2( y );
            const Real zNrm2 = Nrm2( z );
            if( commRank == 0 )
                Output
                ("iter ",numIts,":\n",Indent(),
                 "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
                 "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
                 "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
                 "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
                 "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
                 "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
                 "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
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

        // r_mu := x o z
        // -------------
        rmu = z;
        DiagonalScale( LEFT, NORMAL, x, rmu );

        if( ctrl.system == FULL_KKT )
        {
            // Construct the KKT system
            // ------------------------
            KKT( Q, A, x, z, J );
            KKTRHS( rc, rb, rmu, z, d );

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
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandSolution( m, n, d, dxAff, dyAff, dzAff );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the KKT system
            // ------------------------
            AugmentedKKT( Q, A, x, z, J );
            AugmentedKKTRHS( x, rc, rb, rmu, d );

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
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandAugmentedSolution( x, z, rmu, d, dxAff, dyAff, dzAff );
        }
        else
            LogicError("Invalid KKT system choice");

        if( ctrl.checkResiduals && ctrl.print )
        {
            dxError = rb;
            Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
            Real dxErrorNrm2 = Nrm2( dxError );

            dyError = rc;
            Hemv( LOWER,     Real(1), Q, dxAff, Real(1), dyError );
            Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
            dyError -= dzAff;
            Real dyErrorNrm2 = Nrm2( dyError );

            Real rmuNrm2 = Nrm2( rmu );
            dzError = rmu;
            prod = dzAff;
            DiagonalScale( LEFT, NORMAL, x, prod );
            dzError += prod;
            prod = dxAff;
            DiagonalScale( LEFT, NORMAL, z, prod );
            dzError += prod;
            Real dzErrorNrm2 = Nrm2( dzError );

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2)  = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2)  = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_mu ||_2) = ",
                 dzErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = pos_orth::MaxStep( x, dxAff, Real(1) );
        Real alphaAffDual = pos_orth::MaxStep( z, dzAff, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and dx are used as temporaries
        dx = x;
        dz = z;
        Axpy( alphaAffPri,  dxAff, dx );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(dx,dz) / degree;
        if( ctrl.print && commRank == 0 )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma = centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output("sigma=",sigma);

        // Compute the combined direction
        // ==============================
        rb *= 1-sigma;
        rc *= 1-sigma;
        Shift( rmu, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu := dxAff o dzAff
            // ---------------------
            // NOTE: dz is used as a temporary
            dz = dzAff;
            DiagonalScale( LEFT, NORMAL, dxAff, dz );
            rmu += dz;
        }

        if( ctrl.system == FULL_KKT )
        {
            // Construct the new KKT RHS
            // -------------------------
            KKTRHS( rc, rb, rmu, z, d );

            // Solve for the direction
            // -----------------------
            try { ldl::SolveAfter( J, dSub, p, d, false ); }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandSolution( m, n, d, dx, dy, dz );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Construct the new KKT RHS
            // -------------------------
            AugmentedKKTRHS( x, rc, rb, rmu, d );

            // Solve for the direction
            // -----------------------
            try { ldl::SolveAfter( J, dSub, p, d, false ); }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }
        else
            LogicError("Invalid KKT system choice");
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = pos_orth::MaxStep( x, dx, 1/ctrl.maxStepRatio );
        Real alphaDual = pos_orth::MaxStep( z, dz, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  dx, x );
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
        DiagonalSolve( LEFT, NORMAL, dCol, x );
        DiagonalSolve( LEFT, NORMAL, dRow, y );
        DiagonalScale( LEFT, NORMAL, dCol, z );
    }
}

template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& QPre,
  const SparseMatrix<Real>& APre, 
  const Matrix<Real>& bPre,
  const Matrix<Real>& cPre,
        Matrix<Real>& x,
        Matrix<Real>& y, 
        Matrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_CSE

    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const bool standardShift = true;

    // Equilibrate the QP by diagonally scaling A
    auto Q = QPre;
    auto A = APre;
    auto b = bPre;
    auto c = cPre;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int degree = n;
    Matrix<Real> dRow, dCol;
    if( ctrl.outerEquil )
    {
        RuizEquil( A, dRow, dCol, ctrl.print );

        DiagonalSolve( LEFT, NORMAL, dRow, b );
        DiagonalSolve( LEFT, NORMAL, dCol, c );
        // TODO: Replace with SymmetricDiagonalSolve
        {
            DiagonalSolve( LEFT, NORMAL, dCol, Q );
            DiagonalSolve( RIGHT, NORMAL, dCol, Q );
        }
        if( ctrl.primalInit )
            DiagonalScale( LEFT, NORMAL, dCol, x );
        if( ctrl.dualInit )
        {
            DiagonalScale( LEFT, NORMAL, dRow, y );
            DiagonalSolve( LEFT, NORMAL, dCol, z );
        }
    }
    else
    {
        Ones( dRow, m, 1 );
        Ones( dCol, n, 1 );
    }

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real twoNormEstQ = HermitianTwoNormEstimate( Q, ctrl.basisSize );
    const Real twoNormEstA = TwoNormEstimate( A, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstQ + twoNormEstA + 1;
    if( ctrl.print )
    {
        Output("|| Q ||_2 estimate: ",twoNormEstQ);
        Output("|| A ||_2 estimate: ",twoNormEstA);
        Output("|| b ||_2 = ",bNrm2);
        Output("|| c ||_2 = ",cNrm2);
    }

    vector<Int> map, invMap;
    ldl::NodeInfo info;
    ldl::Separator rootSep;
    // The initialization involves an augmented KKT system, and so we can
    // only reuse the factorization metadata if the this IPM is using the
    // augmented formulation
    // TODO: Add permanent regularization and cache J metadata
    if( ctrl.system == AUGMENTED_KKT )
    {
        Initialize
        ( Q, A, b, c, x, y, z, map, invMap, rootSep, info,
          ctrl.primalInit, ctrl.dualInit, standardShift, ctrl.solveCtrl );
    }  
    else
    {
        vector<Int> augMap, augInvMap;
        ldl::NodeInfo augInfo;
        ldl::Separator augRootSep;
        Initialize
        ( Q, A, b, c, x, y, z, augMap, augInvMap, augRootSep, augInfo,
          ctrl.primalInit, ctrl.dualInit, standardShift, ctrl.solveCtrl );
    }

    Matrix<Real> regTmp;
    if( ctrl.system == FULL_KKT )
    {
        regTmp.Resize( m+2*n, 1 );
        for( Int i=0; i<m+2*n; ++i )
        {
            if( i < n )        regTmp(i) =  ctrl.reg0Tmp*ctrl.reg0Tmp;
            else if( i < n+m ) regTmp(i) = -ctrl.reg1Tmp*ctrl.reg1Tmp;
            else               regTmp(i) = -ctrl.reg2Tmp*ctrl.reg2Tmp;
        }
    }
    else if( ctrl.system == AUGMENTED_KKT )
    {
        regTmp.Resize( n+m, 1 );
        for( Int i=0; i<n+m; ++i )
        {
            if( i < n ) regTmp(i) =  ctrl.reg0Tmp*ctrl.reg0Tmp;
            else        regTmp(i) = -ctrl.reg1Tmp*ctrl.reg1Tmp;
        }
    }
    regTmp *= origTwoNormEst;

    SparseMatrix<Real> J, JOrig;
    ldl::Front<Real> JFront;
    Matrix<Real> d, 
                 w,
                 rc,    rb,    rmu, 
                 dxAff, dyAff, dzAff,
                 dx,    dy,    dz;

    Real relError = 1;
    Matrix<Real> dInner;
    Matrix<Real> dxError, dyError, dzError, prod;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = pos_orth::NumOutside( x );
        const Int zNumNonPos = pos_orth::NumOutside( z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure and scaling point
        // =============================================
        const Real mu = Dot(x,z) / n;
        pos_orth::NesterovTodd( x, z, w );
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
        const Real dualObj = -xTQx/2 - Dot(b,y); 
        const Real objConv = Abs(primObj-dualObj) / (1+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        rb *= -1;
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (1+bNrm2);
        Axpy( -ctrl.reg1Perm*ctrl.reg1Perm, y, rb );
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        rc = c;
        Multiply( NORMAL,    Real(1), Q, x, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        rc -= z;
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        Axpy( ctrl.reg0Perm*ctrl.reg0Perm, x, rc );
        // Now check the pieces
        // --------------------
        relError = Max(Max(objConv,rbConv),rcConv);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( x );
            const Real yNrm2 = Nrm2( y );
            const Real zNrm2 = Nrm2( z );
            Output
            ("iter ",numIts,":\n",Indent(),
             "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
             "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
             "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
             "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
             "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
             "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
             "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
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

        // r_mu := x o z
        // -------------
        rmu = z; 
        DiagonalScale( LEFT, NORMAL, x, rmu );

        if( ctrl.system == FULL_KKT || ctrl.system == AUGMENTED_KKT )
        {
            // Form the KKT system
            // -------------------
            if( ctrl.system == FULL_KKT )
            {
                KKT
                ( Q, A, ctrl.reg0Perm, ctrl.reg1Perm, ctrl.reg2Perm, x, z,
                  JOrig, false );
                KKTRHS( rc, rb, rmu, z, d );
            }
            else
            {
                AugmentedKKT
                ( Q, A, ctrl.reg0Perm, ctrl.reg1Perm, x, z, JOrig, false );
                // TODO: Incorporate ctrl.reg2Perm?
                AugmentedKKTRHS( x, rc, rb, rmu, d );
            }

            // Solve for the direction
            // -----------------------
            try
            {
                J = JOrig;

                UpdateDiagonal( J, Real(1), regTmp );
                if( wMaxNorm >= ctrl.ruizEquilTol )
                    SymmetricRuizEquil
                    ( J, dInner, ctrl.ruizMaxIter, ctrl.print );
                else if( wMaxNorm >= ctrl.diagEquilTol )
                    SymmetricDiagonalEquil( J, dInner, ctrl.print );
                else
                    Ones( dInner, J.Height(), 1 );

                if( numIts == 0 && 
                    (ctrl.system == FULL_KKT || 
                     (ctrl.primalInit && ctrl.dualInit) ) )
                {
                    NestedDissection( J.LockedGraph(), map, rootSep, info );
                    InvertMap( map, invMap );
                }
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
            if( ctrl.system == FULL_KKT )
                ExpandSolution( m, n, d, dxAff, dyAff, dzAff );
            else
                ExpandAugmentedSolution( x, z, rmu, d, dxAff, dyAff, dzAff );
        }
        else
            LogicError("Invalid KKT system choice");

        if( ctrl.checkResiduals && ctrl.print )
        {
            dxError = rb;
            Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
            Real dxErrorNrm2 = Nrm2( dxError );

            dyError = rc;
            Multiply( NORMAL,    Real(1), Q, dxAff, Real(1), dyError );
            Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
            dyError -= dzAff;
            Real dyErrorNrm2 = Nrm2( dyError );

            Real rmuNrm2 = Nrm2( rmu );
            dzError = rmu;
            prod = dzAff;
            DiagonalScale( LEFT, NORMAL, x, prod );
            dzError += prod;
            prod = dxAff;
            DiagonalScale( LEFT, NORMAL, z, prod );
            dzError += prod;
            Real dzErrorNrm2 = Nrm2( dzError );

            Output
            ("|| dxError ||_2 / (1 + || r_b ||_2)  = ",
             dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
             "|| dyError ||_2 / (1 + || r_c ||_2)  = ",
             dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
             "|| dzError ||_2 / (1 + || r_mu ||_2) = ",
             dzErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = pos_orth::MaxStep( x, dxAff, Real(1) );
        Real alphaAffDual = pos_orth::MaxStep( z, dzAff, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and dx are used as temporaries
        dx = x;
        dz = z;
        Axpy( alphaAffPri,  dxAff, dx );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(dx,dz) / degree;
        if( ctrl.print )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma = centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output("sigma=",sigma);

        // Compute the combined direction
        // ==============================
        rb *= 1-sigma;
        rc *= 1-sigma;
        Shift( rmu, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dxAff o dzAff
            // ---------------------
            // NOTE: dz is being used as a temporary
            dz = dzAff;
            DiagonalScale( LEFT, NORMAL, dxAff, dz );
            rmu += dz;
        }

        if( ctrl.system == FULL_KKT )
        {
            // Form the new KKT RHS
            // --------------------
            KKTRHS( rc, rb, rmu, z, d );
            // Solve for the direction
            // -----------------------
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
            ExpandSolution( m, n, d, dx, dy, dz );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Form the new KKT RHS
            // --------------------
            AugmentedKKTRHS( x, rc, rb, rmu, d );
            // Solve for the direction
            // -----------------------
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
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }
        else
            LogicError("Invalid KKT system choice");
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = pos_orth::MaxStep( x, dx, 1/ctrl.maxStepRatio );
        Real alphaDual = pos_orth::MaxStep( z, dz, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  dx, x );
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
        DiagonalSolve( LEFT, NORMAL, dCol, x );
        DiagonalSolve( LEFT, NORMAL, dRow, y );
        DiagonalScale( LEFT, NORMAL, dCol, z );
    }
}

template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& QPre,
  const DistSparseMatrix<Real>& APre, 
  const DistMultiVec<Real>& bPre,
  const DistMultiVec<Real>& cPre,
        DistMultiVec<Real>& x,
        DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_CSE

    const bool stepLengthSigma = true;
    function<Real(Real,Real,Real,Real)> centralityRule;
    if( stepLengthSigma )
        centralityRule = StepLengthCentrality<Real>;
    else
        centralityRule = MehrotraCentrality<Real>;
    const bool standardShift = true;

    mpi::Comm comm = APre.Comm();
    const int commRank = mpi::Rank(comm);
    Timer timer;

    // Equilibrate the QP by diagonally scaling A
    auto Q = QPre;
    auto A = APre;
    auto b = bPre;
    auto c = cPre;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int degree = n;
    DistMultiVec<Real> dRow(comm), dCol(comm);
    if( ctrl.outerEquil )
    {
        if( commRank == 0 && ctrl.time )
            timer.Start();
        RuizEquil( A, dRow, dCol, ctrl.print );
        if( commRank == 0 && ctrl.time )
            Output("RuizEquil: ",timer.Stop()," secs");

        DiagonalSolve( LEFT, NORMAL, dRow, b );
        DiagonalSolve( LEFT, NORMAL, dCol, c );
        // TODO: Replace with SymmetricDiagonalSolve
        {
            DiagonalSolve( LEFT, NORMAL, dCol, Q );
            DiagonalSolve( RIGHT, NORMAL, dCol, Q );
        }
        if( ctrl.primalInit )
            DiagonalScale( LEFT, NORMAL, dCol, x );
        if( ctrl.dualInit )
        {
            DiagonalScale( LEFT, NORMAL, dRow, y );
            DiagonalSolve( LEFT, NORMAL, dCol, z );
        }
    }
    else
    {
        Ones( dRow, m, 1 );
        Ones( dCol, n, 1 );
    }

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real twoNormEstQ = HermitianTwoNormEstimate( Q, ctrl.basisSize );
    const Real twoNormEstA = TwoNormEstimate( A, ctrl.basisSize );
    const Real origTwoNormEst = twoNormEstQ + twoNormEstA + 1;
    if( ctrl.print )
    {
        const double imbalanceQ = Q.Imbalance();
        const double imbalanceA = A.Imbalance();
        if( commRank == 0 )
        {
            Output("|| Q ||_2 estimate: ",twoNormEstQ);
            Output("|| A ||_2 estimate: ",twoNormEstA);
            Output("|| b ||_2 = ",bNrm2);
            Output("|| c ||_2 = ",cNrm2);
            Output("Imbalance factor of Q: ",imbalanceQ);
            Output("Imbalance factor of A: ",imbalanceA);
        }
    }

    DistMap map, invMap;
    ldl::DistNodeInfo info;
    ldl::DistSeparator rootSep;
    vector<Int> mappedSources, mappedTargets, colOffs;
    // The initialization involves an augmented KKT system, and so we can
    // only reuse the factorization metadata if the this IPM is using the
    // augmented formulation
    // TODO: Add permanent regularization and cache J metadata
    if( commRank == 0 && ctrl.time )
        timer.Start();
    if( ctrl.system == AUGMENTED_KKT )
    {
        Initialize
        ( Q, A, b, c, x, y, z, map, invMap, rootSep, info,
          mappedSources, mappedTargets, colOffs,
          ctrl.primalInit, ctrl.dualInit, standardShift, ctrl.solveCtrl ); 
    }  
    else
    {
        DistMap augMap, augInvMap;
        ldl::DistNodeInfo augInfo;
        ldl::DistSeparator augRootSep;
        vector<Int> augMappedSources, augMappedTargets, augColOffs;
        Initialize
        ( Q, A, b, c, x, y, z, augMap, augInvMap, augRootSep, augInfo,
          augMappedSources, augMappedTargets, augColOffs,
          ctrl.primalInit, ctrl.dualInit, standardShift, ctrl.solveCtrl );
    }
    if( commRank == 0 && ctrl.time )
        Output("Init: ",timer.Stop()," secs");

    DistMultiVec<Real> regTmp(comm);
    if( ctrl.system == FULL_KKT )
    {
        regTmp.Resize( m+2*n, 1 );
        for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
        {
            const Int i = regTmp.GlobalRow(iLoc);
            if( i < n )
              regTmp.SetLocal( iLoc, 0,  ctrl.reg0Tmp*ctrl.reg0Tmp );
            else if( i < n+m )
              regTmp.SetLocal( iLoc, 0, -ctrl.reg1Tmp*ctrl.reg1Tmp );
            else
              regTmp.SetLocal( iLoc, 0, -ctrl.reg2Tmp*ctrl.reg2Tmp );
        }
    }
    else if( ctrl.system == AUGMENTED_KKT )
    {
        regTmp.Resize( n+m, 1 );
        for( Int iLoc=0; iLoc<regTmp.LocalHeight(); ++iLoc )
        {
            const Int i = regTmp.GlobalRow(iLoc);
            if( i < n ) regTmp.SetLocal( iLoc, 0,  ctrl.reg0Tmp*ctrl.reg0Tmp );
            else        regTmp.SetLocal( iLoc, 0, -ctrl.reg1Tmp*ctrl.reg1Tmp );
        }
    }
    regTmp *= origTwoNormEst;

    DistGraphMultMeta metaOrig, meta;
    DistSparseMatrix<Real> J(comm), JOrig(comm);
    ldl::DistFront<Real> JFront;
    DistMultiVec<Real> d(comm), 
                       w(comm),
                       rc(comm),    rb(comm),    rmu(comm), 
                       dxAff(comm), dyAff(comm), dzAff(comm),
                       dx(comm),    dy(comm),    dz(comm);

    Real relError = 1;
    DistMultiVec<Real> dInner(comm);
    DistMultiVec<Real> dxError(comm), dyError(comm), dzError(comm), prod(comm);
    ldl::DistMultiVecNodeMeta dmvMeta;
    const Int indent = PushIndent();
    for( Int numIts=0; numIts<=ctrl.maxIts; ++numIts )
    {
        // Ensure that x and z are in the cone
        // ===================================
        const Int xNumNonPos = pos_orth::NumOutside( x );
        const Int zNumNonPos = pos_orth::NumOutside( z );
        if( xNumNonPos > 0 || zNumNonPos > 0 )
            LogicError
            (xNumNonPos," entries of x were nonpositive and ",
             zNumNonPos," entries of z were nonpositive");

        // Compute the duality measure and scaling point
        // =============================================
        const Real mu = Dot(x,z) / degree;
        pos_orth::NesterovTodd( x, z, w );
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
        const Real dualObj = -xTQx/2 - Dot(b,y); 
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
        Multiply( NORMAL,    Real(1), Q, x, Real(1), rc );
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        rc -= z;
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (1+cNrm2);
        // Now check the pieces
        // --------------------
        relError = Max(Max(objConv,rbConv),rcConv);
        if( ctrl.print )
        {
            const Real xNrm2 = Nrm2( x );
            const Real yNrm2 = Nrm2( y );
            const Real zNrm2 = Nrm2( z );
            if( commRank == 0 )
                Output
                ("iter ",numIts,":\n",Indent(),
                 "  ||  x  ||_2 = ",xNrm2,"\n",Indent(),
                 "  ||  y  ||_2 = ",yNrm2,"\n",Indent(),
                 "  ||  z  ||_2 = ",zNrm2,"\n",Indent(),
                 "  || r_b ||_2 = ",rbNrm2,"\n",Indent(),
                 "  || r_c ||_2 = ",rcNrm2,"\n",Indent(),
                 "  || r_b ||_2 / (1 + || b ||_2) = ",rbConv,"\n",Indent(),
                 "  || r_c ||_2 / (1 + || c ||_2) = ",rcConv,"\n",Indent(),
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

        // r_mu := x o z
        // -------------
        rmu = z; 
        DiagonalScale( LEFT, NORMAL, x, rmu );

        if( ctrl.system == FULL_KKT || ctrl.system == AUGMENTED_KKT )
        {
            // Form the KKT system
            // -------------------
            if( ctrl.system == FULL_KKT )
            {
                KKT
                ( Q, A, ctrl.reg0Perm, ctrl.reg1Perm, ctrl.reg2Perm, x, z,
                  JOrig, false );
                KKTRHS( rc, rb, rmu, z, d );
            }
            else
            {
                AugmentedKKT
                ( Q, A, ctrl.reg0Perm, ctrl.reg1Perm, x, z, JOrig, false );
                AugmentedKKTRHS( x, rc, rb, rmu, d );
            }

            // Solve for the direction
            // -----------------------
            try
            {
                // Cache the metadata for the finalized JOrig
                if( numIts == 0 )
                    metaOrig = JOrig.InitializeMultMeta();
                else
                    JOrig.LockedDistGraph().multMeta = metaOrig;
                J = JOrig;

                UpdateDiagonal( J, Real(1), regTmp );
                // Cache the metadata for the finalized J
                if( numIts == 0 )
                {
                    if( ctrl.print )
                    {
                        const double imbalanceJ = J.Imbalance();
                        if( commRank == 0 )
                            Output("Imbalance factor of J: ",imbalanceJ);
                    }

                    meta = J.InitializeMultMeta();
                    if( ctrl.system == FULL_KKT || 
                        (ctrl.primalInit && ctrl.dualInit) )
                    {
                        if( commRank == 0 && ctrl.time )
                            timer.Start();
                        NestedDissection
                        ( J.LockedDistGraph(), map, rootSep, info );
                        if( commRank == 0 && ctrl.time )
                            Output("ND: ",timer.Stop()," secs");
                        InvertMap( map, invMap );
                    }
                }
                else
                    J.LockedDistGraph().multMeta = meta;

                if( commRank == 0 && ctrl.time )
                    timer.Start();
                if( wMaxNorm >= ctrl.ruizEquilTol )
                    SymmetricRuizEquil
                    ( J, dInner, ctrl.ruizMaxIter, ctrl.print );
                else if( wMaxNorm >= ctrl.diagEquilTol )
                    SymmetricDiagonalEquil( J, dInner, ctrl.print );
                else
                    Ones( dInner, J.Height(), 1 );
                if( commRank == 0 && ctrl.time )
                    Output("Equilibration: ",timer.Stop()," secs");

                JFront.Pull
                ( J, map, rootSep, info, 
                  mappedSources, mappedTargets, colOffs );

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
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            if( ctrl.system == FULL_KKT )
                ExpandSolution( m, n, d, dxAff, dyAff, dzAff );
            else
                ExpandAugmentedSolution( x, z, rmu, d, dxAff, dyAff, dzAff );
        }
        else
            LogicError("Invalid KKT system choice");

        if( ctrl.checkResiduals && ctrl.print )
        {
            dxError = rb;
            Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
            Real dxErrorNrm2 = Nrm2( dxError );

            dyError = rc;
            Multiply( NORMAL,    Real(1), Q, dxAff, Real(1), dyError );
            Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
            dyError -= dzAff;
            Real dyErrorNrm2 = Nrm2( dyError );

            Real rmuNrm2 = Nrm2( rmu );
            dzError = rmu;
            prod = dzAff;
            DiagonalScale( LEFT, NORMAL, x, prod );
            dzError += prod;
            prod = dxAff;
            DiagonalScale( LEFT, NORMAL, z, prod );
            dzError += prod;
            Real dzErrorNrm2 = Nrm2( dzError );

            if( commRank == 0 )
                Output
                ("|| dxError ||_2 / (1 + || r_b ||_2)  = ",
                 dxErrorNrm2/(1+rbNrm2),"\n",Indent(),
                 "|| dyError ||_2 / (1 + || r_c ||_2)  = ",
                 dyErrorNrm2/(1+rcNrm2),"\n",Indent(),
                 "|| dzError ||_2 / (1 + || r_mu ||_2) = ",
                 dzErrorNrm2/(1+rmuNrm2));
        }

        // Compute a centrality parameter
        // ==============================
        Real alphaAffPri = pos_orth::MaxStep( x, dxAff, Real(1) );
        Real alphaAffDual = pos_orth::MaxStep( z, dzAff, Real(1) );
        if( ctrl.forceSameStep )
            alphaAffPri = alphaAffDual = Min(alphaAffPri,alphaAffDual);
        if( ctrl.print && commRank == 0 )
            Output
            ("alphaAffPri = ",alphaAffPri,", alphaAffDual = ",alphaAffDual);
        // NOTE: dz and dx are used as temporaries
        dx = x;
        dz = z;
        Axpy( alphaAffPri,  dxAff, dx );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(dx,dz) / degree;
        if( ctrl.print )
            Output("muAff = ",muAff,", mu = ",mu);
        const Real sigma = centralityRule(mu,muAff,alphaAffPri,alphaAffDual);
        if( ctrl.print )
            Output("sigma=",sigma);

        // Solve for the combined direction
        // ================================
        rb *= 1-sigma;
        rc *= 1-sigma;
        Shift( rmu, -sigma*mu );
        if( ctrl.mehrotra )
        {
            // r_mu += dxAff o dzAff
            // ---------------------
            // NOTE: dz is used as a temporary
            dz = dzAff;
            DiagonalScale( LEFT, NORMAL, dxAff, dz );
            rmu += dz;
        }

        if( ctrl.system == FULL_KKT )
        {
            // Form the KKT system
            // -------------------
            KKTRHS( rc, rb, rmu, z, d );
            // Solve for the direction
            // -----------------------
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
                    Output("Corrector: ",timer.Stop()," secs");
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandSolution( m, n, d, dx, dy, dz );
        }
        else if( ctrl.system == AUGMENTED_KKT )
        {
            // Form the KKT system
            // -------------------
            AugmentedKKTRHS( x, rc, rb, rmu, d );
            // Solve for the direction
            // -----------------------
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
                    Output("Corrector: ",timer.Stop()," secs");
            }
            catch(...)
            {
                if( relError <= ctrl.minTol )
                    break;
                else
                    RuntimeError
                    ("Could not achieve minimum tolerance of ",ctrl.minTol);
            }
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }
        else
            LogicError("Invalid KKT system choice");
        // TODO: Residual checks

        // Update the current estimates
        // ============================
        Real alphaPri = pos_orth::MaxStep( x, dx, 1/ctrl.maxStepRatio );
        Real alphaDual = pos_orth::MaxStep( z, dz, 1/ctrl.maxStepRatio );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.forceSameStep )
            alphaPri = alphaDual = Min(alphaPri,alphaDual);
        if( ctrl.print && commRank == 0 )
            Output("alphaPri = ",alphaPri,", alphaDual = ",alphaDual);
        Axpy( alphaPri,  dx, x );
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
        DiagonalSolve( LEFT, NORMAL, dCol, x );
        DiagonalSolve( LEFT, NORMAL, dRow, y );
        DiagonalScale( LEFT, NORMAL, dCol, z );
    }
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const Matrix<Real>& Q, \
    const Matrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const ElementalMatrix<Real>& Q, \
    const ElementalMatrix<Real>& A, \
    const ElementalMatrix<Real>& b, \
    const ElementalMatrix<Real>& c, \
          ElementalMatrix<Real>& x, \
          ElementalMatrix<Real>& y, \
          ElementalMatrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& Q, \
    const SparseMatrix<Real>& A, \
    const Matrix<Real>& b, \
    const Matrix<Real>& c, \
          Matrix<Real>& x, \
          Matrix<Real>& y, \
          Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& Q, \
    const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, \
    const DistMultiVec<Real>& c, \
          DistMultiVec<Real>& x, \
          DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#define EL_ENABLE_DOUBLEDOUBLE
#define EL_ENABLE_QUADDOUBLE
#define EL_ENABLE_QUAD
#define EL_ENABLE_BIGFLOAT
#include <El/macros/Instantiate.h>

} // namespace direct
} // namespace qp
} // namespace El
