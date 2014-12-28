/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "./util.hpp"

namespace El {
namespace qp {
namespace primal {

// TODO: Experiment with running a few iterations of IPF until the residuals
//       are sufficiently low before switching to Mehrotra's scheme

// The following solves a quadratic program in "primal" conic form:
//
//   min 1/2 x^T Q x + c^T x
//   s.t. A x = b, x >= 0,
//
// as opposed to the more general "dual" conic form:
//
//   min 1/2 x^T Q x + c^T x
//   s.t. A x = b, G x + s = h, s >= 0,
//
// using a Mehrotra Predictor-Corrector scheme.
//

// TODO: Incorporate/extend recent lp::primal modifications

template<typename Real>
void Mehrotra
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c,
  Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::primal::Mehrotra"))    
    const Int m = A.Height();
    const Int n = A.Width();

    const Real bNrm2 = Nrm2( b );

    Matrix<Real> J, d, 
                 rmu,   rb,    rc, 
                 dxAff, dyAff, dzAff,
                 dx,    dy,    dz;
    Matrix<Real> dSub;
    Matrix<Int> p;
#ifndef EL_RELEASE
    Matrix<Real> dxError, dyError, dzError;
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
        // |c^T x + x^T Q x - b^T y| / (1 + |c^T x + 1/2 x^T Q x|) <= tol ?
        // ----------------------------------------------------------------
        Hemv( LOWER, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj = Dot(c,x) + xTQx/2;
        const Real dualObj = Dot(b,y) - xTQx/2; 
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        Scale( Real(-1), rb );
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c + Q x ||_2) <= tol ?
        // --------------------------------------------
        // NOTE: d currently contains Q x
        Axpy( Real(1), c, d );
        const Real objGradNrm2 = Nrm2( d );
        rc = d;
        Scale( Real(-1), rc );
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Axpy( Real(1), z, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+objGradNrm2);
        // Now check the pieces
        // --------------------
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;
        else if( ctrl.print )
            std::cout << " iter " << numIts << ":\n"
                      << "  |primObj - dualObj| / (1 + |primObj|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)         = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c + Q x ||_2)   = "
                      << rcConv << std::endl;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // r_mu := X Z e
        // =============
        rmu.Resize( n, 1 );
        for( Int i=0; i<n; ++i )
            rmu.Set( i, 0, x.Get(i,0)*z.Get(i,0) );

        // Compute the affine search direction
        // ===================================
        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // -----------------------------
            KKT( Q, A, x, z, J );
            KKTRHS( rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            LU( J, p );
            lu::SolveAfter( NORMAL, J, p, d );
            ExpandKKTSolution( m, n, d, dxAff, dyAff, dzAff );
        }
        else // ctrl.system == AUGMENTED_KKT
        {
            // Construct the "augmented" KKT system
            // ------------------------------------
            AugmentedKKT( Q, A, x, z, J );
            AugmentedKKTRHS( x, rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            LDL( J, dSub, p, false );
            ldl::SolveAfter( J, dSub, p, d, false );
            ExpandAugmentedSolution( x, z, rmu, d, dxAff, dyAff, dzAff );
        }

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        Real rmuNrm2 = Nrm2( rmu ); 
        dzError = rmu;
        for( Int i=0; i<n; ++i )
        {
            const Real xi = x.Get(i,0);
            const Real zi = z.Get(i,0);
            const Real dxi = dxAff.Get(i,0);
            const Real dzi = dzAff.Get(i,0);
            dzError.Update( i, 0, xi*dzi + zi*dxi );
        }
        Real dzErrorNrm2 = Nrm2( dzError );

        dyError = dzAff;
        Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Hemv( LOWER, Real(-1), Q, dx, Real(1), dyError );
        Axpy( Real(1), rc, dyError );
        Real dyErrorNrm2 = Nrm2( dyError );

        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        Real dxErrorNrm2 = Nrm2( dxError );

        if( ctrl.print )
            std::cout << "  || dxAffError ||_2 / (1 + || r_b ||_2)  = " 
                      << dxErrorNrm2/(1+rbNrm2) << "\n"
                      << "  || dyAffError ||_2 / (1 + || r_c ||_2)  = " 
                      << dyErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dzAffError ||_2 / (1 + || r_mu ||_2) = " 
                      << dzErrorNrm2/(1+rmuNrm2) << std::endl;
#endif

        // Compute the maximum affine [0,1]-step which preserves positivity
        // ================================================================
        Real alphaAffPri = 1; 
        for( Int i=0; i<n; ++i )
        {
            const Real xi = x.Get(i,0);
            const Real dxAffi = dxAff.Get(i,0);
            if( dxAffi < Real(0) )
                alphaAffPri = Min(alphaAffPri,-xi/dxAffi);
        }
        Real alphaAffDual = 1; 
        for( Int i=0; i<n; ++i )
        {
            const Real zi = z.Get(i,0);
            const Real dzAffi = dzAff.Get(i,0);
            if( dzAffi < Real(0) )
                alphaAffDual = Min(alphaAffDual,-zi/dzAffi);
        }
        if( ctrl.print )
            std::cout << "  alphaAffPri = " << alphaAffPri 
                      << ", alphaAffDual = " << alphaAffDual << std::endl;

        // Compute what the new duality measure would become
        // =================================================
        const Real mu = Dot(x,z) / n;
        // NOTE: dz and dx are used as temporaries
        dx = x;
        dz = z;
        Axpy( alphaAffPri,  dxAff, dx );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(dx,dz) / n;

        // Compute a centrality parameter using Mehotra's formula
        // ======================================================
        // TODO: Allow the user to override this function
        const Real sigma = Pow(muAff/mu,Real(3)); 
        if( ctrl.print )
            std::cout << "  muAff = " << muAff 
                      << ", mu = " << mu 
                      << ", sigma = " << sigma << std::endl;

        // Solve for the centering-corrector 
        // =================================
        Zeros( rc, n, 1 );
        Zeros( rb, m, 1 );
        for( Int i=0; i<n; ++i )
            rmu.Set( i, 0, dxAff.Get(i,0)*dzAff.Get(i,0) - sigma*mu );
        if( ctrl.system == FULL_KKT )
        {
            // Construct the new full KKT RHS
            // ------------------------------
            KKTRHS( rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            lu::SolveAfter( NORMAL, J, p, d );
            ExpandKKTSolution( m, n, d, dx, dy, dz );
        }
        else // ctrl.system == AUGMENTED_KKT
        {
            // Construct the new "augmented" KKT RHS
            // -------------------------------------
            AugmentedKKTRHS( x, rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            ldl::SolveAfter( J, dSub, p, d, false );
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }

        // TODO: Residual checks for center-corrector

        // Add in the affine search direction
        // ==================================
        Axpy( Real(1), dxAff, dx );
        Axpy( Real(1), dyAff, dy );
        Axpy( Real(1), dzAff, dz );

        // Compute the max positive [0,1/maxStepRatio] step length
        // =======================================================
        Real alphaPri = 1/ctrl.maxStepRatio;
        for( Int i=0; i<n; ++i )
        {
            const Real xi = x.Get(i,0);
            const Real dxi = dx.Get(i,0);
            if( dxi < Real(0) )
                alphaPri = Min(alphaPri,-xi/dxi);
        }
        Real alphaDual = 1/ctrl.maxStepRatio;
        for( Int i=0; i<n; ++i )
        {
            const Real zi = z.Get(i,0);
            const Real dzi = dz.Get(i,0);
            if( dzi < Real(0) )
                alphaDual = Min(alphaDual,-zi/dzi);
        }
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.print )
            std::cout << "  alphaPri = " << alphaPri 
                      << ", alphaDual = " << alphaDual << std::endl;

        // Update the current estimates
        // ============================
        Axpy( alphaPri,  dx, x );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z ); 
    }
}

template<typename Real>
void Mehrotra
( const AbstractDistMatrix<Real>& QPre, const AbstractDistMatrix<Real>& APre, 
  const AbstractDistMatrix<Real>& b,    const AbstractDistMatrix<Real>& c,
  AbstractDistMatrix<Real>& xPre, AbstractDistMatrix<Real>& y, 
  AbstractDistMatrix<Real>& zPre,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::primal::Mehrotra"))

    ProxyCtrl control;
    control.colConstrain = true;
    control.rowConstrain = true;
    control.colAlign = 0;
    control.rowAlign = 0;
    auto QPtr = ReadProxy<Real,MC,MR>(&QPre,control);      auto& Q = *QPtr;
    auto APtr = ReadProxy<Real,MC,MR>(&APre,control);      auto& A = *APtr;
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,control); auto& x = *xPtr;
    auto zPtr = ReadWriteProxy<Real,MC,MR>(&zPre,control); auto& z = *zPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const Int commRank = A.Grid().Rank();

    const Real bNrm2 = Nrm2( b );

    DistMatrix<Real> 
        J(grid), d(grid), 
        rmu(grid),   rb(grid),    rc(grid), 
        dxAff(grid), dyAff(grid), dzAff(grid),
        dx(grid),    dy(grid),    dz(grid);
    dx.AlignWith( x );
    dz.AlignWith( x );
    dxAff.AlignWith( x );
    dzAff.AlignWith( x );
    rmu.AlignWith( x );
    DistMatrix<Real> dSub(grid);
    DistMatrix<Int> p(grid);
#ifndef EL_RELEASE
    DistMatrix<Real> dxError(grid), dyError(grid), dzError(grid);
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
        // |c^T x + x^T Q x - b^T y| / (1 + |c^T x + 1/2 x^T Q x|) <= tol ?
        // ----------------------------------------------------------------
        Zeros( d, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj = Dot(c,x) + xTQx/2;
        const Real dualObj = Dot(b,y) - xTQx/2; 
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        Scale( Real(-1), rb );
        Gemv( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c + Q x ||_2) <= tol ?
        // --------------------------------------------
        // NOTE: d currently contains Q x
        Axpy( Real(1), c, d );
        const Real objGradNrm2 = Nrm2( d );
        rc = d;
        Scale( Real(-1), rc );
        Gemv( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Axpy( Real(1), z, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+objGradNrm2);
        // Now check the pieces
        // --------------------
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;
        else if( ctrl.print && commRank == 0 )
            std::cout << " iter " << numIts << ":\n"
                      << "  |primObj - dualObj| / (1 + |primObj|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)         = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c + Q x ||_2)   = "
                      << rcConv << std::endl;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // r_mu := X Z e
        // =============
        // TODO: Find a more convenient syntax for expressing this operation
        rmu.Resize( n, 1 );
        if( rmu.IsLocalCol(0) )
            for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
                rmu.SetLocal( iLoc, 0, x.GetLocal(iLoc,0)*z.GetLocal(iLoc,0) );

        // Compute the affine search direction
        // ===================================
        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // -----------------------------
            KKT( Q, A, x, z, J );
            KKTRHS( rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            LU( J, p );
            lu::SolveAfter( NORMAL, J, p, d );
            ExpandKKTSolution( m, n, d, dxAff, dyAff, dzAff );
        }
        else // ctrl.system == AUGMENTED_KKT
        {
            // Construct the "augmented" KKT system
            // ------------------------------------
            AugmentedKKT( Q, A, x, z, J );
            AugmentedKKTRHS( x, rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            LDL( J, dSub, p, false );
            ldl::SolveAfter( J, dSub, p, d, false );
            ExpandAugmentedSolution( x, z, rmu, d, dxAff, dyAff, dzAff );
        }

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        Real rmuNrm2 = Nrm2( rmu ); 
        // TODO: Find a more convenient syntax for expressing this operation
        dzError = rmu;
        if( dzError.IsLocalCol(0) )
        {
            for( Int iLoc=0; iLoc<dzError.LocalHeight(); ++iLoc )
            {
                const Real xi = x.GetLocal(iLoc,0);
                const Real zi = z.GetLocal(iLoc,0);
                const Real dxi = dxAff.GetLocal(iLoc,0);
                const Real dzi = dzAff.GetLocal(iLoc,0);
                dzError.UpdateLocal( iLoc, 0, xi*dzi + zi*dxi );
            }
        }
        Real dzErrorNrm2 = Nrm2( dzError );

        dyError = dzAff;
        Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Hemv( LOWER, Real(-1), Q, dxAff, Real(1), dyError );
        Axpy( Real(1), rc, dyError );
        Real dyErrorNrm2 = Nrm2( dyError );

        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        Real dxErrorNrm2 = Nrm2( dxError );

        if( ctrl.print && commRank == 0 )
            std::cout << "  || dxAffError ||_2 / (1 + || r_b ||_2)  = " 
                      << dxErrorNrm2/(1+rbNrm2) << "\n"
                      << "  || dyAffError ||_2 / (1 + || r_c ||_2)  = " 
                      << dyErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dzAffError ||_2 / (1 + || r_mu ||_2) = " 
                      << dzErrorNrm2/(1+rmuNrm2) << std::endl;
#endif

        // Compute the maximum affine [0,1]-step which preserves positivity
        // ================================================================
        Real alphaAffPri = 1; 
        // TODO: Find a more convenient way to represent this operation
        if( x.IsLocalCol(0) )
        {
            for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
            {
                const Real xi = x.GetLocal(iLoc,0);
                const Real dxAffi = dxAff.GetLocal(iLoc,0);
                if( dxAffi < Real(0) )
                    alphaAffPri = Min(alphaAffPri,-xi/dxAffi);
            }
        }
        alphaAffPri = mpi::AllReduce( alphaAffPri, mpi::MIN, x.DistComm() );
        Real alphaAffDual = 1; 
        // TODO: Find a more convenient way to represent this operation
        if( z.IsLocalCol(0) )
        {
            for( Int iLoc=0; iLoc<z.LocalHeight(); ++iLoc )
            {
                const Real zi = z.GetLocal(iLoc,0);
                const Real dzAffi = dzAff.GetLocal(iLoc,0);
                if( dzAffi < Real(0) )
                    alphaAffDual = Min(alphaAffDual,-zi/dzAffi);
            }
        }
        alphaAffDual = mpi::AllReduce( alphaAffDual, mpi::MIN, z.DistComm() );
        if( ctrl.print && commRank == 0 )
            std::cout << "  alphaAffPri = " << alphaAffPri 
                      << ", alphaAffDual = " << alphaAffDual << std::endl;

        // Compute what the new duality measure would become
        // =================================================
        const Real mu = Dot(x,z) / n;
        // NOTE: dz and dx are used as temporaries
        dx = x;
        dz = z;
        Axpy( alphaAffPri,  dxAff, dx );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(dx,dz) / n;

        // Compute a centrality parameter using Mehotra's formula
        // ======================================================
        // TODO: Allow the user to override this function
        const Real sigma = Pow(muAff/mu,Real(3)); 
        if( ctrl.print && commRank == 0 )
            std::cout << "  muAff = " << muAff 
                      << ", mu = " << mu 
                      << ", sigma = " << sigma << std::endl;

        // Solve for the centering-corrector 
        // =================================
        Zeros( rc, n, 1 );
        Zeros( rb, m, 1 );
        // TODO: Find a more convenient means of expressing this operation
        if( dxAff.IsLocalCol(0) )
            for( Int iLoc=0; iLoc<dxAff.LocalHeight(); ++iLoc )
                rmu.SetLocal
                ( iLoc, 0, 
                  dxAff.GetLocal(iLoc,0)*dzAff.GetLocal(iLoc,0) - sigma*mu );
        if( ctrl.system == FULL_KKT )
        {
            // Construct the new full KKT RHS
            // ------------------------------
            KKTRHS( rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            lu::SolveAfter( NORMAL, J, p, d );
            ExpandKKTSolution( m, n, d, dx, dy, dz );
        }
        else // ctrl.system == AUGMENTED_KKT
        {
            // Construct the new "augmented" KKT RHS
            // -------------------------------------
            AugmentedKKTRHS( x, rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            ldl::SolveAfter( J, dSub, p, d, false );
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }

        // TODO: Residual checks for center-corrector

        // Add in the affine search direction
        // ==================================
        Axpy( Real(1), dxAff, dx );
        Axpy( Real(1), dyAff, dy );
        Axpy( Real(1), dzAff, dz );

        // Compute the max positive [0,1/maxStepRatio] step length
        // =======================================================
        Real alphaPri = 1/ctrl.maxStepRatio;
        // TODO: Find a more convenient means of expressing this operation
        if( x.IsLocalCol(0) )
        {
            for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
            {
                const Real xi = x.GetLocal(iLoc,0);
                const Real dxi = dx.GetLocal(iLoc,0);
                if( dxi < Real(0) )
                    alphaPri = Min(alphaPri,-xi/dxi);
            }
        }
        alphaPri = mpi::AllReduce( alphaPri, mpi::MIN, x.DistComm() );
        Real alphaDual = 1/ctrl.maxStepRatio;
        // TODO: Find a more convenient means of expressing this operation
        if( z.IsLocalCol(0) )
        {
            for( Int iLoc=0; iLoc<z.LocalHeight(); ++iLoc )
            {
                const Real zi = z.GetLocal(iLoc,0);
                const Real dzi = dz.GetLocal(iLoc,0);
                if( dzi < Real(0) )
                    alphaDual = Min(alphaDual,-zi/dzi);
            }
        }
        alphaDual = mpi::AllReduce( alphaDual, mpi::MIN, z.DistComm() );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.print && commRank == 0 )
            std::cout << "  alphaPri = " << alphaAffPri 
                      << ", alphaDual = " << alphaAffDual << std::endl;

        // Update the current estimates
        // ============================
        Axpy( alphaPri,  dx, x );
        Axpy( alphaDual, dy, y ); 
        Axpy( alphaDual, dz, z );
    }
}

template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::primal::Mehrotra"))    
    LogicError("Sequential sparse-direct solvers not yet supported");
}

template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  DistMultiVec<Real>& x, DistMultiVec<Real>& y, DistMultiVec<Real>& z,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::primal::Mehrotra"))    

    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);
    const Real epsilon = lapack::MachineEpsilon<Real>();

    const Real bNrm2 = Nrm2( b );

    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, invMap;
    DistSparseMatrix<Real> J(comm);
    DistSymmFrontTree<Real> JFrontTree;

    DistMultiVec<Real> d(comm), 
                       rmu(comm),   rb(comm),    rc(comm), 
                       dxAff(comm), dyAff(comm), dzAff(comm),
                       dx(comm),    dy(comm),    dz(comm);
    DistNodalMultiVec<Real> dNodal;

    DistMultiVec<Real> regCand(comm), reg(comm);
    DistNodalMultiVec<Real> regCandNodal, regNodal;

#ifndef EL_RELEASE
    DistMultiVec<Real> dxError(comm), dyError(comm), dzError(comm);
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
        // |c^T x + x^T Q x - b^T y| / (1 + |c^T x + 1/2 x^T Q x|) <= tol ?
        // ----------------------------------------------------------------
        Zeros( d, n, 1 );
        // NOTE: This assumes that Q is explicitly Hermitian
        Multiply( NORMAL, Real(1), Q, x, Real(0), d );
        const Real xTQx = Dot(x,d);
        const Real primObj = Dot(c,x) + xTQx/2;
        const Real dualObj = Dot(b,y) - xTQx/2; 
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        rb = b;
        Scale( Real(-1), rb );
        Multiply( NORMAL, Real(1), A, x, Real(1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c + Q x ||_2) <= tol ?
        // --------------------------------------------
        // NOTE: d currently contains Q x
        Axpy( Real(1), c, d );
        const Real objGradNrm2 = Nrm2( d );
        rc = d;
        Scale( Real(-1), rc );
        Multiply( TRANSPOSE, Real(1), A, y, Real(1), rc );
        Axpy( Real(1), z, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+objGradNrm2);
        // Now check the pieces
        // --------------------
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;
        else if( ctrl.print && commRank == 0 )
            std::cout << " iter " << numIts << ":\n"
                      << "  |primObj - dualObj| / (1 + |primObj|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)         = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c + Q x ||_2)   = "
                      << rcConv << std::endl;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // r_mu := X Z e
        // =============
        rmu.Resize( n, 1 );
        for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
            rmu.SetLocal( iLoc, 0, x.GetLocal(iLoc,0)*z.GetLocal(iLoc,0) );

        // Compute the affine search direction
        // ===================================
        //const Real minReductionFactor = 2;
        //const Int maxRefineIts = 10;
        // ctrl.system == AUGMENTED_KKT
        {
            // Construct the "normal" KKT system
            // ---------------------------------
            // TODO: Add default regularization
            AugmentedKKT( Q, A, x, z, J, false );
            AugmentedKKTRHS( x, rmu, rc, rb, d );
            const Real pivTol = MaxNorm(J)*epsilon;
            const Real regMagPrimal = Pow(epsilon,Real(0.75));
            const Real regMagLagrange = Pow(epsilon,Real(0.5));
            regCand.Resize( n+m, 1 );
            for( Int iLoc=0; iLoc<regCand.LocalHeight(); ++iLoc )
            {
                const Int i = regCand.FirstLocalRow() + iLoc;
                if( i < n )
                    regCand.SetLocal( iLoc, 0, -regMagPrimal );
                else
                    regCand.SetLocal( iLoc, 0, regMagLagrange );
            }
            // Do not use any a priori regularization
            Zeros( reg, m+n, 1 );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            if( numIts == 0 )
            {
                NestedDissection( J.LockedDistGraph(), map, sepTree, info );
                map.FormInverse( invMap );
            }
            JFrontTree.Initialize( J, map, sepTree, info );
            regCandNodal.Pull( invMap, info, regCand );
            regNodal.Pull( invMap, info, reg );
            RegularizedLDL
            ( info, JFrontTree, pivTol, regCandNodal, regNodal, LDL_1D );
            regNodal.Push( invMap, info, reg );
            // NOTE: Need to modify iterative refinement procedure
            /*
            SolveWithIterativeRefinement
            ( J, invMap, info, JFrontTree, d, 
              minReductionFactor, maxRefineIts );
            */
            dNodal.Pull( invMap, info, d );
            Solve( info, JFrontTree, dNodal );
            dNodal.Push( invMap, info, d );
            ExpandAugmentedSolution( x, z, rmu, d, dxAff, dyAff, dzAff );
        }

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        Real rmuNrm2 = Nrm2( rmu ); 
        dzError = rmu;
        for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
        {
            const Real xi = x.GetLocal(iLoc,0);
            const Real zi = z.GetLocal(iLoc,0);
            const Real dxi = dxAff.GetLocal(iLoc,0);
            const Real dzi = dzAff.GetLocal(iLoc,0);
            dzError.UpdateLocal( iLoc, 0, xi*dzi + zi*dxi );
        }
        Real dzErrorNrm2 = Nrm2( dzError );

        dyError = dzAff;
        Multiply( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Multiply( NORMAL, Real(-1), Q, dxAff, Real(1), dyError );
        Axpy( Real(1), rc, dyError );
        Real dyErrorNrm2 = Nrm2( dyError );

        dxError = rb;
        Multiply( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        Real dxErrorNrm2 = Nrm2( dxError );

        if( ctrl.print && commRank == 0 )
            std::cout << "  || dxAffError ||_2 / (1 + || r_b ||_2)  = " 
                      << dxErrorNrm2/(1+rbNrm2) << "\n"
                      << "  || dyAffError ||_2 / (1 + || r_c ||_2)  = " 
                      << dyErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dzAffError ||_2 / (1 + || r_mu ||_2) = " 
                      << dzErrorNrm2/(1+rmuNrm2) << std::endl;
#endif

        // Compute the maximum affine [0,1]-step which preserves positivity
        // ================================================================
        Real alphaAffPri = 1; 
        for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
        {
            const Real xi = x.GetLocal(iLoc,0);
            const Real dxAffi = dxAff.GetLocal(iLoc,0);
            if( dxAffi < Real(0) )
                alphaAffPri = Min(alphaAffPri,-xi/dxAffi);
        }
        alphaAffPri = mpi::AllReduce( alphaAffPri, mpi::MIN, comm );
        Real alphaAffDual = 1; 
        for( Int iLoc=0; iLoc<z.LocalHeight(); ++iLoc )
        {
            const Real zi = z.GetLocal(iLoc,0);
            const Real dzAffi = dzAff.GetLocal(iLoc,0);
            if( dzAffi < Real(0) )
                alphaAffDual = Min(alphaAffDual,-zi/dzAffi);
        }
        alphaAffDual = mpi::AllReduce( alphaAffDual, mpi::MIN, comm );
        if( ctrl.print && commRank == 0 )
            std::cout << "  alphaAffPri = " << alphaAffPri 
                      << ", alphaAffDual = " << alphaAffDual << std::endl;

        // Compute what the new duality measure would become
        // =================================================
        const Real mu = Dot(x,z) / n;
        // NOTE: dz and dx are used as temporaries
        dx = x;
        dz = z;
        Axpy( alphaAffPri,  dxAff, dx );
        Axpy( alphaAffDual, dzAff, dz );
        const Real muAff = Dot(dx,dz) / n;

        // Compute a centrality parameter using Mehotra's formula
        // ======================================================
        // TODO: Allow the user to override this function
        const Real sigma = Pow(muAff/mu,Real(3)); 
        if( ctrl.print && commRank == 0 )
            std::cout << "  muAff = " << muAff 
                      << ", mu = " << mu 
                      << ", sigma = " << sigma << std::endl;

        // Solve for the centering-corrector 
        // =================================
        Zeros( rc, n, 1 );
        Zeros( rb, m, 1 );
        for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
            rmu.SetLocal
            ( iLoc, 0, 
              dxAff.GetLocal(iLoc,0)*dzAff.GetLocal(iLoc,0) - sigma*mu );
        // ctrl.system == AUGMENTED_KKT
        {
            // Construct the new "normal" KKT RHS
            // ----------------------------------
            AugmentedKKTRHS( x, rmu, rc, rb, d );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            // TODO: Iterative refinement
            dNodal.Pull( invMap, info, d );
            Solve( info, JFrontTree, dNodal );
            dNodal.Push( invMap, info, d );
            ExpandAugmentedSolution( x, z, rmu, d, dx, dy, dz );
        }

        // TODO: Residual checks for center-corrector

        // Add in the affine search direction
        // ==================================
        Axpy( Real(1), dxAff, dx );
        Axpy( Real(1), dyAff, dy );
        Axpy( Real(1), dzAff, dz );

        // Compute the max positive [0,1/maxStepRatio] step length
        // =======================================================
        Real alphaPri = 1/ctrl.maxStepRatio;
        for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
        {
            const Real xi = x.GetLocal(iLoc,0);
            const Real dxi = dx.GetLocal(iLoc,0);
            if( dxi < Real(0) )
                alphaPri = Min(alphaPri,-xi/dxi);
        }
        alphaPri = mpi::AllReduce( alphaPri, mpi::MIN, comm );
        Real alphaDual = 1/ctrl.maxStepRatio;
        for( Int iLoc=0; iLoc<z.LocalHeight(); ++iLoc )
        {
            const Real zi = z.GetLocal(iLoc,0);
            const Real dzi = dz.GetLocal(iLoc,0);
            if( dzi < Real(0) )
                alphaDual = Min(alphaDual,-zi/dzi);
        }
        alphaDual = mpi::AllReduce( alphaDual, mpi::MIN, comm );
        alphaPri = Min(ctrl.maxStepRatio*alphaPri,Real(1));
        alphaDual = Min(ctrl.maxStepRatio*alphaDual,Real(1));
        if( ctrl.print && commRank == 0 )
            std::cout << "  alphaPri = " << alphaPri 
                      << ", alphaDual = " << alphaDual << std::endl;

        // Update the current estimates
        // ============================
        Axpy( alphaPri,  dx, x );
        Axpy( alphaDual, dy, y ); 
        Axpy( alphaDual, dz, z );
    }
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    AbstractDistMatrix<Real>& x, AbstractDistMatrix<Real>& y, \
    AbstractDistMatrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    Matrix<Real>& x, Matrix<Real>& y, Matrix<Real>& z, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    DistMultiVec<Real>& x, DistMultiVec<Real>& y, DistMultiVec<Real>& z, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace primal
} // namespace qp
} // namespace El
