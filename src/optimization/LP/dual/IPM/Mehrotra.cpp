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
namespace lp {
namespace dual {

// The following solves a linear program in "dual" conic form:
//
//   min c^T x
//   s.t. A x = b, G x + s = h, s >= 0,
//
// as opposed to the more specific "primal" conic form:
//
//   min c^T x
//   s.t. A x = b, x >= 0,
//
// using a Mehrotra Predictor-Corrector scheme.
//

template<typename Real>
void Mehrotra
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,       Matrix<Real>& y, 
        Matrix<Real>& z,       Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::Mehrotra"))    
    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    // TODO: Use the initialization suggested by Vandenberghe
    if( !ctrl.initialized )
    {
        Gaussian( x, n, 1 );
        Zeros( y, m, 1 );
        Uniform( z, k, 1, Real(0.5), Real(0.49) );
        s = h;
        Gemv( NORMAL, Real(-1), G, x, Real(1), s );
        LowerClip( s, Real(1) );
    }

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
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;
        else if( ctrl.print )
            std::cout << " iter " << numIts << ":\n"
                      << "  |primal - dual| / (1 + |primal|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)   = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c ||_2)   = "
                      << rcConv << "\n"
                      << "  || r_h ||_2 / (1 + || h ||_2)   = "
                      << rhConv << std::endl;

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
        KKT( A, G, s, z, J );
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        LDL( J, dSub, p, false );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        // TODO: dmuError

        if( ctrl.print )
            std::cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                      << dxErrorNrm2/(1+rbNrm2) << "\n"
                      << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                      << dyErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dzError ||_2 / (1 + || r_mu ||_2) = "
                      << dzErrorNrm2/(1+rhNrm2) << std::endl;
#endif

        // Compute the max affine [0,1]-step which keeps s and z in the cone
        // =================================================================
        const Real alphaAffPri = MaxStepInPositiveCone( s, dsAff, Real(1) );
        const Real alphaAffDual = MaxStepInPositiveCone( z, dzAff, Real(1) );
        if( ctrl.print )
            std::cout << "  alphaAffPri = " << alphaAffPri
                      << ", alphaAffDual = " << alphaAffDual << std::endl;

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
            std::cout << "  muAff = " << muAff
                      << ", mu = " << mu
                      << ", sigma = " << sigma << std::endl;

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
            std::cout << "  alphaPri = " << alphaPri
                      << ", alphaDual = " << alphaDual << std::endl;

        // Update the current estimates
        // ============================
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
    }
}

template<typename Real>
void Mehrotra
( const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& GPre,
  const AbstractDistMatrix<Real>& b,    const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
        AbstractDistMatrix<Real>& xPre,       AbstractDistMatrix<Real>& yPre, 
        AbstractDistMatrix<Real>& zPre,       AbstractDistMatrix<Real>& sPre,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::Mehrotra"))    

    ProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = 0;
    proxCtrl.rowAlign = 0;
    auto APtr = ReadProxy<Real,MC,MR>(&APre,proxCtrl);      auto& A = *APtr;
    auto GPtr = ReadProxy<Real,MC,MR>(&GPre,proxCtrl);      auto& G = *GPtr;
    // NOTE: The following do not need to be read proxies when !ctrl.initialized
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,proxCtrl); auto& x = *xPtr;
    auto yPtr = ReadWriteProxy<Real,MC,MR>(&yPre,proxCtrl); auto& y = *yPtr;
    auto zPtr = ReadWriteProxy<Real,MC,MR>(&zPre,proxCtrl); auto& z = *zPtr;
    auto sPtr = ReadWriteProxy<Real,MC,MR>(&sPre,proxCtrl); auto& s = *sPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Int k = G.Height();
    const Grid& grid = A.Grid();
    const Int commRank = A.Grid().Rank();

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

    // TODO: Use the initialization suggested by Vandenberghe
    if( !ctrl.initialized )
    {
        Gaussian( x, n, 1 );
        Zeros( y, m, 1 );
        Uniform( z, k, 1, Real(0.5), Real(0.49) );
        Copy( h, s );
        Gemv( NORMAL, Real(-1), G, x, Real(1), s );
        LowerClip( s, Real(1) );
    }

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
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;
        else if( ctrl.print && commRank == 0 )
            std::cout << " iter " << numIts << ":\n"
                      << "  |primal - dual| / (1 + |primal|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)   = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c ||_2)   = "
                      << rcConv << "\n"
                      << "  || r_h ||_2 / (1 + || h ||_2)   = "
                      << rhConv << std::endl;

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
        KKT( A, G, s, z, J );
        KKTRHS( rc, rb, rh, rmu, z, d );
        // Compute the proposed step from the KKT system
        // ---------------------------------------------
        LDL( J, dSub, p, false );
        ldl::SolveAfter( J, dSub, p, d, false );
        ExpandSolution( m, n, d, rmu, s, z, dxAff, dyAff, dzAff, dsAff );
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dyAff, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dzAff, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dxAff, Real(1), dzError );
        Axpy( Real(1), dsAff, dzError );
        const Real dzErrorNrm2 = Nrm2( dzError );

        // TODO: dmuError

        if( ctrl.print && commRank == 0 )
            std::cout << "  || dxError ||_2 / (1 + || r_b ||_2) = "
                      << dxErrorNrm2/(1+rbNrm2) << "\n"
                      << "  || dyError ||_2 / (1 + || r_c ||_2) = "
                      << dyErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dzError ||_2 / (1 + || r_mu ||_2) = "
                      << dzErrorNrm2/(1+rhNrm2) << std::endl;
#endif
 
        // Compute the max affine [0,1]-step which keeps s and z in the cone
        // =================================================================
        const Real alphaAffPri = MaxStepInPositiveCone( s, dsAff, Real(1) );
        const Real alphaAffDual = MaxStepInPositiveCone( z, dzAff, Real(1) );
        if( ctrl.print && commRank == 0 )
            std::cout << "  alphaAffPri = " << alphaAffPri
                      << ", alphaAffDual = " << alphaAffDual << std::endl;

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
            std::cout << "  muAff = " << muAff
                      << ", mu = " << mu
                      << ", sigma = " << sigma << std::endl;

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
            std::cout << "  alphaPri = " << alphaPri
                      << ", alphaDual = " << alphaDual << std::endl;

        // Update the current estimates
        // ============================
        Axpy( alphaPri,  dx, x );
        Axpy( alphaPri,  ds, s );
        Axpy( alphaDual, dy, y );
        Axpy( alphaDual, dz, z );
    }
}

template<typename Real>
void Mehrotra
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
        Matrix<Real>& x,             Matrix<Real>& y, 
        Matrix<Real>& z,             Matrix<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::Mehrotra"))    
    LogicError("Sequential sparse-direct solvers not yet supported");
}

template<typename Real>
void Mehrotra
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
        DistMultiVec<Real>& x,           DistMultiVec<Real>& y, 
        DistMultiVec<Real>& z,           DistMultiVec<Real>& s,
  const MehrotraCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::Mehrotra"))    
    LogicError("This routine is not yet written");
}

#define PROTO(Real) \
  template void Mehrotra \
  ( const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,       Matrix<Real>& y, \
          Matrix<Real>& z,       Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
          AbstractDistMatrix<Real>& x,       AbstractDistMatrix<Real>& y, \
          AbstractDistMatrix<Real>& z,       AbstractDistMatrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
    const Matrix<Real>& h, \
          Matrix<Real>& x,             Matrix<Real>& y, \
          Matrix<Real>& z,             Matrix<Real>& s, \
    const MehrotraCtrl<Real>& ctrl ); \
  template void Mehrotra \
  ( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
          DistMultiVec<Real>& x,           DistMultiVec<Real>& y, \
          DistMultiVec<Real>& z,           DistMultiVec<Real>& s, \
    const MehrotraCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace dual
} // namespace lp
} // namespace El
