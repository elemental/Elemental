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
// using a simple Infeasible Path Following (IPF) scheme. This routine
// should only be used for academic purposes, as the Mehrotra alternative
// typically requires an order of magnitude fewer iterations.

template<typename Real>
void IPF
( const Matrix<Real>& A, const Matrix<Real>& G,
  const Matrix<Real>& b, const Matrix<Real>& c,
  const Matrix<Real>& h,
  Matrix<Real>& x, Matrix<Real>& y, 
  Matrix<Real>& z, Matrix<Real>& s,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::IPF"))    
    const Int m = A.Height();
    const Int n = A.Width();

    const Real bNrm2 = Nrm2( b );
    const Real cNrm2 = Nrm2( c );
    const Real hNrm2 = Nrm2( h );

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
        // |c^T x - (-b^T y + h^T z)| / (1 + |c^T x|) <= tol ?
        // ---------------------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = -Dot(b,y) + Dot(h,z);
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
        Gemv( NORMAL, Real(1), G, x, Real(1), rb );
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

        // Form the residual for the scaled equation, S Z - sigma mu e
        // ===========================================================
        const Real mu = Dot(s,z) / n;
        rmu.Resize( n, 1 );
        for( Int i=0; i<n; ++i )
            rmu.Set( i, 0, s.Get(i,0)*z.Get(i,0) - ctrl.centering*mu );

        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // =============================
            KKT( A, G, s, z, J );
            KKTRHS( rc, rb, rh, rmu, z, d );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, d );
            ExpandSolution( m, n, d, rmu, s, z, dx, dy, dz, ds );
        }
        else
            LogicError("Invalid KKT system choice");
#ifndef EL_RELEASE
        // Sanity checks
        // =============
        dyError = rc;
        Gemv( TRANSPOSE, Real(1), A, dy, Real(1), dyError );
        Gemv( TRANSPOSE, Real(1), G, dz, Real(1), dyError );
        const Real dyErrorNrm2 = Nrm2( dyError );

        dxError = rb;
        Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        dzError = rh;
        Gemv( NORMAL, Real(1), G, dx, Real(1), dzError );
        Axpy( Real(1), ds, dzError );
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

        // Take a step in the computed direction
        // =====================================
        const Real alpha =
          IPFLineSearch
          ( A, G, b, c, h, x, y, z, s, dx, dy, dz, ds,
            ctrl.tol*(1+bNrm2), ctrl.tol*(1+cNrm2), ctrl.tol*(1+hNrm2),
            ctrl.lineSearchCtrl );
        if( ctrl.print )
            std::cout << "  alpha = " << alpha << std::endl;
        Axpy( alpha, dx, x );
        Axpy( alpha, dy, y );
        Axpy( alpha, dz, z );
        Axpy( alpha, ds, s );
    }
}

template<typename Real>
void IPF
( const AbstractDistMatrix<Real>& APre, const AbstractDistMatrix<Real>& GPre,
  const AbstractDistMatrix<Real>& b,    const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& h,
  AbstractDistMatrix<Real>& xPre, AbstractDistMatrix<Real>& y, 
  AbstractDistMatrix<Real>& zPre, AbstractDistMatrix<Real>& s,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::IPF"))    

    /*
    ProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = 0;
    proxCtrl.rowAlign = 0;
    auto APtr = ReadProxy<Real,MC,MR>(&APre,proxCtrl);      auto& A = *APtr;
    auto GPtr = ReadProxy<Real,MC,MR>(&GPre,proxCtrl);      auto& G = *GPtr;
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,proxCtrl); auto& x = *xPtr;
    auto zPtr = ReadWriteProxy<Real,MC,MR>(&zPre,proxCtrl); auto& z = *zPtr;
    */

    LogicError("This routine is not yet written");
}

template<typename Real>
void IPF
( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G,
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& h,
  Matrix<Real>& x, Matrix<Real>& y, 
  Matrix<Real>& z, Matrix<Real>& s,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::IPF"))    
    LogicError("Sequential sparse-direct solvers not yet supported");
}

template<typename Real>
void IPF
( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G,
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& h,
  DistMultiVec<Real>& x, DistMultiVec<Real>& y, 
  DistMultiVec<Real>& z, DistMultiVec<Real>& s,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lp::dual::IPF"))    
    LogicError("This routine is not yet written");
}

#define PROTO(Real) \
  template void IPF \
  ( const Matrix<Real>& A, const Matrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& h, \
    Matrix<Real>& x, Matrix<Real>& y, \
    Matrix<Real>& z, Matrix<Real>& s, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const AbstractDistMatrix<Real>& A, const AbstractDistMatrix<Real>& G, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& h, \
    AbstractDistMatrix<Real>& x, AbstractDistMatrix<Real>& y, \
    AbstractDistMatrix<Real>& z, AbstractDistMatrix<Real>& s, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const SparseMatrix<Real>& A, const SparseMatrix<Real>& G, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    const Matrix<Real>& h, \
    Matrix<Real>& x, Matrix<Real>& y, \
    Matrix<Real>& z, Matrix<Real>& s, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const DistSparseMatrix<Real>& A, const DistSparseMatrix<Real>& G, \
    const DistMultiVec<Real>& b, const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& h, \
    DistMultiVec<Real>& x, DistMultiVec<Real>& y, \
    DistMultiVec<Real>& z, DistMultiVec<Real>& s, \
    const IPFCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace dual
} // namespace lp
} // namespace El
