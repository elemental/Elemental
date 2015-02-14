/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace qp {
namespace direct {

// TODO: Expose maximum number of trials as parameter of IPFLineSearchCtrl

template<typename Real>
Real IPFLineSearch
( const Matrix<Real>& Q,  const Matrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& x,  const Matrix<Real>& y,  
  const Matrix<Real>& z,
  const Matrix<Real>& dx, const Matrix<Real>& dy, 
  const Matrix<Real>& dz,
  Real upperBound,
  Real bTol, Real cTol,
  const qp::IPFLineSearchCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::IPFLineSearch"))
    if( ctrl.gamma <= Real(0) || ctrl.gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( ctrl.beta < Real(1) )
        LogicError("beta must be at least one");
    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || y.Width() != 1 || z.Width() != 1 ||
        dx.Width() != 1 || dy.Width() != 1 || dz.Width() != 1 )
        LogicError("{b,c,x,y,z,dx,dy,dz} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-infinity" neighborhood of the central path, i.e.,
    //  (a) || r_b(alpha) ||_2 <= || r_b ||_2 beta mu(alpha) / mu,
    //  (b) || r_c(alpha) ||_2 <= || r_c ||_2 beta mu(alpha) / mu,
    //  (c) x(alpha), z(alpha) > 0, and, for all i,
    //  (d) x_i(alpha) z_i(alpha) >= gamma mu(alpha),
    // where 
    //    x(alpha) = x + alpha dx,
    //    y(alpha) = y + alpha dy,
    //    z(alpha) = z + alpha dz,
    //    r_b(alpha) = A x(alpha) - b, and
    //    r_c(alpha) = Q x(alpha) + A^T y(alpha) - z(alpha) + c,
    // and the Armijo condition,
    //   mu(alpha) <= (1 - alpha/\psi) mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    Matrix<Real> Q_x, Q_dx, A_x, A_dx, AT_y, AT_dy, rb, rc;
    Zeros( Q_x,   n, 1 );
    Zeros( Q_dx,  n, 1 );
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_y,  n, 1 );
    Zeros( AT_dy, n, 1 );
    Hemv( LOWER,     Real(1), Q, x,  Real(0), Q_x   );
    Hemv( LOWER,     Real(1), Q, dx, Real(0), Q_dx  );
    Gemv( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Gemv( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Gemv( TRANSPOSE, Real(1), A, y,  Real(0), AT_y  );
    Gemv( TRANSPOSE, Real(1), A, dy, Real(0), AT_dy );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = Q_x;
    Axpy( Real( 1), AT_y, rc );
    Axpy( Real(-1), z,    rc );
    Axpy( Real( 1), c,    rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,z) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = upperBound;
    Matrix<Real> x_alpha, z_alpha, rb_alpha, rc_alpha;
    for( Int k=0; k<100; ++k, alpha=alpha/ctrl.stepRatio )
    {
        // x(alpha) = x + alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // z(alpha) = z + alpha dz
        // ^^^^^^^^^^^^^^^^^^^^^^^ 
        z_alpha = z;
        Axpy( alpha, dz, z_alpha );

        // mu(alpha) = x(alpha)^T z / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,z_alpha) / n;
        if( mu_alpha > (1-alpha/ctrl.psi)*mu )
        {
            if( ctrl.print )
                cout << "Armijo condition not satisfied" << endl;
            continue;
        }

        // Check 
        //    x(alpha), z(alpha) > 0, and 
        //    x_i(alpha) z_i(alpha) >= gamma mu(alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        bool balanced = true;
        for( Int i=0; i<n; ++i )
        {
            const Real xi_alpha = x_alpha.Get(i,0);
            const Real zi_alpha = z_alpha.Get(i,0);
            if( xi_alpha <= Real(0) || zi_alpha <= Real(0) )
                balanced = false;
            if( xi_alpha*zi_alpha < ctrl.gamma*mu_alpha )
                balanced = false;
        }
        if( !balanced )
        {
            if( ctrl.print )
                cout << "  unbalanced entries" << endl;
            continue;
        }
        // Check ||r_b(alpha)||_2 <= Max(bTol,||r_b||_2 beta mu(alpha)/mu)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > Max(bTol,rbNrm2*ctrl.beta*mu_alpha/mu) )
        {
            if( ctrl.print )
                cout << "  r_b failure: " << rb_alphaNrm2 << " > Max("
                     << bTol << "," << rbNrm2*ctrl.beta*mu_alpha/mu << ")"
                     << endl;
            continue;
        }
        // Check ||r_c(alpha)||_2 <= Max(cTol,||r_c||_2 beta mu(alpha)/mu)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy(  alpha, Q_dx,  rc_alpha );
        Axpy(  alpha, AT_dy, rc_alpha );
        Axpy( -alpha, dz,    rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > Max(cTol,rcNrm2*ctrl.beta*mu_alpha/mu) )
        {
            if( ctrl.print )
                cout << "  r_c failure: " << rc_alphaNrm2 << " > Max("
                     << cTol << "," << rcNrm2*ctrl.beta*mu_alpha/mu << ")"
                     << endl;
        }
        else
            break;
    }
    return alpha;
}

template<typename Real>
Real IPFLineSearch
( const AbstractDistMatrix<Real>& Q,    const AbstractDistMatrix<Real>& APre, 
  const AbstractDistMatrix<Real>& b,    const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& x,    const AbstractDistMatrix<Real>& y,
  const AbstractDistMatrix<Real>& z,
  const AbstractDistMatrix<Real>& dx,   const AbstractDistMatrix<Real>& dy, 
  const AbstractDistMatrix<Real>& dz,
  Real upperBound,
  Real bTol, Real cTol,
  const qp::IPFLineSearchCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::IPFLineSearch"))
    if( ctrl.gamma <= Real(0) || ctrl.gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( ctrl.beta < Real(1) )
        LogicError("beta must be at least one");

    auto APtr = ReadProxy<Real,MC,MR>(&APre);
    auto& A = *APtr;

    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || y.Width() != 1 || z.Width() != 1 ||
        dx.Width() != 1 || dy.Width() != 1 || dz.Width() != 1 )
        LogicError("{b,c,x,y,z,dx,dy,dz} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const Int commRank = grid.Rank();

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-infinity" neighborhood of the central path, i.e.,
    //  (a) || r_b(alpha) ||_2 <= || r_b ||_2 beta mu(alpha) / mu,
    //  (b) || r_c(alpha) ||_2 <= || r_c ||_2 beta mu(alpha) / mu,
    //  (c) x(alpha), z(alpha) > 0, and, for all i,
    //  (d) x_i(alpha) z_i(alpha) >= gamma mu(alpha),
    // where 
    //    x(alpha) = x + alpha dx,
    //    y(alpha) = y + alpha dy,
    //    z(alpha) = z + alpha dz,
    //    r_b(alpha) = A x(alpha) - b, and
    //    r_c(alpha) = Q x(alpha) + A^T y(alpha) - z(alpha) + c,
    // and the Armijo condition,
    //   mu(alpha) <= (1 - alpha/\psi) mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    DistMatrix<Real> Q_x(grid), Q_dx(grid), A_x(grid), A_dx(grid), 
                     AT_y(grid), AT_dy(grid), 
                     rb(grid), rc(grid);
    Zeros( Q_x,   n, 1 );
    Zeros( Q_dx,  n, 1 );
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_y,  n, 1 );
    Zeros( AT_dy, n, 1 );
    Hemv( LOWER,     Real(1), Q, x,  Real(0), Q_x   );
    Hemv( LOWER,     Real(1), Q, dx, Real(0), Q_dx  );
    Gemv( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Gemv( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Gemv( TRANSPOSE, Real(1), A, y,  Real(0), AT_y  );
    Gemv( TRANSPOSE, Real(1), A, dy, Real(0), AT_dy );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = Q_x;
    Axpy( Real( 1), AT_y, rc );
    Axpy( Real(-1), z,    rc );
    Axpy( Real( 1), c,    rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,z) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = upperBound;
    DistMatrix<Real> x_alpha(grid), z_alpha(grid), 
                     rb_alpha(grid), rc_alpha(grid);
    x_alpha.Align(0,0);
    z_alpha.Align(0,0);
    for( Int k=0; k<100; ++k, alpha=alpha/ctrl.stepRatio )
    {
        // x(alpha) = x + alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // z(alpha) = z + alpha dz
        // ^^^^^^^^^^^^^^^^^^^^^^^ 
        z_alpha = z;
        Axpy( alpha, dz, z_alpha );

        // mu(alpha) = x(alpha)^T z / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,z_alpha) / n;
        if( mu_alpha > (1-alpha/ctrl.psi)*mu )
        {
            if( ctrl.print && commRank == 0 )
                cout << "Armijo condition not satisfied" << endl;
            continue;
        }

        // Check 
        //    x(alpha), z(alpha) > 0, and 
        //    x_i(alpha) z_i(alpha) >= gamma mu(alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        byte locallyBalanced = true;
        if( x_alpha.IsLocalCol(0) )
        {
            for( Int iLoc=0; iLoc<x_alpha.LocalHeight(); ++iLoc )
            {
                const Real xi_alpha = x_alpha.GetLocal(iLoc,0);
                const Real zi_alpha = z_alpha.GetLocal(iLoc,0);
                if( xi_alpha <= Real(0) || zi_alpha <= Real(0) )
                    locallyBalanced = false;
                if( xi_alpha*zi_alpha < ctrl.gamma*mu_alpha )
                    locallyBalanced = false;
            }
        }
        const byte balanced =
            mpi::AllReduce( locallyBalanced, mpi::BINARY_AND, grid.VCComm() );
        if( !balanced )
        {
            if( ctrl.print && commRank == 0 )
                cout << "  unbalanced entries" << endl;
            continue;
        }
        // Check ||r_b(alpha)||_2 <= Max(bTol,||r_b||_2 beta mu(alpha)/mu)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > Max(bTol,rbNrm2*ctrl.beta*mu_alpha/mu) )
        {
            if( ctrl.print && commRank == 0 )
                cout << "  r_b failure: " << rb_alphaNrm2 << " > Max("
                     << bTol << "," << rbNrm2*ctrl.beta*mu_alpha/mu << ")"
                     << endl;
            continue;
        }
        // Check ||r_c(alpha)||_2 <= Max(cTol,||r_c||_2 beta mu(alpha)/mu)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy(  alpha, Q_dx,  rc_alpha );
        Axpy(  alpha, AT_dy, rc_alpha );
        Axpy( -alpha, dz,    rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > Max(cTol,rcNrm2*ctrl.beta*mu_alpha/mu) )
        {
            if( ctrl.print && commRank == 0 )
                cout << "  r_c failure: " << rc_alphaNrm2 << " > Max("
                     << cTol << "," << rcNrm2*ctrl.beta*mu_alpha/mu << ")"
                     << endl;
        }
        else
            break;
    }
    return alpha;
}

template<typename Real>
Real IPFLineSearch
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,       const Matrix<Real>& c,
  const Matrix<Real>& x,       const Matrix<Real>& y,  
  const Matrix<Real>& z,
  const Matrix<Real>& dx,      const Matrix<Real>& dy, 
  const Matrix<Real>& dz,
  Real upperBound,
  Real bTol, Real cTol,
  const qp::IPFLineSearchCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::IPFLineSearch"))
    if( ctrl.gamma <= Real(0) || ctrl.gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( ctrl.beta < Real(1) )
        LogicError("beta must be at least one");
    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || y.Width() != 1 || z.Width() != 1 ||
        dx.Width() != 1 || dy.Width() != 1 || dz.Width() != 1 )
        LogicError("{b,c,x,y,z,dx,dy,dz} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-infinity" neighborhood of the central path, i.e.,
    //  (a) || r_b(alpha) ||_2 <= || r_b ||_2 beta mu(alpha) / mu,
    //  (b) || r_c(alpha) ||_2 <= || r_c ||_2 beta mu(alpha) / mu,
    //  (c) x(alpha), z(alpha) > 0, and, for all i,
    //  (d) x_i(alpha) z_i(alpha) >= gamma mu(alpha),
    // where 
    //    x(alpha) = x + alpha dx,
    //    y(alpha) = y + alpha dy,
    //    z(alpha) = z + alpha dz,
    //    r_b(alpha) = A x(alpha) - b, and
    //    r_c(alpha) = Q x(alpha) + A^T y(alpha) - z(alpha) + c,
    // and the Armijo condition,
    //   mu(alpha) <= (1 - alpha/\psi) mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    Matrix<Real> Q_x, Q_dx, A_x, A_dx, AT_y, AT_dy, rb, rc;
    Zeros( Q_x,   n, 1 );
    Zeros( Q_dx,  n, 1 );
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_y,  n, 1 );
    Zeros( AT_dy, n, 1 );
    // NOTE: This requires Q to be explicitly symmetric
    Multiply( NORMAL,    Real(1), Q, x,  Real(0), Q_x   );
    Multiply( NORMAL,    Real(1), Q, dx, Real(0), Q_dx  );
    Multiply( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Multiply( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Multiply( TRANSPOSE, Real(1), A, y,  Real(0), AT_y  );
    Multiply( TRANSPOSE, Real(1), A, dy, Real(0), AT_dy );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = Q_x;
    Axpy( Real( 1), AT_y, rc );
    Axpy( Real(-1), z,    rc );
    Axpy( Real( 1), c,    rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,z) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = upperBound;
    Matrix<Real> x_alpha, z_alpha, rb_alpha, rc_alpha;
    for( Int k=0; k<100; ++k, alpha=alpha/ctrl.stepRatio )
    {
        // x(alpha) = x + alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // z(alpha) = z + alpha dz
        // ^^^^^^^^^^^^^^^^^^^^^^^ 
        z_alpha = z;
        Axpy( alpha, dz, z_alpha );

        // mu(alpha) = x(alpha)^T z / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,z_alpha) / n;
        if( mu_alpha > (1-alpha/ctrl.psi)*mu )
        {
            if( ctrl.print )
                cout << "Armijo condition not satisfied" << endl;
            continue;
        }

        // Check 
        //    x(alpha), z(alpha) > 0, and 
        //    x_i(alpha) z_i(alpha) >= gamma mu(alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        bool balanced = true;
        for( Int i=0; i<n; ++i )
        {
            const Real xi_alpha = x_alpha.Get(i,0);
            const Real zi_alpha = z_alpha.Get(i,0);
            if( xi_alpha <= Real(0) || zi_alpha <= Real(0) )
                balanced = false;
            if( xi_alpha*zi_alpha < ctrl.gamma*mu_alpha )
                balanced = false;
        }
        if( !balanced )
        {
            if( ctrl.print )
                cout << "  unbalanced entries" << endl;
            continue;
        }
        // Check ||r_b(alpha)||_2 <= Max(bTol,||r_b||_2 beta mu(alpha)/mu)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > Max(bTol,rbNrm2*ctrl.beta*mu_alpha/mu) )
        {
            if( ctrl.print )
                cout << "  r_b failure: " << rb_alphaNrm2 << " > Max("
                     << bTol << "," << rbNrm2*ctrl.beta*mu_alpha/mu << ")"
                     << endl;
            continue;
        }
        // Check ||r_c(alpha)||_2 <= Max(cTol,||r_c||_2 beta mu(alpha)/mu)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy(  alpha, Q_dx,  rc_alpha );
        Axpy(  alpha, AT_dy, rc_alpha );
        Axpy( -alpha, dz,    rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > Max(cTol,rcNrm2*ctrl.beta*mu_alpha/mu) )
        {
            if( ctrl.print )
                cout << "  r_c failure: " << rc_alphaNrm2 << " > Max("
                     << cTol << "," << rcNrm2*ctrl.beta*mu_alpha/mu << ")"
                     << endl;
        }
        else
            break;
    }
    return alpha;
}

template<typename Real>
Real IPFLineSearch
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& x,     const DistMultiVec<Real>& y, 
  const DistMultiVec<Real>& z,
  const DistMultiVec<Real>& dx,    const DistMultiVec<Real>& dy, 
  const DistMultiVec<Real>& dz,
  Real upperBound,
  Real bTol, Real cTol,
  const qp::IPFLineSearchCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("qp::direct::IPFLineSearch"))
    if( ctrl.gamma <= Real(0) || ctrl.gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( ctrl.beta < Real(1) )
        LogicError("beta must be at least one");
    // TODO: Ensure communicators are congruent
    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || y.Width() != 1 || z.Width() != 1 ||
        dx.Width() != 1 || dy.Width() != 1 || dz.Width() != 1 )
        LogicError("{b,c,x,y,z,dx,dy,dz} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();
    const Int nLocal = x.LocalHeight();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-infinity" neighborhood of the central path, i.e.,
    //  (a) || r_b(alpha) ||_2 <= || r_b ||_2 beta mu(alpha) / mu,
    //  (b) || r_c(alpha) ||_2 <= || r_c ||_2 beta mu(alpha) / mu,
    //  (c) x(alpha), z(alpha) > 0, and, for all i,
    //  (d) x_i(alpha) z_i(alpha) >= gamma mu(alpha),
    // where 
    //    x(alpha) = x + alpha dx,
    //    y(alpha) = y + alpha dy,
    //    z(alpha) = z + alpha dz,
    //    r_b(alpha) = A x(alpha) - b, and
    //    r_c(alpha) = Q x(alpha) + A^T y(alpha) - z(alpha) + c,
    // and the Armijo condition,
    //   mu(alpha) <= (1 - alpha/\psi) mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    DistMultiVec<Real> Q_x(comm), Q_dx(comm), A_x(comm), A_dx(comm), 
                       AT_y(comm), AT_dy(comm), 
                       rb(comm), rc(comm);
    Zeros( Q_x,   n, 1 );
    Zeros( Q_dx,  n, 1 );
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_y,  n, 1 );
    Zeros( AT_dy, n, 1 );
    // NOTE: The following require Q to be explicitly symmetric
    Multiply( NORMAL,    Real(1), Q, x,  Real(0), Q_x   );
    Multiply( NORMAL,    Real(1), Q, dx, Real(0), Q_dx  );
    Multiply( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Multiply( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Multiply( TRANSPOSE, Real(1), A, y,  Real(0), AT_y  );
    Multiply( TRANSPOSE, Real(1), A, dy, Real(0), AT_dy );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = Q_x;
    Axpy( Real( 1), AT_y, rc );
    Axpy( Real(-1), z,    rc );
    Axpy( Real( 1), c,    rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,z) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = upperBound;
    DistMultiVec<Real> x_alpha(comm), z_alpha(comm), 
                       rb_alpha(comm), rc_alpha(comm);
    for( Int k=0; k<100; ++k, alpha=alpha/ctrl.stepRatio )
    {
        // x(alpha) = x + alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // z(alpha) = z + alpha dz
        // ^^^^^^^^^^^^^^^^^^^^^^^ 
        z_alpha = z;
        Axpy( alpha, dz, z_alpha );

        // mu(alpha) = x(alpha)^T z / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,z_alpha) / n;
        if( mu_alpha > (1-alpha/ctrl.psi)*mu )
        {
            if( ctrl.print && commRank == 0 )
                cout << "Armijo condition not satisfied" << endl;
            continue;
        }

        // Check 
        //    x(alpha), z(alpha) > 0, and 
        //    x_i(alpha) z_i(alpha) >= gamma mu(alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        byte locallyBalanced = true;
        for( Int iLoc=0; iLoc<nLocal; ++iLoc )
        {
            const Real xi_alpha = x_alpha.GetLocal(iLoc,0);
            const Real zi_alpha = z_alpha.GetLocal(iLoc,0);
            if( xi_alpha <= Real(0) || zi_alpha <= Real(0) )
                locallyBalanced = false;
            if( xi_alpha*zi_alpha < ctrl.gamma*mu_alpha )
                locallyBalanced = false;
        }
        const byte balanced = 
            mpi::AllReduce( locallyBalanced, mpi::BINARY_AND, comm ); 
        if( !balanced )
        {
            if( ctrl.print && commRank == 0 )
                cout << "  unbalanced entries" << endl;
            continue;
        }
        // Check ||r_b(alpha)||_2 <= Max(bTol,||r_b||_2 beta mu(alpha)/mu)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > Max(bTol,rbNrm2*ctrl.beta*mu_alpha/mu) )
        {
            if( ctrl.print && commRank == 0 )
                cout << "  r_b failure: " << rb_alphaNrm2 << " > Max("
                     << bTol << "," << rbNrm2*ctrl.beta*mu_alpha/mu << ")"
                     << endl;
            continue;
        }
        // Check ||r_c(alpha)||_2 <= Max(cTol,||r_c||_2 beta mu(alpha)/mu)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy(  alpha, Q_dx,  rc_alpha );
        Axpy(  alpha, AT_dy, rc_alpha );
        Axpy( -alpha, dz,    rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > Max(cTol,rcNrm2*ctrl.beta*mu_alpha/mu) )
        {
            if( ctrl.print && commRank == 0 )
                cout << "  r_c failure: " << rc_alphaNrm2 << " > Max("
                     << cTol << "," << rcNrm2*ctrl.beta*mu_alpha/mu << ")"
                     << endl;
        }
        else
            break;
    }
    return alpha;
}

#define PROTO(Real) \
  template Real IPFLineSearch \
  ( const Matrix<Real>& Q,  const Matrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    const Matrix<Real>& x,  const Matrix<Real>& y, \
    const Matrix<Real>& z, \
    const Matrix<Real>& dx, const Matrix<Real>& dy, \
    const Matrix<Real>& dz, \
    Real upperBound, \
    Real bTol, Real cTol, \
    const qp::IPFLineSearchCtrl<Real>& ctrl ); \
  template Real IPFLineSearch \
  ( const AbstractDistMatrix<Real>& Q,  const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& x,  const AbstractDistMatrix<Real>& y, \
    const AbstractDistMatrix<Real>& z, \
    const AbstractDistMatrix<Real>& dx, const AbstractDistMatrix<Real>& dy, \
    const AbstractDistMatrix<Real>& dz, \
    Real upperBound, \
    Real bTol, Real cTol, \
    const qp::IPFLineSearchCtrl<Real>& ctrl ); \
  template Real IPFLineSearch \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
    const Matrix<Real>& x,       const Matrix<Real>& y,  \
    const Matrix<Real>& z, \
    const Matrix<Real>& dx,      const Matrix<Real>& dy, \
    const Matrix<Real>& dz, \
    Real upperBound, \
    Real bTol, Real cTol, \
    const qp::IPFLineSearchCtrl<Real>& ctrl ); \
  template Real IPFLineSearch \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& x,     const DistMultiVec<Real>& y, \
    const DistMultiVec<Real>& z, \
    const DistMultiVec<Real>& dx,    const DistMultiVec<Real>& dy, \
    const DistMultiVec<Real>& dz, \
    Real upperBound, \
    Real bTol, Real cTol, \
    const qp::IPFLineSearchCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace direct
} // namespace qp
} // namespace El
