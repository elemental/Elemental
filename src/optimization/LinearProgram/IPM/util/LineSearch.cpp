/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace lin_prog {

template<typename Real>
Real IPFLineSearch
( const Matrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& s,  const Matrix<Real>& x,  const Matrix<Real>& l,
  const Matrix<Real>& ds, const Matrix<Real>& dx, const Matrix<Real>& dl,
  const IPFLineSearchCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPFLineSearch"))
    if( ctrl.gamma <= Real(0) || ctrl.gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( ctrl.beta < Real(1) )
        LogicError("beta must be at least one");
    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || l.Width() != 1 || s.Width() != 1 ||
        dx.Width() != 1 || dl.Width() != 1 || ds.Width() != 1 )
        LogicError("{b,c,x,l,s,dx,dl,ds} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-\infty" neighborhood of the central path, i.e.,
    //  (a) || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu,
    //  (b) || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu,
    //  (c) x(\alpha), s(\alpha) > 0, and, for all i,
    //  (d) x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha),
    // where 
    //    x(\alpha) = x + \alpha dx,
    //    l(\alpha) = l + \alpha dl,
    //    s(\alpha) = s + \alpha ds,
    //    r_b(\alpha) = A x(\alpha) - b, and
    //    r_c(\alpha) = A^T l(\alpha) + s(\alpha) - c,
    // and the Armijo condition,
    //   \mu(\alpha) \le (1 - \alpha/\psi) \mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    Matrix<Real> A_x, A_dx, AT_l, AT_dl, rb, rc;
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_l,  n, 1 );
    Zeros( AT_dl, n, 1 );
    Gemv( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Gemv( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Gemv( TRANSPOSE, Real(1), A, l,  Real(0), AT_l  );
    Gemv( TRANSPOSE, Real(1), A, dl, Real(0), AT_dl );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = AT_l;
    Axpy( Real(1), s, rc );
    Axpy( Real(-1), c, rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,s) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = 1;
    Matrix<Real> x_alpha, s_alpha, rb_alpha, rc_alpha;
    for( Int k=0; k<100; ++k, alpha=alpha/2 )
    {
        // x(\alpha) = x + \alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // s(\alpha) = s + \alpha ds
        // ^^^^^^^^^^^^^^^^^^^^^^^^^ 
        s_alpha = s;
        Axpy( alpha, ds, s_alpha );

        // \mu(\alpha) = x(\alpha)^T s / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,s_alpha) / n;
        if( mu_alpha > (1-alpha/ctrl.psi)*mu )
        {
            if( ctrl.print )
                std::cout << "Armijo condition not satisfied" << std::endl;
            continue;
        }

        // Check 
        //    x(\alpha), s(\alpha) > 0, and 
        //    x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        bool balanced = true;
        for( Int i=0; i<n; ++i )
        {
            const Real xi_alpha = x_alpha.Get(i,0);
            const Real si_alpha = s_alpha.Get(i,0);
            if( xi_alpha <= Real(0) || si_alpha <= Real(0) )
                balanced = false;
            if( xi_alpha*si_alpha < ctrl.gamma*mu_alpha )
                balanced = false;
        }
        if( !balanced )
        {
            if( ctrl.print )
                std::cout << "  unbalanced entries" << std::endl;
            continue;
        }
        // Check || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > rbNrm2*ctrl.beta*mu_alpha/mu )
        {
            if( ctrl.print )
                std::cout << "  r_b failure: " << rb_alphaNrm2 << " > "
                          << rbNrm2*ctrl.beta*mu_alpha/mu << std::endl;
            continue;
        }
        // Check || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy( alpha, AT_dl, rc_alpha );
        Axpy( alpha, ds, rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > rcNrm2*ctrl.beta*mu_alpha/mu )
        {
            if( ctrl.print )
                std::cout << "  r_c failure: " << rc_alphaNrm2 << " > "
                          << rcNrm2*ctrl.beta*mu_alpha/mu << std::endl;
        }
        else
            break;
    }
    return alpha;
}

template<typename Real>
Real IPFLineSearch
( const AbstractDistMatrix<Real>& APre, 
  const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c,
  const AbstractDistMatrix<Real>& s,  const AbstractDistMatrix<Real>& x,  
  const AbstractDistMatrix<Real>& l,
  const AbstractDistMatrix<Real>& ds, const AbstractDistMatrix<Real>& dx, 
  const AbstractDistMatrix<Real>& dl,
  const IPFLineSearchCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPFLineSearch"))
    if( ctrl.gamma <= Real(0) || ctrl.gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( ctrl.beta < Real(1) )
        LogicError("beta must be at least one");

    auto APtr = ReadProxy<Real,MC,MR>(&APre);
    auto& A = *APtr;

    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || l.Width() != 1 || s.Width() != 1 ||
        dx.Width() != 1 || dl.Width() != 1 || ds.Width() != 1 )
        LogicError("{b,c,x,l,s,dx,dl,ds} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const Int commRank = grid.Rank();

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-\infty" neighborhood of the central path, i.e.,
    //  (a) || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu,
    //  (b) || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu,
    //  (c) x(\alpha), s(\alpha) > 0, and, for all i,
    //  (d) x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha),
    // where 
    //    x(\alpha) = x + \alpha dx,
    //    l(\alpha) = l + \alpha dl,
    //    s(\alpha) = s + \alpha ds,
    //    r_b(\alpha) = A x(\alpha) - b, and
    //    r_c(\alpha) = A^T l(\alpha) + s(\alpha) - c,
    // and the Armijo condition,
    //   \mu(\alpha) \le (1 - \alpha/\psi) \mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    DistMatrix<Real> A_x(grid), A_dx(grid), AT_l(grid), AT_dl(grid), 
                     rb(grid), rc(grid);
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_l,  n, 1 );
    Zeros( AT_dl, n, 1 );
    Gemv( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Gemv( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Gemv( TRANSPOSE, Real(1), A, l,  Real(0), AT_l  );
    Gemv( TRANSPOSE, Real(1), A, dl, Real(0), AT_dl );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = AT_l;
    Axpy( Real(1), s, rc );
    Axpy( Real(-1), c, rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,s) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = 1;
    DistMatrix<Real> x_alpha(grid), s_alpha(grid), 
                     rb_alpha(grid), rc_alpha(grid);
    x_alpha.Align(0,0);
    s_alpha.Align(0,0);
    for( Int k=0; k<100; ++k, alpha=alpha/2 )
    {
        // x(\alpha) = x + \alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // s(\alpha) = s + \alpha ds
        // ^^^^^^^^^^^^^^^^^^^^^^^^^ 
        s_alpha = s;
        Axpy( alpha, ds, s_alpha );

        // \mu(\alpha) = x(\alpha)^T s / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,s_alpha) / n;
        if( mu_alpha > (1-alpha/ctrl.psi)*mu )
        {
            if( ctrl.print && commRank == 0 )
                std::cout << "Armijo condition not satisfied" << std::endl;
            continue;
        }

        // Check 
        //    x(\alpha), s(\alpha) > 0, and 
        //    x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        byte locallyBalanced = true;
        if( x_alpha.IsLocalCol(0) )
        {
            for( Int iLoc=0; iLoc<x_alpha.LocalHeight(); ++iLoc )
            {
                const Real xi_alpha = x_alpha.GetLocal(iLoc,0);
                const Real si_alpha = s_alpha.GetLocal(iLoc,0);
                if( xi_alpha <= Real(0) || si_alpha <= Real(0) )
                    locallyBalanced = false;
                if( xi_alpha*si_alpha < ctrl.gamma*mu_alpha )
                    locallyBalanced = false;
            }
        }
        const byte balanced =
            mpi::AllReduce( locallyBalanced, mpi::BINARY_AND, grid.VCComm() );
        if( !balanced )
        {
            if( ctrl.print && commRank == 0 )
                std::cout << "  unbalanced entries" << std::endl;
            continue;
        }
        // Check || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > rbNrm2*ctrl.beta*mu_alpha/mu )
        {
            if( ctrl.print && commRank == 0 )
                std::cout << "  r_b failure: " << rb_alphaNrm2 << " > "
                          << rbNrm2*ctrl.beta*mu_alpha/mu << std::endl;
            continue;
        }
        // Check || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy( alpha, AT_dl, rc_alpha );
        Axpy( alpha, ds, rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > rcNrm2*ctrl.beta*mu_alpha/mu )
        {
            if( ctrl.print && commRank == 0 )
                std::cout << "  r_c failure: " << rc_alphaNrm2 << " > "
                          << rcNrm2*ctrl.beta*mu_alpha/mu << std::endl;
        }
        else
            break;
    }
    return alpha;
}

template<typename Real>
Real IPFLineSearch
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  const Matrix<Real>& s,  const Matrix<Real>& x,  const Matrix<Real>& l,
  const Matrix<Real>& ds, const Matrix<Real>& dx, const Matrix<Real>& dl,
  const IPFLineSearchCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPFLineSearch"))
    if( ctrl.gamma <= Real(0) || ctrl.gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( ctrl.beta < Real(1) )
        LogicError("beta must be at least one");
    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || l.Width() != 1 || s.Width() != 1 ||
        dx.Width() != 1 || dl.Width() != 1 || ds.Width() != 1 )
        LogicError("{b,c,x,l,s,dx,dl,ds} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-\infty" neighborhood of the central path, i.e.,
    //  (a) || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu,
    //  (b) || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu,
    //  (c) x(\alpha), s(\alpha) > 0, and, for all i,
    //  (d) x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha),
    // where 
    //    x(\alpha) = x + \alpha dx,
    //    l(\alpha) = l + \alpha dl,
    //    s(\alpha) = s + \alpha ds,
    //    r_b(\alpha) = A x(\alpha) - b, and
    //    r_c(\alpha) = A^T l(\alpha) + s(\alpha) - c,
    // and the Armijo condition,
    //   \mu(\alpha) \le (1 - \alpha/\psi) \mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    Matrix<Real> A_x, A_dx, AT_l, AT_dl, rb, rc;
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_l,  n, 1 );
    Zeros( AT_dl, n, 1 );
    Multiply( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Multiply( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Multiply( TRANSPOSE, Real(1), A, l,  Real(0), AT_l  );
    Multiply( TRANSPOSE, Real(1), A, dl, Real(0), AT_dl );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = AT_l;
    Axpy( Real(1), s, rc );
    Axpy( Real(-1), c, rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,s) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = 1;
    Matrix<Real> x_alpha, s_alpha, rb_alpha, rc_alpha;
    for( Int k=0; k<100; ++k, alpha=alpha/2 )
    {
        // x(\alpha) = x + \alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // s(\alpha) = s + \alpha ds
        // ^^^^^^^^^^^^^^^^^^^^^^^^^ 
        s_alpha = s;
        Axpy( alpha, ds, s_alpha );

        // \mu(\alpha) = x(\alpha)^T s / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,s_alpha) / n;
        if( mu_alpha > (1-alpha/ctrl.psi)*mu )
        {
            if( ctrl.print )
                std::cout << "Armijo condition not satisfied" << std::endl;
            continue;
        }

        // Check 
        //    x(\alpha), s(\alpha) > 0, and 
        //    x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        bool balanced = true;
        for( Int i=0; i<n; ++i )
        {
            const Real xi_alpha = x_alpha.Get(i,0);
            const Real si_alpha = s_alpha.Get(i,0);
            if( xi_alpha <= Real(0) || si_alpha <= Real(0) )
                balanced = false;
            if( xi_alpha*si_alpha < ctrl.gamma*mu_alpha )
                balanced = false;
        }
        if( !balanced )
        {
            if( ctrl.print )
                std::cout << "  unbalanced entries" << std::endl;
            continue;
        }
        // Check || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > rbNrm2*ctrl.beta*mu_alpha/mu )
        {
            if( ctrl.print )
                std::cout << "  r_b failure: " << rb_alphaNrm2 << " > "
                          << rbNrm2*ctrl.beta*mu_alpha/mu << std::endl;
            continue;
        }
        // Check || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy( alpha, AT_dl, rc_alpha );
        Axpy( alpha, ds, rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > rcNrm2*ctrl.beta*mu_alpha/mu )
        {
            if( ctrl.print )
                std::cout << "  r_c failure: " << rc_alphaNrm2 << " > "
                          << rcNrm2*ctrl.beta*mu_alpha/mu << std::endl;
        }
        else
            break;
    }
    return alpha;
}

template<typename Real>
Real IPFLineSearch
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, 
  const DistMultiVec<Real>& c,
  const DistMultiVec<Real>& s, 
  const DistMultiVec<Real>& x, 
  const DistMultiVec<Real>& l,
  const DistMultiVec<Real>& ds, 
  const DistMultiVec<Real>& dx, 
  const DistMultiVec<Real>& dl,
  const IPFLineSearchCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPFLineSearch"))
    if( ctrl.gamma <= Real(0) || ctrl.gamma >= Real(1) )
        LogicError("gamma must be in (0,1)");
    if( ctrl.beta < Real(1) )
        LogicError("beta must be at least one");
    // TODO: Ensure communicators are congruent
    // TODO: Ensure the dimensions match
    if( b.Width() != 1 || c.Width() != 1 ||
        x.Width() != 1 || l.Width() != 1 || s.Width() != 1 ||
        dx.Width() != 1 || dl.Width() != 1 || ds.Width() != 1 )
        LogicError("{b,c,x,l,s,dx,dl,ds} must be column vectors");
    const Int m = A.Height();
    const Int n = A.Width();
    const Int nLocal = x.LocalHeight();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    // Continually halve the step size until it is small enough so that the new
    // state lies within the "-\infty" neighborhood of the central path, i.e.,
    //  (a) || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu,
    //  (b) || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu,
    //  (c) x(\alpha), s(\alpha) > 0, and, for all i,
    //  (d) x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha),
    // where 
    //    x(\alpha) = x + \alpha dx,
    //    l(\alpha) = l + \alpha dl,
    //    s(\alpha) = s + \alpha ds,
    //    r_b(\alpha) = A x(\alpha) - b, and
    //    r_c(\alpha) = A^T l(\alpha) + s(\alpha) - c,
    // and the Armijo condition,
    //   \mu(\alpha) \le (1 - \alpha/\psi) \mu,
    // is also satisfied.
    // ===============================================
    // Setup
    // -----
    DistMultiVec<Real> A_x(comm), A_dx(comm), AT_l(comm), AT_dl(comm), 
                       rb(comm), rc(comm);
    Zeros( A_x,   m, 1 );
    Zeros( A_dx,  m, 1 );
    Zeros( AT_l,  n, 1 );
    Zeros( AT_dl, n, 1 );
    Multiply( NORMAL,    Real(1), A, x,  Real(0), A_x   );
    Multiply( NORMAL,    Real(1), A, dx, Real(0), A_dx  );
    Multiply( TRANSPOSE, Real(1), A, l,  Real(0), AT_l  );
    Multiply( TRANSPOSE, Real(1), A, dl, Real(0), AT_dl );
    rb = A_x; 
    Axpy( Real(-1), b, rb );
    rc = AT_l;
    Axpy( Real(1), s, rc );
    Axpy( Real(-1), c, rc );
    const Real rbNrm2 = Nrm2( rb );
    const Real rcNrm2 = Nrm2( rc );
    const Real mu = Dot(x,s) / n;
    // Perform the line search using the cached data
    // ---------------------------------------------
    Real alpha = 1;
    DistMultiVec<Real> x_alpha(comm), s_alpha(comm), 
                       rb_alpha(comm), rc_alpha(comm);
    for( Int k=0; k<100; ++k, alpha=alpha/2 )
    {
        // x(\alpha) = x + \alpha dx
        // ^^^^^^^^^^^^^^^^^^^^^^^^^
        x_alpha = x;
        Axpy( alpha, dx, x_alpha );

        // s(\alpha) = s + \alpha ds
        // ^^^^^^^^^^^^^^^^^^^^^^^^^ 
        s_alpha = s;
        Axpy( alpha, ds, s_alpha );

        // \mu(\alpha) = x(\alpha)^T s / n
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        const Real mu_alpha = Dot(x_alpha,s_alpha) / n;
        if( mu_alpha > (1-alpha/ctrl.psi)*mu )
        {
            if( ctrl.print && commRank == 0 )
                std::cout << "Armijo condition not satisfied" << std::endl;
            continue;
        }

        // Check 
        //    x(\alpha), s(\alpha) > 0, and 
        //    x_i(\alpha) s_i(\alpha) \ge \gamma \mu(\alpha)
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
        byte locallyBalanced = true;
        for( Int iLoc=0; iLoc<nLocal; ++iLoc )
        {
            const Real xi_alpha = x_alpha.GetLocal(iLoc,0);
            const Real si_alpha = s_alpha.GetLocal(iLoc,0);
            if( xi_alpha <= Real(0) || si_alpha <= Real(0) )
                locallyBalanced = false;
            if( xi_alpha*si_alpha < ctrl.gamma*mu_alpha )
                locallyBalanced = false;
        }
        const byte balanced = 
            mpi::AllReduce( locallyBalanced, mpi::BINARY_AND, comm ); 
        if( !balanced )
        {
            if( ctrl.print && commRank == 0 )
                std::cout << "  unbalanced entries" << std::endl;
            continue;
        }
        // Check || r_b(\alpha) ||_2 \le || r_b ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rb_alpha = rb;
        Axpy( alpha, A_dx, rb_alpha );
        const Real rb_alphaNrm2 = Nrm2( rb_alpha );
        if( rb_alphaNrm2 > rbNrm2*ctrl.beta*mu_alpha/mu )
        {
            if( ctrl.print && commRank == 0 )
                std::cout << "  r_b failure: " << rb_alphaNrm2 << " > "
                          << rbNrm2*ctrl.beta*mu_alpha/mu << std::endl;
            continue;
        }
        // Check || r_c(\alpha) ||_2 \le || r_c ||_2 \beta \mu(\alpha) / \mu
        // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^ 
        rc_alpha = rc;
        Axpy( alpha, AT_dl, rc_alpha );
        Axpy( alpha, ds, rc_alpha );
        const Real rc_alphaNrm2 = Nrm2( rc_alpha );
        if( rc_alphaNrm2 > rcNrm2*ctrl.beta*mu_alpha/mu )
        {
            if( ctrl.print && commRank == 0 )
                std::cout << "  r_c failure: " << rc_alphaNrm2 << " > "
                          << rcNrm2*ctrl.beta*mu_alpha/mu << std::endl;
        }
        else
            break;
    }
    return alpha;
}

#define PROTO(Real) \
  template Real IPFLineSearch \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    const Matrix<Real>& s,  const Matrix<Real>& x,  const Matrix<Real>& l, \
    const Matrix<Real>& ds, const Matrix<Real>& dx, const Matrix<Real>& dl, \
    const IPFLineSearchCtrl<Real>& ctrl ); \
  template Real IPFLineSearch \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c, \
    const AbstractDistMatrix<Real>& s,  const AbstractDistMatrix<Real>& x, \
    const AbstractDistMatrix<Real>& l, \
    const AbstractDistMatrix<Real>& ds, const AbstractDistMatrix<Real>& dx, \
    const AbstractDistMatrix<Real>& dl, \
    const IPFLineSearchCtrl<Real>& ctrl ); \
  template Real IPFLineSearch \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    const Matrix<Real>& s,  const Matrix<Real>& x,  const Matrix<Real>& l, \
    const Matrix<Real>& ds, const Matrix<Real>& dx, const Matrix<Real>& dl, \
    const IPFLineSearchCtrl<Real>& ctrl ); \
  template Real IPFLineSearch \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c, \
    const DistMultiVec<Real>& s,  const DistMultiVec<Real>& x, \
    const DistMultiVec<Real>& l, \
    const DistMultiVec<Real>& ds, const DistMultiVec<Real>& dx, \
    const DistMultiVec<Real>& dl, \
    const IPFLineSearchCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
