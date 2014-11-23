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
void IPF
( const Matrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l,
  Real muTol, Real rbTol, Real rcTol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    
    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> J, y, rb, rc, dx, dl, ds;
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // mu = x^T s / n
        // ---------------
        const Real mu = Dot(x,s) / n;
        // || r_b ||_2 = || A x - b ||_2
        // -----------------------------
        rb = b;
        Gemv( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        // || r_c ||_2 = || A^T l + s - c ||_2
        // -----------------------------------
        rc = c;
        Gemv( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        // Now check the pieces
        // --------------------
        if( mu <= muTol && rbNrm2 <= rbTol && rcNrm2 <= rcTol )
            break;
        else if( print )
            std::cout << " iter " << numIts << ":\n"
                      << "  || r_b ||_2  = " << rbNrm2
                      << "  || r_c ||_2  = " << rcNrm2
                      << "  mu = " << mu << std::endl;

        // Construct the reduced KKT system, J dl = y
        // ==========================================
        FormSystem( A, b, c, s, x, l, sigma*mu, J, y );

        // Compute the proposed step from the KKT system
        // =============================================
        SolveSystem( m, n, J, y, ds, dx, dl );

#ifndef RELEASE
          // Sanity checks
          // =============
          Matrix<Real> dsError, dxError, dlError;

          dsError.Resize( n, 1 );
          for( Int i=0; i<n; ++i )
          {
              const Real xi = x.Get(i,0);
              const Real si = s.Get(i,0);
              dsError.Set( i, 0, xi*si - sigma*mu );
          }
          const Real rmuNrm2 = Nrm2( dsError );
          for( Int i=0; i<n; ++i )
          {
              const Real xi = x.Get(i,0);
              const Real si = s.Get(i,0);
              const Real dxi = dx.Get(i,0);
              const Real dsi = ds.Get(i,0);
              dsError.Update( i, 0, xi*dsi + si*dxi );
          }
          const Real dsErrorNrm2 = Nrm2( dsError );

          dlError = ds;
          Gemv( TRANSPOSE, Real(1), A, dl, Real(1), dlError );
          Axpy( Real(1), rc, dlError );
          const Real dlErrorNrm2 = Nrm2( dlError );

          dxError = rb;
          Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
          const Real dxErrorNrm2 = Nrm2( dxError );
  
          if( print )
              std::cout << "  || dsError ||_2 / || r_mu ||_2 = " 
                        << dsErrorNrm2/rmuNrm2 << "\n"
                        << "  || dxError ||_2 / || r_c ||_2 = " 
                        << dxErrorNrm2/rcNrm2 << "\n"
                        << "  || dlError ||_2 / || r_b ||_2 = " 
                        << dlErrorNrm2/rbNrm2 << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha =
          IPFLineSearch
          ( A, b, c, s, x, l, ds, dx, dl, gamma, beta, psi, print );
        if( print )
            std::cout << "  alpha = " << alpha << std::endl;

        // Update the state by stepping a distance of alpha
        // ================================================
        Axpy( alpha, ds, s );
        Axpy( alpha, dx, x );
        Axpy( alpha, dl, l );
    }
}

template<typename Real>
void IPF
( const AbstractDistMatrix<Real>& A, 
  const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c,
  AbstractDistMatrix<Real>& s, AbstractDistMatrix<Real>& x, 
  AbstractDistMatrix<Real>& l,
  Real muTol, Real rbTol, Real rcTol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    
    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const Int commRank = A.Grid().Rank();

    DistMatrix<Real> J(grid), y(grid), rb(grid), rc(grid), 
                     dx(grid), dl(grid), ds(grid);
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // mu = x^T s / n
        // ---------------
        const Real mu = Dot(x,s) / n;
        // || r_b ||_2 = || A x - b ||_2
        // -----------------------------
        rb = b;
        Gemv( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        // || r_c ||_2 = || A^T l + s - c ||_2
        // -----------------------------------
        rc = c;
        Gemv( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        // Now check the pieces
        // --------------------
        if( mu <= muTol && rbNrm2 <= rbTol && rcNrm2 <= rcTol )
            break;
        else if( print && commRank == 0 )
            std::cout << " iter " << numIts << ":\n"
                      << "  || r_b ||_2  = " << rbNrm2
                      << "  || r_c ||_2  = " << rcNrm2
                      << "  mu = " << mu << std::endl;

        // Construct the reduced KKT system, J dl = y
        // ==========================================
        FormSystem( A, b, c, s, x, l, sigma*mu, J, y );

        // Compute the proposed step from the KKT system
        // =============================================
        SolveSystem( m, n, J, y, ds, dx, dl );

#ifndef RELEASE
          // Sanity checks
          // =============
          DistMatrix<Real> dsError(grid), dxError(grid), dlError(grid);

          dsError.Resize( n, 1 );
          for( Int i=0; i<n; ++i )
          {
              const Real xi = x.Get(i,0);
              const Real si = s.Get(i,0);
              dsError.Set( i, 0, xi*si - sigma*mu );
          }
          const Real rmuNrm2 = Nrm2( dsError );
          for( Int i=0; i<n; ++i )
          {
              const Real xi = x.Get(i,0);
              const Real si = s.Get(i,0);
              const Real dxi = dx.Get(i,0);
              const Real dsi = ds.Get(i,0);
              dsError.Update( i, 0, xi*dsi + si*dxi );
          }
          const Real dsErrorNrm2 = Nrm2( dsError );

          dlError = ds;
          Gemv( TRANSPOSE, Real(1), A, dl, Real(1), dlError );
          Axpy( Real(1), rc, dlError );
          const Real dlErrorNrm2 = Nrm2( dlError );

          dxError = rb;
          Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
          const Real dxErrorNrm2 = Nrm2( dxError );
  
          if( print && commRank == 0 )
              std::cout << "  || dsError ||_2 / || r_mu ||_2 = " 
                        << dsErrorNrm2/rmuNrm2 << "\n"
                        << "  || dxError ||_2 / || r_c ||_2 = " 
                        << dxErrorNrm2/rcNrm2 << "\n"
                        << "  || dlError ||_2 / || r_b ||_2 = " 
                        << dlErrorNrm2/rbNrm2 << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha =
          IPFLineSearch
          ( A, b, c, s, x, l, ds, dx, dl, gamma, beta, psi, print );
        if( print && commRank == 0 )
            std::cout << "  alpha = " << alpha << std::endl;

        // Update the state by stepping a distance of alpha
        // ================================================
        Axpy( alpha, ds, s );
        Axpy( alpha, dx, x );
        Axpy( alpha, dl, l );
    }
}

template<typename Real>
void IPF
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l,
  Real muTol, Real rbTol, Real rcTol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    
    const Int n = A.Width();
    SparseMatrix<Real> J;
    Matrix<Real> y, rb, rc, dx, dl, ds;
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // mu = x^T s / n
        // ---------------
        const Real mu = Dot(x,s) / n;
        // || r_b ||_2 = || A x - b ||_2
        // -----------------------------
        rb = b;
        Multiply( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        // || r_c ||_2 = || A^T l + s - c ||_2
        // -----------------------------------
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        // Now check the pieces
        // --------------------
        if( mu <= muTol && rbNrm2 <= rbTol && rcNrm2 <= rcTol )
            break;
        else if( print )
            std::cout << " iter " << numIts << ":\n"
                      << "  || r_b ||_2  = " << rbNrm2
                      << "  || r_c ||_2  = " << rcNrm2
                      << "  mu = " << mu << std::endl;

        // Construct the reduced KKT system, J dl = y
        // ==========================================
        FormNormalSystem( A, b, c, s, x, l, sigma*mu, J, y );

        // Compute the proposed step from the KKT system
        // =============================================
        SolveNormalSystem( A, b, c, s, x, l, sigma*mu, J, y, ds, dx, dl );

#ifndef RELEASE
          // Sanity checks
          // =============
          Matrix<Real> dsError, dxError, dlError;

          dsError.Resize( n, 1 );
          for( Int i=0; i<n; ++i )
          {
              const Real xi = x.Get(i,0);
              const Real si = s.Get(i,0);
              dsError.Set( i, 0, xi*si - sigma*mu );
          }
          const Real rmuNrm2 = Nrm2( dsError );
          for( Int i=0; i<n; ++i )
          {
              const Real xi = x.Get(i,0);
              const Real si = s.Get(i,0);
              const Real dxi = dx.Get(i,0);
              const Real dsi = ds.Get(i,0);
              dsError.Update( i, 0, xi*dsi + si*dxi );
          }
          const Real dsErrorNrm2 = Nrm2( dsError );

          dlError = ds;
          Multiply( TRANSPOSE, Real(1), A, dl, Real(1), dlError );
          Axpy( Real(1), rc, dlError );
          const Real dlErrorNrm2 = Nrm2( dlError );

          dxError = rb;
          Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
          const Real dxErrorNrm2 = Nrm2( dxError );
  
          if( print )
              std::cout << "  || dsError ||_2 / || r_mu ||_2 = " 
                        << dsErrorNrm2/rmuNrm2 << "\n"
                        << "  || dxError ||_2 / || r_c ||_2 = " 
                        << dxErrorNrm2/rcNrm2 << "\n"
                        << "  || dlError ||_2 / || r_b ||_2 = " 
                        << dlErrorNrm2/rbNrm2 << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha =
          IPFLineSearch
          ( A, b, c, s, x, l, ds, dx, dl, gamma, beta, psi, print );
        if( print )
            std::cout << "  alpha = " << alpha << std::endl;

        // Update the state by stepping a distance of alpha
        // ================================================
        Axpy( alpha, ds, s );
        Axpy( alpha, dx, x );
        Axpy( alpha, dl, l );
    }
}

template<typename Real>
void IPF
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c,
  DistMultiVec<Real>& s, DistMultiVec<Real>& x, DistMultiVec<Real>& l,
  Real muTol, Real rbTol, Real rcTol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    
    // TODO: Input checks
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    DistSparseMatrix<Real> J(comm);
    DistMultiVec<Real> y(comm), rb(comm), rc(comm), 
                       dx(comm), dl(comm), ds(comm);
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // mu = x^T s / n
        // ---------------
        const Real mu = Dot(x,s) / n;
        // || r_b ||_2 = || A x - b ||_2
        // -----------------------------
        rb = b; 
        Multiply( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        // || r_c ||_2 = || A^T l + s - c ||_2
        // -----------------------------------
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        // Now check the pieces
        // --------------------
        if( mu <= muTol && rbNrm2 <= rbTol && rcNrm2 <= rcTol )
            break;
        else if( print && commRank == 0 )
            std::cout << "  iter " << numIts << ": \n"
                      << "  || r_b ||_2 = " << rbNrm2 
                      << "  || r_c ||_2 = " << rcNrm2 
                      << "  mu = " << mu << std::endl;

        // Construct the reduced KKT system, J dl = y
        // ==========================================
        FormNormalSystem( A, b, c, s, x, l, sigma*mu, J, y );
       
        // Compute the proposed step from the KKT system
        // =============================================
        SolveNormalSystem( A, b, c, s, x, l, sigma*mu, J, y, ds, dx, dl );

#ifndef RELEASE
          // Sanity checks
          // =============
          DistMultiVec<Real> dsError(comm), dxError(comm), dlError(comm);

          dsError.Resize( n, 1 );
          for( Int iLoc=0; iLoc<dsError.LocalHeight(); ++iLoc )
          {
              const Real xi = x.GetLocal(iLoc,0);
              const Real si = s.GetLocal(iLoc,0);
              dsError.SetLocal( iLoc, 0, xi*si - sigma*mu );
          }
          const Real rmuNrm2 = Nrm2( dsError );
          for( Int iLoc=0; iLoc<dsError.LocalHeight(); ++iLoc )
          {
              const Real xi = x.GetLocal(iLoc,0);
              const Real si = s.GetLocal(iLoc,0);
              const Real dxi = dx.GetLocal(iLoc,0);
              const Real dsi = ds.GetLocal(iLoc,0);
              dsError.UpdateLocal( iLoc, 0, xi*dsi + si*dxi );
          }
          const Real dsErrorNrm2 = Nrm2( dsError );

          dlError = ds;
          Multiply( TRANSPOSE, Real(1), A, dl, Real(1), dlError );
          Axpy( Real(1), rc, dlError );
          const Real dlErrorNrm2 = Nrm2( dlError );

          dxError = rb;
          Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
          const Real dxErrorNrm2 = Nrm2( dxError );

          if( print && commRank == 0 )
              std::cout << "  || dsError ||_2 / || r_mu ||_2 = " 
                        << dsErrorNrm2/rmuNrm2 << "\n"
                        << "  || dxError ||_2 / || r_c ||_2 = " 
                        << dxErrorNrm2/rcNrm2 << "\n"
                        << "  || dlError ||_2 / || r_b ||_2 = " 
                        << dlErrorNrm2/rbNrm2 << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha = 
          IPFLineSearch 
          ( A, b, c, s, x, l, ds, dx, dl, gamma, beta, psi, print );
        if( print && commRank == 0 )
            std::cout << "  alpha = " << alpha << std::endl;

        // Update the state by stepping a distance of alpha
        // ================================================
        Axpy( alpha, ds, s );
        Axpy( alpha, dx, x );
        Axpy( alpha, dl, l );
    }
}

#define PROTO(Real) \
  template void IPF \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l, \
    Real muTol, Real rbTol, Real rcTol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print ); \
  template void IPF \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c, \
    AbstractDistMatrix<Real>& s, AbstractDistMatrix<Real>& x, \
    AbstractDistMatrix<Real>& l, \
    Real muTol, Real rbTol, Real rcTol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print ); \
  template void IPF \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l, \
    Real muTol, Real rbTol, Real rcTol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print ); \
  template void IPF \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c, \
    DistMultiVec<Real>& s, DistMultiVec<Real>& x, DistMultiVec<Real>& l, \
    Real muTol, Real rbTol, Real rcTol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
