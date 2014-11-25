/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"
#include "../LinearProgram.hpp"

namespace El {
namespace lin_prog {

template<typename Real>
void IPF
( const Matrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l,
  Real tol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print, System system )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    

    // TODO: Check that x and s are strictly positive

    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> J, y, rmu, rb, rc, ds, dx, dl;
#ifndef EL_RELEASE
    Matrix<Real> dsError, dxError, dlError;
#endif
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // |c^T x - b^T l| / (1 + |c^T x|) <= tol ?
        // ----------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = Dot(b,l); 
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        const Real bNrm2 = Nrm2( b );
        rb = b;
        Gemv( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        const Real cNrm2 = Nrm2( c );
        rc = c;
        Gemv( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // Now check the pieces
        // --------------------
        if( objConv <= tol && rbConv <= tol && rcConv <= tol )
            break;
        else if( print )
            std::cout << " iter " << numIts << ":\n"
                      << "  |c^T x - b^T l| / (1 + |c^T x|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)   = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c ||_2)   = "
                      << rcConv << std::endl;

        // Compute the duality measure and r_mu = X S e - tau e
        // ====================================================
        const Real mu = Dot(x,s) / n;
        rmu.Resize( n, 1 );
        for( Int i=0; i<n; ++i )
            rmu.Set( i, 0, x.Get(i,0)*s.Get(i,0) - sigma*mu );

        if( system == FULL_KKT )
        {
            // Construct the full KKT system
            // =============================
            KKT( A, s, x, J );
            KKTRHS( rmu, rc, rb, y );

            // Compute the proposed step from the KKT system
            // =============================================
            GaussianElimination( J, y );
            ExpandKKTSolution( m, n, y, ds, dx, dl );
        }
        else if( system == AUGMENTED_KKT )
        {
            // Construct the reduced KKT system
            // ================================
            AugmentedKKT( A, s, x, J );
            AugmentedKKTRHS( x, rmu, rc, rb, y );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, y );
            ExpandAugmentedSolution( s, x, rmu, y, ds, dx, dl );
        }
        else if( system == NORMAL_KKT )
        {
            // Construct the reduced KKT system
            // ================================
            NormalKKT( A, s, x, J );
            NormalKKTRHS( A, s, x, rmu, rc, rb, dl );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, dl );
            ExpandNormalSolution( A, c, s, x, rmu, rc, dl, ds, dx );
        }

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        const Real rmuNrm2 = Nrm2( rmu );
        dsError = rmu;
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
            std::cout << "  || dsError ||_2 / (1 + || r_mu ||_2) = " 
                      << dsErrorNrm2/(1+rmuNrm2) << "\n"
                      << "  || dxError ||_2 / (1 + || r_c ||_2) = " 
                      << dxErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dlError ||_2 / (1 + || r_b ||_2) = " 
                      << dlErrorNrm2/(1+rbNrm2) << std::endl;
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
  Real tol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print, System system )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    

    // TODO: Check that x and s are strictly positive

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const Int commRank = A.Grid().Rank();

    DistMatrix<Real> J(grid), y(grid), rmu(grid), rb(grid), rc(grid), 
                     ds(grid), dx(grid), dl(grid);
#ifndef EL_RELEASE
    DistMatrix<Real> dsError(grid), dxError(grid), dlError(grid);
#endif
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // |c^T x - b^T l| / (1 + |c^T x|) <= tol ?
        // ----------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = Dot(b,l);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        const Real bNrm2 = Nrm2( b );
        rb = b;
        Gemv( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        const Real cNrm2 = Nrm2( c );
        rc = c;
        Gemv( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // Now check the pieces
        // --------------------
        if( objConv <= tol && rbConv <= tol && rcConv <= tol )
            break;
        else if( print && commRank == 0 )
            std::cout << " iter " << numIts << ":\n"
                      << "  |c^T x - b^T l| / (1 + |c^T x|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)   = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c ||_2)   = "
                      << rcConv << std::endl;

        // Compute the duality measure and r_mu = X S e - tau e
        // ====================================================
        const Real mu = Dot(x,s) / n;
        rmu.Resize( n, 1 );
        for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
            rmu.SetLocal
            ( iLoc, 0, x.GetLocal(iLoc,0)*s.GetLocal(iLoc,0) - sigma*mu );

        if( system == FULL_KKT )
        {
            // Construct the full KKT system
            // =============================
            KKT( A, s, x, J );
            KKTRHS( rmu, rc, rb, y );

            // Compute the proposed step from the KKT system
            // =============================================
            GaussianElimination( J, y );
            ExpandKKTSolution( m, n, y, ds, dx, dl );
        }
        else if( system == AUGMENTED_KKT )
        {
            // Construct the reduced KKT system
            // ================================
            AugmentedKKT( A, s, x, J );
            AugmentedKKTRHS( x, rmu, rc, rb, y );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, y );
            ExpandAugmentedSolution( s, x, rmu, y, ds, dx, dl );
        }
        else if( system == NORMAL_KKT )
        {
            // Construct the reduced KKT system
            // ================================
            NormalKKT( A, s, x, J );
            NormalKKTRHS( A, s, x, rmu, rc, rb, dl );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, dl );
            ExpandNormalSolution( A, c, s, x, rmu, rc, dl, ds, dx );
        }

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        const Real rmuNrm2 = Nrm2( rmu );
        dsError = rmu;
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
            std::cout << "  || dsError ||_2 / (1 + || r_mu ||_2) = " 
                      << dsErrorNrm2/(1+rmuNrm2) << "\n"
                      << "  || dxError ||_2 / (1 + || r_c ||_2) = " 
                      << dxErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dlError ||_2 / (1 + || r_b ||_2) = " 
                      << dlErrorNrm2/(1+rbNrm2) << std::endl;
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

// TODO: Cache the symbolic analysis
template<typename Real>
void IPF
( const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l,
  Real tol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    

    // TODO: Check that x and s are strictly positive

    const Int n = A.Width();
    SparseMatrix<Real> J;
    Matrix<Real> rmu, rc, rb, ds, dx, dl;
#ifndef EL_RELEASE
    Matrix<Real> dsError, dxError, dlError;
#endif
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // |c^T x - b^T l| / (1 + |c^T x|) <= tol ?
        // ----------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = Dot(b,l);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        const Real bNrm2 = Nrm2( b );
        rb = b;
        Multiply( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        const Real cNrm2 = Nrm2( c );
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // Now check the pieces
        // --------------------
        if( objConv <= tol && rbConv <= tol && rcConv <= tol )
            break;
        else if( print )
            std::cout << " iter " << numIts << ":\n"
                      << "  |c^T x - b^T l| / (1 + |c^T x|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)   = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c ||_2)   = "
                      << rcConv << std::endl;

        // Compute the duality measure and r_mu = X S e - tau e
        // ====================================================
        const Real mu = Dot(x,s) / n;
        rmu.Resize( n, 1 );      
        for( Int i=0; i<n; ++i )
            rmu.Set( i, 0, x.Get(i,0)*s.Get(i,0) - sigma*mu );

        // Construct the reduced KKT system, J dl = y
        // ==========================================
        NormalKKT( A, s, x, J );
        NormalKKTRHS( A, s, x, rmu, rc, rb, dl );

        // Compute the proposed step from the KKT system
        // =============================================
        //SymmetricSolve( J, dl );
        LogicError("Sequential SymmetricSolve not yet supported");
        ExpandNormalSolution( A, c, s, x, rmu, rc, dl, ds, dx );

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        const Real rmuNrm2 = Nrm2( rmu );
        dsError = rmu;
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
            std::cout << "  || dsError ||_2 / (1 + || r_mu ||_2) = " 
                      << dsErrorNrm2/(1+rmuNrm2) << "\n"
                      << "  || dxError ||_2 / (1 + || r_c ||_2) = " 
                      << dxErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dlError ||_2 / (1 + || r_b ||_2) = " 
                      << dlErrorNrm2/(1+rbNrm2) << std::endl;
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

// TODO: Cache the symbolic analysis
template<typename Real>
void IPF
( const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c,
  DistMultiVec<Real>& s, DistMultiVec<Real>& x, DistMultiVec<Real>& l,
  Real tol, Int maxIts,
  Real sigma, Real gamma, Real beta, Real psi, bool print )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::IPF"))    

    // TODO: Check that x and s are strictly positive

    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, invMap;
    DistSparseMatrix<Real> J(comm);
    DistSymmFrontTree<Real> JFrontTree;
    DistNodalMultiVec<Real> dlNodal;

    DistMultiVec<Real> rmu(comm), rb(comm), rc(comm), 
                       ds(comm), dx(comm), dl(comm);
#ifndef EL_RELEASE
    DistMultiVec<Real> dsError(comm), dxError(comm), dlError(comm);
#endif
    for( Int numIts=0; numIts<maxIts; ++numIts )
    {
        // Check for convergence
        // =====================
        // |c^T x - b^T l| / (1 + |c^T x|) <= tol ?
        // ----------------------------------------
        const Real primObj = Dot(c,x);
        const Real dualObj = Dot(b,l);
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        const Real bNrm2 = Nrm2( b );
        rb = b;
        Multiply( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c ||_2) <= tol ?
        // --------------------------------------
        const Real cNrm2 = Nrm2( c );
        rc = c;
        Multiply( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+cNrm2);
        // Now check the pieces
        // --------------------
        if( objConv <= tol && rbConv <= tol && rcConv <= tol )
            break;
        else if( print && commRank == 0 )
            std::cout << " iter " << numIts << ":\n"
                      << "  |c^T x - b^T l| / (1 + |c^T x|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)   = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c ||_2)   = "
                      << rcConv << std::endl;

        // Compute the duality measure and r_mu = X S e - tau e
        // ====================================================
        const Real mu = Dot(x,s) / n;
        rmu.Resize( n, 1 );
        for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
        {
            const Real xi = x.GetLocal(iLoc,0);
            const Real si = s.GetLocal(iLoc,0);
            rmu.SetLocal( iLoc, 0, xi*si - sigma*mu );
        }

        // Construct the reduced KKT system, J dl = y
        // ==========================================
        // NOTE: Explicit symmetry is currently required for both METIS and
        //       for the frontal tree initialization
        NormalKKT( A, s, x, J, false );
        NormalKKTRHS( A, s, x, rmu, rc, rb, dl );

        // Compute the proposed step from the KKT system
        // =============================================
        if( numIts == 0 )
        {
            NestedDissection( J.LockedDistGraph(), map, sepTree, info );
            map.FormInverse( invMap );
        }
        JFrontTree.Initialize( J, map, sepTree, info );
        LDL( info, JFrontTree, LDL_INTRAPIV_1D ); 
        dlNodal.Pull( invMap, info, dl );
        Solve( info, JFrontTree, dlNodal );
        dlNodal.Push( invMap, info, dl );
        ExpandNormalSolution( A, c, s, x, rmu, rc, dl, ds, dx );

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        const Real rmuNrm2 = Nrm2( rmu );
        dsError = rmu;
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
            std::cout << "  || dsError ||_2 / (1 + || r_mu ||_2) = " 
                      << dsErrorNrm2/(1+rmuNrm2) << "\n"
                      << "  || dxError ||_2 / (1 + || r_c ||_2) = " 
                      << dxErrorNrm2/(1+rcNrm2) << "\n"
                      << "  || dlError ||_2 / (1 + || r_b ||_2) = " 
                      << dlErrorNrm2/(1+rbNrm2) << std::endl;
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
    Real tol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print, System system ); \
  template void IPF \
  ( const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b,  const AbstractDistMatrix<Real>& c, \
    AbstractDistMatrix<Real>& s, AbstractDistMatrix<Real>& x, \
    AbstractDistMatrix<Real>& l, \
    Real tol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print, System system ); \
  template void IPF \
  ( const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l, \
    Real tol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print ); \
  template void IPF \
  ( const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,  const DistMultiVec<Real>& c, \
    DistMultiVec<Real>& s, DistMultiVec<Real>& x, DistMultiVec<Real>& l, \
    Real tol, Int maxIts, \
    Real sigma, Real gamma, Real beta, Real psi, bool print );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
