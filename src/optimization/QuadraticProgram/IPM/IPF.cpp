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
namespace quad_prog {

template<typename Real>
void IPF
( const Matrix<Real>& Q, const Matrix<Real>& A, 
  const Matrix<Real>& b, const Matrix<Real>& c,
  Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::IPF"))    

    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> J, y, rmu, rb, rc, ds, dx, dl;
#ifndef EL_RELEASE
    Matrix<Real> dsError, dxError, dlError;
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Check that no entries of x or s are non-positive
        // ================================================
        Int numNonPos_x = 0;
        for( Int i=0; i<x.Height(); ++i )
            if( x.Get(i,0) <= Real(0) )
                ++numNonPos_x;
        Int numNonPos_s = 0;
        for( Int i=0; i<s.Height(); ++i )
            if( s.Get(i,0) <= Real(0) )
                ++numNonPos_s;
        if( numNonPos_x > 0 || numNonPos_s > 0 )
            LogicError
            (numNonPos_x," entries of x were nonpositive and ",
             numNonPos_s," entries of s were nonpositive");

        // Check for convergence
        // =====================
        // |c^T x + x^T Q x - b^T l| / (1 + |c^T x + 1/2 x^T Q x|) <= tol ?
        // ----------------------------------------------------------------
        Zeros( y, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), y );
        const Real xTQx = Dot(x,y);
        const Real primObj = Dot(c,x) + xTQx/2;
        const Real dualObj = Dot(b,l) - xTQx/2; 
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        const Real bNrm2 = Nrm2( b );
        rb = b;
        Gemv( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c + Q x ||_2) <= tol ?
        // --------------------------------------------
        // NOTE: y currently contains Q x
        Axpy( Real(1), c, y );
        const Real objGradNrm2 = Nrm2( y );
        rc = y;
        Gemv( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
        const Real rcNrm2 = Nrm2( rc );
        const Real rcConv = rcNrm2 / (Real(1)+objGradNrm2);
        // Now check the pieces
        // --------------------
        if( objConv <= ctrl.tol && rbConv <= ctrl.tol && rcConv <= ctrl.tol )
            break;
        else if( ctrl.print )
            std::cout << " iter " << numIts << ":\n"
                      << "  | primObj - dualObj | / (1 + |primObj|) = "
                      << objConv << "\n"
                      << "  || r_b ||_2 / (1 + || b ||_2)           = "
                      << rbConv << "\n"
                      << "  || r_c ||_2 / (1 + || c + Q x ||_2)     = "
                      << rcConv << std::endl;

        // Raise an exception after an unacceptable number of iterations
        // =============================================================
        if( numIts == ctrl.maxIts )
            RuntimeError
            ("Maximum number of iterations (",ctrl.maxIts,") exceeded");

        // Compute the duality measure and r_mu = X S e - tau e
        // ====================================================
        const Real mu = Dot(x,s) / n;
        rmu.Resize( n, 1 );
        for( Int i=0; i<n; ++i )
            rmu.Set( i, 0, x.Get(i,0)*s.Get(i,0) - ctrl.centering*mu );

        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // =============================
            KKT( Q, A, s, x, J );
            KKTRHS( rmu, rc, rb, y );

            // Compute the proposed step from the KKT system
            // =============================================
            GaussianElimination( J, y );
            ExpandKKTSolution( m, n, y, ds, dx, dl );
        }
        else // ctrl.system == AUGMENTED_KKT
        {
            // Construct the reduced KKT system
            // ================================
            AugmentedKKT( Q, A, s, x, J );
            AugmentedKKTRHS( x, rmu, rc, rb, y );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, y );
            ExpandAugmentedSolution( s, x, rmu, y, ds, dx, dl );
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
        Hemv( LOWER, Real(-1), Q, dx, Real(1), dlError );
        Axpy( Real(1), rc, dlError );
        const Real dlErrorNrm2 = Nrm2( dlError );

        dxError = rb;
        Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        if( ctrl.print )
            std::cout << "  || dsError ||_2 / (1 + || r_mu ||_2) = " 
                      << dsErrorNrm2/(1+rmuNrm2) << "\n"
                      << "  || dxError ||_2 / (1 + || r_b ||_2)  = " 
                      << dxErrorNrm2/(1+rbNrm2) << "\n"
                      << "  || dlError ||_2 / (1 + || r_c ||_2)  = " 
                      << dlErrorNrm2/(1+rcNrm2) << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha =
          IPFLineSearch( Q, A, b, c, s, x, l, ds, dx, dl, ctrl.lineSearchCtrl );
        if( ctrl.print )
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
( const AbstractDistMatrix<Real>& QPre, const AbstractDistMatrix<Real>& APre, 
  const AbstractDistMatrix<Real>& b,    const AbstractDistMatrix<Real>& c,
  AbstractDistMatrix<Real>& sPre, AbstractDistMatrix<Real>& xPre, 
  AbstractDistMatrix<Real>& l, const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::IPF"))    

    ProxyCtrl proxCtrl;
    proxCtrl.colConstrain = true;
    proxCtrl.rowConstrain = true;
    proxCtrl.colAlign = 0;
    proxCtrl.rowAlign = 0;
    auto QPtr = ReadProxy<Real,MC,MR>(&QPre,proxCtrl);      auto& Q = *QPtr;
    auto APtr = ReadProxy<Real,MC,MR>(&APre,proxCtrl);      auto& A = *APtr;
    auto sPtr = ReadWriteProxy<Real,MC,MR>(&sPre,proxCtrl); auto& s = *sPtr;
    auto xPtr = ReadWriteProxy<Real,MC,MR>(&xPre,proxCtrl); auto& x = *xPtr;

    const Int m = A.Height();
    const Int n = A.Width();
    const Grid& grid = A.Grid();
    const Int commRank = A.Grid().Rank();

    DistMatrix<Real> J(grid), y(grid), rmu(grid), rb(grid), rc(grid), 
                     ds(grid), dx(grid), dl(grid);
    ds.AlignWith( x );
    dx.AlignWith( x );
#ifndef EL_RELEASE
    DistMatrix<Real> dsError(grid), dxError(grid), dlError(grid);
    dsError.AlignWith( ds );
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Check that no entries of x or s are non-positive
        // ================================================
        Int numNonPos_x = 0;
        if( x.IsLocalCol(0) )
            for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
                if( x.GetLocal(iLoc,0) <= Real(0) )
                    ++numNonPos_x;
        numNonPos_x = mpi::AllReduce( numNonPos_x, x.DistComm() );
        Int numNonPos_s = 0;
        if( s.IsLocalCol(0) )
            for( Int iLoc=0; iLoc<s.LocalHeight(); ++iLoc )
                if( s.GetLocal(iLoc,0) <= Real(0) )
                    ++numNonPos_s;
        numNonPos_s = mpi::AllReduce( numNonPos_s, s.DistComm() );
        if( numNonPos_x > 0 || numNonPos_s > 0 )
            LogicError
            (numNonPos_x," entries of x were nonpositive and ",
             numNonPos_s," entries of s were nonpositive");

        // Check for convergence
        // =====================
        // |c^T x + x^T Q x - b^T l| / (1 + |c^T x + 1/2 x^T Q x|) <= tol ?
        // ----------------------------------------------------------------
        Zeros( y, n, 1 );
        Hemv( LOWER, Real(1), Q, x, Real(0), y );
        const Real xTQx = Dot(x,y);
        const Real primObj = Dot(c,x) + xTQx/2;
        const Real dualObj = Dot(b,l) - xTQx/2;
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        const Real bNrm2 = Nrm2( b );
        rb = b;
        Gemv( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c + Q x ||_2) <= tol ?
        // --------------------------------------------
        // NOTE: y currently contains Q x
        Axpy( Real(1), c, y );
        const Real objGradNrm2 = Nrm2( y );
        rc = y;
        Gemv( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
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

        // Compute the duality measure and r_mu = X S e - tau e
        // ====================================================
        const Real mu = Dot(x,s) / n;
        rmu.Resize( n, 1 );
        if( rmu.IsLocalCol(0) )
            for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
                rmu.SetLocal
                ( iLoc, 0, 
                  x.GetLocal(iLoc,0)*s.GetLocal(iLoc,0) - ctrl.centering*mu );

        if( ctrl.system == FULL_KKT )
        {
            // Construct the full KKT system
            // =============================
            KKT( Q, A, s, x, J );
            KKTRHS( rmu, rc, rb, y );

            // Compute the proposed step from the KKT system
            // =============================================
            GaussianElimination( J, y );
            ExpandKKTSolution( m, n, y, ds, dx, dl );
        }
        else // ctrl.system == AUGMENTED_KKT
        {
            // Construct the reduced KKT system
            // ================================
            AugmentedKKT( Q, A, s, x, J );
            AugmentedKKTRHS( x, rmu, rc, rb, y );

            // Compute the proposed step from the KKT system
            // =============================================
            SymmetricSolve( LOWER, NORMAL, J, y );
            ExpandAugmentedSolution( s, x, rmu, y, ds, dx, dl );
        }

#ifndef EL_RELEASE
        // Sanity checks
        // =============
        const Real rmuNrm2 = Nrm2( rmu );
        // TODO: Find a more convenient syntax for expressing this operation
        dsError = rmu;
        if( dsError.IsLocalCol(0) )
        {
            for( Int iLoc=0; iLoc<dsError.LocalHeight(); ++iLoc )
            {
                const Real xi = x.GetLocal(iLoc,0);
                const Real si = s.GetLocal(iLoc,0);
                const Real dxi = dx.GetLocal(iLoc,0);
                const Real dsi = ds.GetLocal(iLoc,0);
                dsError.UpdateLocal( iLoc, 0, xi*dsi + si*dxi );
            }
        }
        const Real dsErrorNrm2 = Nrm2( dsError );

        dlError = ds;
        Gemv( TRANSPOSE, Real(1), A, dl, Real(1), dlError );
        Hemv( LOWER, Real(-1), Q, dx, Real(1), dlError );
        Axpy( Real(1), rc, dlError );
        const Real dlErrorNrm2 = Nrm2( dlError );

        dxError = rb;
        Gemv( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        if( ctrl.print && commRank == 0 )
            std::cout << "  || dsError ||_2 / (1 + || r_mu ||_2) = " 
                      << dsErrorNrm2/(1+rmuNrm2) << "\n"
                      << "  || dxError ||_2 / (1 + || r_b ||_2)  = " 
                      << dxErrorNrm2/(1+rbNrm2) << "\n"
                      << "  || dlError ||_2 / (1 + || r_c ||_2)  = " 
                      << dlErrorNrm2/(1+rcNrm2) << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha =
          IPFLineSearch( Q, A, b, c, s, x, l, ds, dx, dl, ctrl.lineSearchCtrl );
        if( ctrl.print && commRank == 0 )
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
( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, 
  const Matrix<Real>& b,       const Matrix<Real>& c,
  Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::IPF"))    
    LogicError("Sequential sparse-direct solvers not yet supported");
}

template<typename Real>
void IPF
( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, 
  const DistMultiVec<Real>& b, const DistMultiVec<Real>& c,
  DistMultiVec<Real>& s, DistMultiVec<Real>& x, DistMultiVec<Real>& l,
  const IPFCtrl<Real>& ctrl )
{
    DEBUG_ONLY(CallStackEntry cse("quad_prog::IPF"))    

    const Int m = A.Height();
    const Int n = A.Width();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);
    const Real epsilon = lapack::MachineEpsilon<Real>();

    DistSymmInfo info;
    DistSeparatorTree sepTree;
    DistMap map, invMap;
    DistSparseMatrix<Real> J(comm);
    DistSymmFrontTree<Real> JFrontTree;

    DistMultiVec<Real> y(comm),
                       rmu(comm), rb(comm), rc(comm), 
                       ds(comm),  dx(comm), dl(comm);
    DistNodalMultiVec<Real> yNodal;

    DistMultiVec<Real> regCand(comm), reg(comm);
    DistNodalMultiVec<Real> regCandNodal, regNodal;

#ifndef EL_RELEASE
    DistMultiVec<Real> dsError(comm), dxError(comm), dlError(comm);
#endif
    for( Int numIts=0; ; ++numIts )
    {
        // Check that no entries of x or s are non-positive
        // ================================================
        Int numNonPos_x = 0;
        for( Int iLoc=0; iLoc<x.LocalHeight(); ++iLoc )
            if( x.GetLocal(iLoc,0) <= Real(0) )
                ++numNonPos_x;
        numNonPos_x = mpi::AllReduce( numNonPos_x, comm );
        Int numNonPos_s = 0;
        for( Int iLoc=0; iLoc<s.LocalHeight(); ++iLoc )
            if( s.GetLocal(iLoc,0) <= Real(0) )
                ++numNonPos_s;
        numNonPos_s = mpi::AllReduce( numNonPos_s, comm );
        if( numNonPos_x > 0 || numNonPos_s > 0 )
            LogicError
            (numNonPos_x," entries of x were nonpositive and ",
             numNonPos_s," entries of s were nonpositive");

        // Check for convergence
        // =====================
        // |c^T x + x^T Q x - b^T l| / (1 + |c^T x + 1/2 x^T Q x|) <= tol ?
        // ----------------------------------------------------------------
        Zeros( y, n, 1 );
        // NOTE: This assumes that Q is explicitly Hermitian
        Multiply( NORMAL, Real(1), Q, x, Real(0), y );
        const Real xTQx = Dot(x,y);
        const Real primObj = Dot(c,x) + xTQx/2;
        const Real dualObj = Dot(b,l) - xTQx/2;
        const Real objConv = Abs(primObj-dualObj) / (Real(1)+Abs(primObj));
        // || r_b ||_2 / (1 + || b ||_2) <= tol ?
        // --------------------------------------
        const Real bNrm2 = Nrm2( b );
        rb = b;
        Multiply( NORMAL, Real(1), A, x, Real(-1), rb );
        const Real rbNrm2 = Nrm2( rb );
        const Real rbConv = rbNrm2 / (Real(1)+bNrm2);
        // || r_c ||_2 / (1 + || c + Q x ||_2) <= tol ?
        // --------------------------------------------
        // NOTE: y currently contains Q x
        Axpy( Real(1), c, y );
        const Real objGradNrm2 = Nrm2( y );
        rc = y;
        Multiply( TRANSPOSE, Real(1), A, l, Real(-1), rc );
        Axpy( Real(1), s, rc );
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

        // Compute the duality measure and r_mu = X S e - tau e
        // ====================================================
        const Real mu = Dot(x,s) / n;
        rmu.Resize( n, 1 );
        for( Int iLoc=0; iLoc<rmu.LocalHeight(); ++iLoc )
        {
            const Real xi = x.GetLocal(iLoc,0);
            const Real si = s.GetLocal(iLoc,0);
            rmu.SetLocal( iLoc, 0, xi*si - ctrl.centering*mu );
        }

        // Compute the search direction
        // ============================
        const Real minReductionFactor = 2;
        const Int maxRefineIts = 10;
        // ctrl.system == AUGMENTED_KKT
        {
            // Construct the "normal" KKT system
            // ---------------------------------
            // TODO: Add default regularization
            AugmentedKKT( Q, A, s, x, J, false );
            AugmentedKKTRHS( x, rmu, rc, rb, y );
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
            // TODO: Iterative refinement
            /*
            SolveWithIterativeRefinement
            ( J, invMap, info, JFrontTree, y, 
              minReductionFactor, maxRefineIts );
            */
            yNodal.Pull( invMap, info, y );
            Solve( info, JFrontTree, yNodal );
            yNodal.Push( invMap, info, y );
            ExpandAugmentedSolution( s, x, rmu, y, ds, dx, dl );
        }

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
        Multiply( NORMAL, Real(-1), Q, dx, Real(1), dlError );
        Axpy( Real(1), rc, dlError );
        const Real dlErrorNrm2 = Nrm2( dlError );

        dxError = rb;
        Multiply( NORMAL, Real(1), A, dx, Real(1), dxError );
        const Real dxErrorNrm2 = Nrm2( dxError );

        if( ctrl.print && commRank == 0 )
            std::cout << "  || dsError ||_2 / (1 + || r_mu ||_2) = " 
                      << dsErrorNrm2/(1+rmuNrm2) << "\n"
                      << "  || dxError ||_2 / (1 + || r_b ||_2)  = " 
                      << dxErrorNrm2/(1+rbNrm2) << "\n"
                      << "  || dlError ||_2 / (1 + || r_c ||_2)  = " 
                      << dlErrorNrm2/(1+rcNrm2) << std::endl;
#endif

        // Decide on the step length
        // =========================
        const Real alpha = 
          IPFLineSearch( Q, A, b, c, s, x, l, ds, dx, dl, ctrl.lineSearchCtrl );
        if( ctrl.print && commRank == 0 )
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
  ( const Matrix<Real>& Q, const Matrix<Real>& A, \
    const Matrix<Real>& b, const Matrix<Real>& c, \
    Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const AbstractDistMatrix<Real>& Q, const AbstractDistMatrix<Real>& A, \
    const AbstractDistMatrix<Real>& b, const AbstractDistMatrix<Real>& c, \
    AbstractDistMatrix<Real>& s, AbstractDistMatrix<Real>& x, \
    AbstractDistMatrix<Real>& l, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const SparseMatrix<Real>& Q, const SparseMatrix<Real>& A, \
    const Matrix<Real>& b,       const Matrix<Real>& c, \
    Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l, \
    const IPFCtrl<Real>& ctrl ); \
  template void IPF \
  ( const DistSparseMatrix<Real>& Q, const DistSparseMatrix<Real>& A, \
    const DistMultiVec<Real>& b,     const DistMultiVec<Real>& c, \
    DistMultiVec<Real>& s, DistMultiVec<Real>& x, DistMultiVec<Real>& l, \
    const IPFCtrl<Real>& ctrl );

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace quad_prog
} // namespace El
