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

// TODO: Add std::function<Real(Real)> to generalize sigma = (mu_aff/mu)^3
template<typename Real>
void MPC
( const Matrix<Real>& A, 
  const Matrix<Real>& b,  const Matrix<Real>& c,
  Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l,
  Real tol, Int maxIts,
  Real maxPriStepRatio, Real maxDualStepRatio, 
  bool print, System system )
{
    DEBUG_ONLY(CallStackEntry cse("lin_prog::MPC"))    

    // TODO: Check that x and s are strictly positive

    const Int m = A.Height();
    const Int n = A.Width();
    Matrix<Real> J, y, rb, rc, 
                 dsAff, dxAff, dlAff,
                 dsCC,  dxCC,  dlCC;
#ifndef RELEASE
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

        // Compute the affine search direction
        // ===================================
        if( system == LIN_PROG_FULL )
        {
            // Construct the reduced KKT system, J dl = y
            // ------------------------------------------
            FormSystem( A, b, c, s, x, l, Real(0), J, y );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            SolveSystem( m, n, J, y, dsAff, dxAff, dlAff );
        }
        else if( system == LIN_PROG_AUGMENTED )
        {
            // Construct the reduced KKT system, J dl = y
            // ------------------------------------------
            FormAugmentedSystem( A, b, c, s, x, l, Real(0), J, y );

            // Compute the proposed step from the KKT system
            // ---------------------------------------------
            SolveAugmentedSystem( s, x, Real(0), J, y, dsAff, dxAff, dlAff );
        }
        else
            LogicError("Unsupported system option");

#ifndef RELEASE
        // Sanity checks
        // =============
        dsError.Resize( n, 1 );
        for( Int i=0; i<n; ++i )
        {
            const Real xi = x.Get(i,0);
            const Real si = s.Get(i,0);
            dsError.Set( i, 0, xi*si );
        }
        Real rmuNrm2 = Nrm2( dsError );
        for( Int i=0; i<n; ++i )
        {
            const Real xi = x.Get(i,0);
            const Real si = s.Get(i,0);
            const Real dxi = dxAff.Get(i,0);
            const Real dsi = dsAff.Get(i,0);
            dsError.Update( i, 0, xi*dsi + si*dxi );
        }
        Real dsErrorNrm2 = Nrm2( dsError );

        dlError = dsAff;
        Gemv( TRANSPOSE, Real(1), A, dlAff, Real(1), dlError );
        Axpy( Real(1), rc, dlError );
        Real dlErrorNrm2 = Nrm2( dlError );

        dxError = rb;
        Gemv( NORMAL, Real(1), A, dxAff, Real(1), dxError );
        Real dxErrorNrm2 = Nrm2( dxError );

        if( print )
            std::cout << "  || dsAffError ||_2 / || r_mu ||_2 = " 
                      << dsErrorNrm2/rmuNrm2 << "\n"
                      << "  || dxAffError ||_2 / || r_c ||_2 = " 
                      << dxErrorNrm2/rcNrm2 << "\n"
                      << "  || dlAffError ||_2 / || r_b ||_2 = " 
                      << dlErrorNrm2/rbNrm2 << std::endl;
#endif

        // Compute the maximum steps which preserve positivity
        // ===================================================
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
            const Real si = s.Get(i,0);
            const Real dsAffi = dsAff.Get(i,0);
            if( dsAffi < Real(0) )
                alphaAffDual = Min(alphaAffDual,-si/dsAffi);
        }

        // Compute what the new duality measure would become
        // =================================================
        const Real mu = Dot(x,s) / n;
        // NOTE: dsCC and dxCC are used as temporaries
        dxCC = x;
        dsCC = s;
        Axpy( alphaAffPri, dxAff, dxCC );
        Axpy( alphaAffDual, dsAff, dsCC );
        const Real muAff = Dot(dxCC,dsCC) / n;

        // Compute a centrality parameter using Mehotra's formula
        // ======================================================
        // TODO: Allow the user to override this function
        const Real sigma = Pow(muAff/mu,Real(3)); 

        // Solve for the centering-corrector 
        // =================================
        // TODO
        LogicError("This routine is not yet finished");

        // TODO ...
    }
}

#define PROTO(Real) \
  template void MPC \
  ( const Matrix<Real>& A, \
    const Matrix<Real>& b,  const Matrix<Real>& c, \
    Matrix<Real>& s, Matrix<Real>& x, Matrix<Real>& l, \
    Real tol, Int maxIts, \
    Real maxPriStepRatio, Real maxDualStepRatio, \
    bool print, System system ); 

#define EL_NO_INT_PROTO
#define EL_NO_COMPLEX_PROTO
#include "El/macros/Instantiate.h"

} // namespace lin_prog
} // namespace El
