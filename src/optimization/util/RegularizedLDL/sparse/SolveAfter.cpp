/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#include "El.hpp"

namespace El {
namespace reg_ldl {

// TODO: Switch to returning the relative residual of the refined solution

template<typename F>
Int RegularizedSolveAfter
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmInfo& info,
  const DistSymmFrontTree<F>& AFact,       DistMultiVec<F>& b,
  Base<F> minReductionFactor,              Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_ldl::RegularizedSolveAfter"))
    mpi::Comm comm = b.Comm();
    const Int commRank = mpi::Rank(comm);

    DistMultiVec<F> bOrig(comm);
    bOrig = b;

    // Compute the initial guess
    // =========================
    DistMultiVec<F> x(comm);
    DistNodalMultiVec<F> xNodal;
    xNodal.Pull( invMap, info, b );
    Solve( info, AFact, xNodal );
    xNodal.Push( invMap, info, x );

    Int refineIt = 0;
    if( maxRefineIts > 0 )
    {
        DistMultiVec<F> dx(comm), xCand(comm), xRegScaled(comm);
        xRegScaled = x;
        DiagonalScale( NORMAL, reg, xRegScaled );
        Axpy( F(1), xRegScaled, b );
        Multiply( NORMAL, F(-1), A, x, F(1), b );
        Base<F> errorNorm = Nrm2( b );
        if( progress && commRank == 0 )
            std::cout << "    original error norm: " << errorNorm << std::endl;
        for( ; refineIt<maxRefineIts; ++refineIt )
        {
            // Compute the proposed update to the solution
            // -------------------------------------------
            xNodal.Pull( invMap, info, b );
            Solve( info, AFact, xNodal );
            xNodal.Push( invMap, info, dx );
            xCand = x;
            Axpy( F(1), dx, xCand );

            // If the proposed update lowers the residual, accept it
            // -----------------------------------------------------
            b = bOrig;
            xRegScaled = xCand;
            DiagonalScale( NORMAL, reg, xRegScaled );
            Axpy( F(-1), xRegScaled, b );
            Multiply( NORMAL, F(-1), A, xCand, F(1), b );
            Base<F> newErrorNorm = Nrm2( b );
            if( progress && commRank == 0 )
                std::cout << "    reduced by factor " 
                          << errorNorm/newErrorNorm << std::endl;
            if( minReductionFactor*newErrorNorm < errorNorm )
            {
                x = xCand;
                errorNorm = newErrorNorm;
            }
            else if( newErrorNorm < errorNorm )
            {
                x = xCand;
                errorNorm = newErrorNorm;
                break;
            }
            else
                break;
        }
    }
    // Store the final result
    // ======================
    b = x;
    return refineIt;
}

template<typename F>
Int LGMRESSolveAfter
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmInfo& info,
  const DistSymmFrontTree<F>& AFact,       DistMultiVec<F>& b,
  Base<F> minReductionFactor,              Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_ldl::LGMRESSolveAfter"))
    typedef Base<F> Real;
    const Int n = A.Height();
    mpi::Comm comm = A.Comm();
    const Int commRank = mpi::Rank(comm);

    // TODO: Make these exposed control parameters
    const Int k = 10; // restart parameter for GMRES(k)
    const Real relTol = Pow(lapack::MachineEpsilon<Real>(),Real(0.5)); 

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm); 
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( origResidNorm == Real(0) )
        return 0;

    Int iter=0;
    Int maxLargeRefines=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, z;
    DistMultiVec<F> x0(comm), q(comm), V(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            std::cout << "  Starting GMRES iteration " << iter << std::endl;
        Zeros( cs, k, 1 );
        Zeros( sn, k, 1 );
        Zeros( H,  k, k );
        Zeros( V, n, k );
        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        const Int nLoc = V.LocalHeight();
        Zeros( q, n, 1 );
        
        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        Int refineIts = RegularizedSolveAfter
        ( A, reg, invMap, info, AFact, w, 
          minReductionFactor, maxRefineIts, progress );
        maxLargeRefines = Max( refineIts, maxLargeRefines );

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( IR(0,nLoc), IR(0,1) );
        v0Loc = w.Matrix();
        Scale( Real(1)/beta, v0Loc ); 

        // z := beta e_0
        // =============
        Zeros( z, k+1, 1 );
        z.Set( 0, 0, beta );

        // Run one round of GMRES(k)
        // =========================
        for( Int j=0; j<k; ++j )
        {
            // w := A v_j
            // ----------
            q.Matrix() = VLoc( IR(0,nLoc), IR(j,j+1) );
            Multiply( NORMAL, F(1), A, q, F(0), w );

            // w := inv(M) w
            // -------------
            Int refineIts = RegularizedSolveAfter
            ( A, reg, invMap, info, AFact, w, 
              minReductionFactor, maxRefineIts, progress );
            maxLargeRefines = Max( refineIts, maxLargeRefines );

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( IR(0,nLoc), IR(i,i+1) );
                H.Set( i, j, Dot(q,w) ); 
              
                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( std::isnan(delta) )
                RuntimeError("Arnoldi step produced a NaN");
            if( delta == Real(0) )
                RuntimeError("Lucky breakdowns not yet handled");
            if( j+1 != k )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( IR(0,nLoc), IR(j+1,j+2) );
                v_jp1Loc = w.Matrix();
                Scale( Real(1)/delta, v_jp1Loc );
            }

            // Apply existing rotations to the new column of H
            // -----------------------------------------------
            for( Int i=0; i<j; ++i )
            {
                const Real c = cs.Get(i,0);
                const F s = sn.Get(i,0);
                const F sConj = Conj(s);
                const F eta_i_j = H.Get(i,j);
                const F eta_ip1_j = H.Get(i+1,j);
                H.Set( i,   j,  c    *eta_i_j + s*eta_ip1_j );
                H.Set( i+1, j, -sConj*eta_i_j + c*eta_ip1_j );
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, z, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( std::isnan(RealPart(eta_j_j))   || 
                std::isnan(ImagPart(eta_j_j))   ||
                std::isnan(RealPart(eta_jp1_j)) || 
                std::isnan(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was NaN");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( std::isnan(c) || 
                std::isnan(RealPart(s)) || std::isnan(ImagPart(s)) || 
                std::isnan(RealPart(rho)) || std::isnan(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a NaN");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s); 
            const F zeta_j = z.Get(j,0);
            const F zeta_jp1 = z.Get(j+1,0); 
            z.Set( j,   0,  c    *zeta_j + s*zeta_jp1 );
            z.Set( j+1, 0, -sConj*zeta_j + c*zeta_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto zT = z( IR(0,j+1), IR(0,1) );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = zT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                const F eta_i = y.Get(i,0);
                Axpy( eta_i, VLoc( IR(0,nLoc), IR(i,i+1) ), x.Matrix() );
            }

            // w := b - A x
            // ------------
            w = b;
            Multiply( NORMAL, F(-1), A, x, F(1), w );

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( std::isnan(residNorm) )
                RuntimeError("Residual norm was NaN");
            const Real relResidNorm = residNorm/origResidNorm; 
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    std::cout << "  converged with relative tolerance: "
                              << relResidNorm << std::endl;
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    std::cout << "  finished iteration " << iter << " with "
                              << "relResidNorm=" << relResidNorm << std::endl;
            }
            ++iter;
        }
    }
    b = x;
    return maxLargeRefines;
}

// TODO: Add FGMRES (and possibly RGMRES)

template<typename F>
Int SolveAfter
( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg,
  const DistMap& invMap,             const DistSymmInfo& info,
  const DistSymmFrontTree<F>& AFact,       DistMultiVec<F>& b,
  Base<F> minReductionFactor,              Int maxRefineIts,
  bool progress )
{
    DEBUG_ONLY(CallStackEntry cse("reg_ldl::SolveAfter"))
    return LGMRESSolveAfter
    ( A, reg, invMap, info, AFact, b, 
      minReductionFactor, maxRefineIts, progress );
}

#define PROTO(F) \
  template Int SolveAfter \
  ( const DistSparseMatrix<F>& A,      const DistMultiVec<Base<F>>& reg, \
    const DistMap& invMap,             const DistSymmInfo& info, \
    const DistSymmFrontTree<F>& AFact,       DistMultiVec<F>& b, \
    Base<F> minReductionFactor,              Int maxRefineIts, \
    bool progress );

#define EL_NO_INT_PROTO
#include "El/macros/Instantiate.h"

} // namespace reg_ldl
} // namespace El
