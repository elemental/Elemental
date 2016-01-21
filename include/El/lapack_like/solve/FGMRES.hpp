/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SOLVE_FGMRES_HPP
#define EL_SOLVE_FGMRES_HPP

// The pseudocode for Flexible GMRES can be found in "Algorithm 2.2" in
//   Youcef Saad
//   "A flexible inner-outer preconditioned GMRES algorithm"
//   SIAM J. Sci. Comput., Vol. 14, No. 2, pp. 461--469, 1993.

// Add support for promotion?

namespace El {

namespace fgmres {

// TODO: Add support for an initial guess
template<typename F,class ApplyAType,class PrecondType>
inline Int Single
( const ApplyAType& applyA,
  const PrecondType& precond,
        Matrix<F>& b,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("fgmres::Single");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )

    // A z_j and A x_0
    const bool saveProducts = true;
    const bool time = false;

    typedef Base<F> Real;
    const Int n = b.Height();
    Timer iterTimer;

    // x := 0
    // ======
    Matrix<F> x;
    Zeros( x, n, 1 );

    Matrix<F> Ax0;
    if( saveProducts )
    {
        // A x_0 := 0
        // ==========
        Zeros( Ax0, n, 1 );
    }

    // w := b (= b - A x_0)
    // ====================
    auto w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress )
        Output("origResidNorm: ",origResidNorm);
    if( origResidNorm == Real(0) )
        return 0;

    // TODO: Constrain the maximum number of iterations

    Int iter=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    Matrix<F> x0, V, Z, AZ, q;
    while( !converged )
    {
        if( progress )
            Output("Starting FGMRES iteration ",iter);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        Zeros( Z, n, restart );
        if( saveProducts )
            Zeros( AZ, n, restart );

        // x0 := x
        // =======
        x0 = x;
        if( saveProducts && iter != 0 )
            Ax0 = q;

        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0 = V( ALL, IR(0) );
        v0 = w;
        v0 *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // ===============================
        for( Int j=0; j<restart; ++j )
        {
            if( progress )
                Output("Starting inner FGMRES iteration ",j);
            if( time )
                iterTimer.Start();
            const Int innerIndent = PushIndent();

            // z_j := inv(M) v_j
            // =================
            auto vj = V( ALL, IR(j) );
            auto zj = Z( ALL, IR(j) );
            zj = vj;
            precond( zj );

            // w := A z_j
            // ----------
            applyA( F(1), zj, F(0), w );
            if( saveProducts )
            {
                auto Azj = AZ( ALL, IR(j) );
                Azj = w;
            }

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                auto vi = V( ALL, IR(i) );
                H.Set( i, j, Dot(vi,w) );

                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), vi, w );
            }
            const Real delta = Nrm2( w );
            if( !limits::IsFinite(delta) )
                RuntimeError("Arnoldi step produced a non-finite number");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto vjp1 = V( ALL, IR(j+1) );
                vjp1 = w;
                vjp1 *= 1/delta;
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
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( !limits::IsFinite(RealPart(eta_j_j))   ||
                !limits::IsFinite(ImagPart(eta_j_j))   ||
                !limits::IsFinite(RealPart(eta_jp1_j)) ||
                !limits::IsFinite(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was not finite");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( !limits::IsFinite(c) ||
                !limits::IsFinite(RealPart(s)) || !limits::IsFinite(ImagPart(s)) ||
                !limits::IsFinite(RealPart(rho)) || !limits::IsFinite(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a non-finite number");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s);
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0);
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            auto Zj = Z( ALL, IR(0,j+1) );
            auto yj = y( IR(0,j+1), ALL );
            Gemv( NORMAL, F(1), Zj, yj, F(1), x );

            // w := b - A x
            // ------------
            w = b;
            if( saveProducts )
            {
                // q := Ax = Ax0 + A Z_j y_j
                // ^^^^^^^^^^^^^^^^^^^^^^^^^
                q = Ax0;
                auto AZj = AZ( ALL, IR(0,j+1) );
                Gemv( NORMAL, F(1), AZj, yj, F(1), q );

                // w := b - A x
                // ^^^^^^^^^^^^
                w -= q;
            }
            else
            {
                applyA( F(-1), x, F(1), w );
            }

            if( time )
                Output("iter took ",iterTimer.Stop()," secs");

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( !limits::IsFinite(residNorm) )
                RuntimeError("Residual norm was not finite");
            const Real relResidNorm = residNorm/origResidNorm;
            if( relResidNorm < relTol )
            {
                if( progress )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("FGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return iter;
}

} // namespace fgmres

// TODO: Add support for an initial guess
template<typename F,class ApplyAType,class PrecondType>
inline Int FGMRES
( const ApplyAType& applyA,
  const PrecondType& precond,
        Matrix<F>& B,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("FGMRES"))
    Int mostIts = 0;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        const Int its =
          fgmres::Single
          ( applyA, precond, b, relTol, restart, maxIts, progress );
        mostIts = Max(mostIts,its);
    }
    return mostIts;
}

namespace fgmres {

// TODO: Add support for an initial guess
template<typename F,class ApplyAType,class PrecondType>
inline Int Single
( const ApplyAType& applyA,
  const PrecondType& precond,
        DistMultiVec<F>& b,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    DEBUG_ONLY(
      CSE cse("fgmres::Single");
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    // Avoid half of the matrix-vector products by keeping the results of
    // A z_j and A x_0
    const bool saveProducts = true;
    const bool time = false;

    typedef Base<F> Real;
    const Int n = b.Height();
    mpi::Comm comm = b.Comm();
    const int commRank = mpi::Rank(comm);
    Timer iterTimer;

    // x := 0
    // ======
    DistMultiVec<F> x(comm);
    Zeros( x, n, 1 );

    DistMultiVec<F> Ax0(comm);
    if( saveProducts )
    {
        // A x_0 := 0
        // ==========
        Zeros( Ax0, n, 1 );
    }

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<F> w(comm);
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress && commRank == 0 )
        Output("origResidNorm: ",origResidNorm);
    if( origResidNorm == Real(0) )
        return 0;

    // TODO: Constrain the maximum number of iterations
    Int iter=0;
    bool converged = false;
    Matrix<Real> cs;
    Matrix<F> sn, H, t;
    DistMultiVec<F> x0(comm), q(comm), V(comm), Z(comm), AZ(comm);
    while( !converged )
    {
        if( progress && commRank == 0 )
            Output("Starting FGMRES iteration ",iter);
        const Int indent = PushIndent();

        // x0 := x
        // =======
        x0 = x;
        if( saveProducts && iter != 0 )
            Ax0 = q;

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        Zeros( Z, n, restart );
        if( saveProducts )
            Zeros( AZ, n, restart );

        // TODO: Extend DistMultiVec so that it can be directly manipulated
        //       rather than requiring access to the local Matrix and staging
        //       through the temporary vector q
        auto& VLoc = V.Matrix();
        auto& ZLoc = Z.Matrix();
        Zeros( q, n, 1 );

        // NOTE: w = b - A x already

        // beta := || w ||_2
        // =================
        const Real beta = Nrm2( w );

        // v0 := w / beta
        // ==============
        auto v0Loc = VLoc( ALL, IR(0) );
        v0Loc = w.Matrix();
        v0Loc *= 1/beta;

        // t := beta e_0
        // =============
        Zeros( t, restart+1, 1 );
        t.Set( 0, 0, beta );

        // Run one round of GMRES(restart)
        // ===============================
        for( Int j=0; j<restart; ++j )
        {
            if( progress && commRank == 0 )
                Output("Starting inner FGMRES iteration ",j);
            if( time && commRank == 0 )
                iterTimer.Start();
            const Int innerIndent = PushIndent();

            // z_j := inv(M) v_j
            // =================
            auto vjLoc = VLoc( ALL, IR(j) );
            auto zjLoc = ZLoc( ALL, IR(j) );
            q.Matrix() = vjLoc;
            precond( q );
            zjLoc = q.Matrix();

            // w := A z_j
            // ----------
            // NOTE: q currently contains z_j
            applyA( F(1), q, F(0), w );
            if( saveProducts )
            {
                auto& AZLoc = AZ.Matrix();
                auto AzjLoc = AZLoc( ALL, IR(j) );
                AzjLoc = w.LockedMatrix();
            }

            // Run the j'th step of Arnoldi
            // ----------------------------
            for( Int i=0; i<=j; ++i )
            {
                // H(i,j) := v_i' w
                // ^^^^^^^^^^^^^^^^
                q.Matrix() = VLoc( ALL, IR(i) );
                H.Set( i, j, Dot(q,w) );

                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H.Get(i,j), q, w );
            }
            const Real delta = Nrm2( w );
            if( !limits::IsFinite(delta) )
                RuntimeError("Arnoldi step produced a non-finite number");
            if( delta == Real(0) )
                restart = j+1;
            if( j+1 != restart )
            {
                // v_{j+1} := w / delta
                // ^^^^^^^^^^^^^^^^^^^^^^^^^^
                auto v_jp1Loc = VLoc( ALL, IR(j+1) );
                v_jp1Loc = w.Matrix();
                v_jp1Loc *= 1/delta;
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
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const F eta_j_j = H.Get(j,j);
            const F eta_jp1_j = delta;
            if( !limits::IsFinite(RealPart(eta_j_j))   ||
                !limits::IsFinite(ImagPart(eta_j_j))   ||
                !limits::IsFinite(RealPart(eta_jp1_j)) ||
                !limits::IsFinite(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was not finite");
            Real c;
            F s;
            F rho = lapack::Givens( eta_j_j, eta_jp1_j, &c, &s );
            if( !limits::IsFinite(c) ||
                !limits::IsFinite(RealPart(s)) || !limits::IsFinite(ImagPart(s)) ||
                !limits::IsFinite(RealPart(rho)) || !limits::IsFinite(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a non-finite number");
            H.Set( j, j, rho );
            cs.Set( j, 0, c );
            sn.Set( j, 0, s );
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const F sConj = Conj(s);
            const F tau_j = t.Get(j,0);
            const F tau_jp1 = t.Get(j+1,0);
            t.Set( j,   0,  c    *tau_j + s*tau_jp1 );
            t.Set( j+1, 0, -sConj*tau_j + c*tau_jp1 );
            // Minimize the residual
            // ^^^^^^^^^^^^^^^^^^^^^
            auto tT = t( IR(0,j+1), ALL );
            auto HTL = H( IR(0,j+1), IR(0,j+1) );
            auto y = tT;
            Trsv( UPPER, NORMAL, NON_UNIT, HTL, y );
            // x := x0 + Zj y
            // ^^^^^^^^^^^^^^
            x = x0;
            auto ZjLoc = ZLoc( ALL, IR(0,j+1) );
            auto yj = y( IR(0,j+1), ALL );
            Gemv( NORMAL, F(1), ZjLoc, yj, F(1), x.Matrix() );

            // w := b - A x
            // ------------
            w = b;
            if( saveProducts )
            {
                // q := Ax = Ax0 + A Z_j y_j
                // ^^^^^^^^^^^^^^^^^^^^^^^^^
                q = Ax0;
                const auto& AZLoc = AZ.LockedMatrix();
                auto AZjLoc = AZLoc( ALL, IR(0,j+1) );
                Gemv( NORMAL, F(1), AZjLoc, yj, F(1), q.Matrix() );

                // w := b - A x
                // ^^^^^^^^^^^^
                w -= q;
            }
            else
            {
                applyA( F(-1), x, F(1), w );
            }

            if( time && commRank == 0 )
                Output("iter took ",iterTimer.Stop()," secs");

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( !limits::IsFinite(residNorm) )
                RuntimeError("Residual norm was not finite");
            const Real relResidNorm = residNorm/origResidNorm;
            if( relResidNorm < relTol )
            {
                if( progress && commRank == 0 )
                    Output("converged with relative tolerance: ",relResidNorm);
                converged = true;
                ++iter;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    Output
                    ("finished iteration ",iter," with relResidNorm=",
                     relResidNorm);
            }
            ++iter;
            if( iter == maxIts )
                RuntimeError("FGMRES did not converge");
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return iter;
}

} // namespace fgmres

// TODO: Add support for an initial guess
template<typename F,class ApplyAType,class PrecondType>
inline Int FGMRES
( const ApplyAType& applyA,
  const PrecondType& precond,
        DistMultiVec<F>& B,
        Base<F> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    DEBUG_ONLY(CSE cse("FGMRES"))
    const Int height = B.Height();
    const Int width = B.Width();

    Int mostIts = 0;
    DistMultiVec<F> u(B.Comm());
    Zeros( u, height, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<width; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        uLoc = bLoc;
        const Int its =
          fgmres::Single
          ( applyA, precond, u, relTol, restart, maxIts, progress );
        bLoc = uLoc;
        mostIts = Max(mostIts,its);
    }
    return mostIts;
}

} // namespace El

#endif // ifndef EL_SOLVE_FGMRES_HPP
