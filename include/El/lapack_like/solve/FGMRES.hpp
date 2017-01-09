/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SOLVE_FGMRES_HPP
#define EL_SOLVE_FGMRES_HPP

// The pseudocode for Flexible GMRES can be found in "Algorithm 2.2" in
//   Yousef Saad
//   "A flexible inner-outer preconditioned GMRES algorithm"
//   SIAM J. Sci. Comput., Vol. 14, No. 2, pp. 461--469, 1993.

// Add support for promotion?

namespace El {

template<typename Real>
struct FGMRESInfo
{
    Int numIts=0;
    Real relTol=1;
    bool metRequestedTol=false;
};

namespace fgmres {

// In what follows, 'applyA' should be a function of the form
//
//   void applyA
//   ( Field alpha, const Matrix<Field>& x, Field beta, Matrix<Field>& y )
//
// and overwrite y := alpha A x + beta y. However, 'precond' should have the
// form
//
//   void precond( Matrix<Field>& b )
//
// and overwrite b with an approximation of inv(A) b.
//

// TODO(poulson): Add support for an initial guess
template<typename Field,class ApplyAType,class PrecondType>
FGMRESInfo<Base<Field>>
Single
( const ApplyAType& applyA,
  const PrecondType& precond,
        Matrix<Field>& b,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<Field> Real;
    const Int n = b.Height();
    Timer iterTimer;
    FGMRESInfo<Real> info;

    // A z_j and A x_0
    const bool saveProducts = true;
    const bool time = false;

    // x := 0
    // ======
    Matrix<Field> x;
    Zeros( x, n, 1 );

    Matrix<Field> Ax0;
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
    {
        info.numIts = 0;
        info.relTol = 0;
        info.metRequestedTol = true;
        return info;
    }

    bool finished = false;
    Matrix<Real> cs;
    Matrix<Field> sn, H, t;
    Matrix<Field> x0, V, Z, AZ, q;
    while( !finished )
    {
        if( progress )
            Output("Starting FGMRES iteration ",info.numIts);
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
        if( saveProducts && info.numIts != 0 )
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
        t(0) = beta;

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
            applyA( Field(1), zj, Field(0), w );
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
                H(i,j) = Dot(vi,w);

                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H(i,j), vi, w );
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
                const Real& c = cs(i);
                const Field& s = sn(i);
                const Field sConj = Conj(s);
                const Field eta_i_j = H(i,j);
                const Field eta_ip1_j = H(i+1,j);
                H(i,  j) =  c    *eta_i_j + s*eta_ip1_j;
                H(i+1,j) = -sConj*eta_i_j + c*eta_ip1_j;
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const Field eta_j_j = H(j,j);
            const Field eta_jp1_j = delta;
            if( !limits::IsFinite(RealPart(eta_j_j))   ||
                !limits::IsFinite(ImagPart(eta_j_j))   ||
                !limits::IsFinite(RealPart(eta_jp1_j)) ||
                !limits::IsFinite(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was not finite");
            Real c;
            Field s;
            Field rho = Givens( eta_j_j, eta_jp1_j, c, s );
            if( !limits::IsFinite(c) ||
                !limits::IsFinite(RealPart(s)) ||
                !limits::IsFinite(ImagPart(s)) ||
                !limits::IsFinite(RealPart(rho)) ||
                !limits::IsFinite(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a non-finite number");
            H(j,j) = rho;
            cs(j) = c;
            sn(j) = s;
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const Field sConj = Conj(s);
            const Field tau_j = t(j);
            const Field tau_jp1 = t(j+1);
            t(j)   =  c    *tau_j + s*tau_jp1;
            t(j+1) = -sConj*tau_j + c*tau_jp1;
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
            Gemv( NORMAL, Field(1), Zj, yj, Field(1), x );

            // w := b - A x
            // ------------
            w = b;
            if( saveProducts )
            {
                // q := Ax = Ax0 + A Z_j y_j
                // ^^^^^^^^^^^^^^^^^^^^^^^^^
                q = Ax0;
                auto AZj = AZ( ALL, IR(0,j+1) );
                Gemv( NORMAL, Field(1), AZj, yj, Field(1), q );

                // w := b - A x
                // ^^^^^^^^^^^^
                w -= q;
            }
            else
            {
                applyA( Field(-1), x, Field(1), w );
            }

            if( time )
                Output("iter took ",iterTimer.Stop()," secs");

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( !limits::IsFinite(residNorm) )
                RuntimeError("Residual norm was not finite");
            info.relTol = residNorm/origResidNorm;
            ++info.numIts;
            if( info.relTol < relTol )
            {
                if( progress )
                    Output("converged with relative tolerance: ",info.relTol);
                finished = true;
                info.metRequestedTol = true;
                break;
            }
            else
            {
                if( progress )
                    Output
                    ("finished iteration ",info.numIts-1,
                     " with relative accuracy ",info.relTol);
            }
            if( info.numIts == maxIts )
            {
                if( progress )
                    Output("WARNING: FGMRES did not converge");
                finished = true;
                info.metRequestedTol = false;
                break;
            }
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return info;
}

} // namespace fgmres

// TODO(poulson): Add support for an initial guess
template<typename Field,class ApplyAType,class PrecondType>
FGMRESInfo<Base<Field>>
FGMRES
( const ApplyAType& applyA,
  const PrecondType& precond,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    EL_DEBUG_CSE
    FGMRESInfo<Base<Field>> info;
    info.relTol = 0;
    info.metRequestedTol = true;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        auto singleInfo =
          fgmres::Single
          ( applyA, precond, b, relTol, restart, maxIts, progress );
        info.numIts = Max( info.numIts, singleInfo.numIts );
        info.relTol = Max( info.relTol, singleInfo.relTol );
        info.metRequestedTol =
          info.metRequestedTol && singleInfo.metRequestedTol;
    }
    return info;
}

namespace fgmres {

// TODO(poulson): Add support for an initial guess
template<typename Field,class ApplyAType,class PrecondType>
FGMRESInfo<Base<Field>>
Single
( const ApplyAType& applyA,
  const PrecondType& precond,
        DistMultiVec<Field>& b,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( b.Width() != 1 )
          LogicError("Expected a single right-hand side");
    )
    typedef Base<Field> Real;
    const Int n = b.Height();
    const Grid& grid = b.Grid();
    const int commRank = grid.Rank();
    Timer iterTimer;
    FGMRESInfo<Real> info;

    // Avoid half of the matrix-vector products by keeping the results of
    // A z_j and A x_0
    const bool saveProducts = true;
    const bool time = false;

    // x := 0
    // ======
    DistMultiVec<Field> x(grid);
    Zeros( x, n, 1 );

    DistMultiVec<Field> Ax0(grid);
    if( saveProducts )
    {
        // A x_0 := 0
        // ==========
        Zeros( Ax0, n, 1 );
    }

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<Field> w(grid);
    w = b;
    const Real origResidNorm = Nrm2( w );
    if( progress && commRank == 0 )
        Output("origResidNorm: ",origResidNorm);
    if( origResidNorm == Real(0) )
    {
        info.numIts = 0;
        info.relTol = 0;
        info.metRequestedTol = true;
        return info;
    }

    bool finished = false;
    Matrix<Real> cs;
    Matrix<Field> sn, H, t;
    DistMultiVec<Field> x0(grid), q(grid), V(grid), Z(grid), AZ(grid);
    while( !finished )
    {
        if( progress && commRank == 0 )
            Output("Starting FGMRES iteration ",info.numIts);
        const Int indent = PushIndent();

        // x0 := x
        // =======
        x0 = x;
        if( saveProducts && info.numIts != 0 )
            Ax0 = q;

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        Zeros( Z, n, restart );
        if( saveProducts )
            Zeros( AZ, n, restart );

        // TODO(poulson): Extend DistMultiVec so that it can be directly
        // manipulated rather than requiring access to the local Matrix and
        // staging through the temporary vector q
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
        t(0) = beta;

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
            applyA( Field(1), q, Field(0), w );
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
                H(i,j) = Dot(q,w);

                // w := w - H(i,j) v_i
                // ^^^^^^^^^^^^^^^^^^^
                Axpy( -H(i,j), q, w );
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
                const Real& c = cs(i);
                const Field& s = sn(i);
                const Field sConj = Conj(s);
                const Field eta_i_j = H(i,j);
                const Field eta_ip1_j = H(i+1,j);
                H(i,  j) =  c    *eta_i_j + s*eta_ip1_j;
                H(i+1,j) = -sConj*eta_i_j + c*eta_ip1_j;
            }

            // Generate and apply a new rotation to both H and the rotated
            // beta*e_0 vector, t, then solve the minimum residual problem
            // -----------------------------------------------------------
            const Field eta_j_j = H(j,j);
            const Field eta_jp1_j = delta;
            if( !limits::IsFinite(RealPart(eta_j_j))   ||
                !limits::IsFinite(ImagPart(eta_j_j))   ||
                !limits::IsFinite(RealPart(eta_jp1_j)) ||
                !limits::IsFinite(ImagPart(eta_jp1_j)) )
                RuntimeError("Either H(j,j) or H(j+1,j) was not finite");
            Real c;
            Field s;
            Field rho = Givens( eta_j_j, eta_jp1_j, c, s );
            if( !limits::IsFinite(c) ||
                !limits::IsFinite(RealPart(s)) ||
                !limits::IsFinite(ImagPart(s)) ||
                !limits::IsFinite(RealPart(rho)) ||
                !limits::IsFinite(ImagPart(rho)) )
                RuntimeError("Givens rotation produced a non-finite number");
            H(j,j) = rho;
            cs(j) = c;
            sn(j) = s;
            // Apply the rotation to the rotated beta*e_0 vector
            // ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            const Field sConj = Conj(s);
            const Field tau_j = t(j);
            const Field tau_jp1 = t(j+1);
            t(j)   =  c    *tau_j + s*tau_jp1;
            t(j+1) = -sConj*tau_j + c*tau_jp1;
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
            Gemv( NORMAL, Field(1), ZjLoc, yj, Field(1), x.Matrix() );

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
                Gemv( NORMAL, Field(1), AZjLoc, yj, Field(1), q.Matrix() );

                // w := b - A x
                // ^^^^^^^^^^^^
                w -= q;
            }
            else
            {
                applyA( Field(-1), x, Field(1), w );
            }

            if( time && commRank == 0 )
                Output("iter took ",iterTimer.Stop()," secs");

            // Residual checks
            // ---------------
            const Real residNorm = Nrm2( w );
            if( !limits::IsFinite(residNorm) )
                RuntimeError("Residual norm was not finite");
            info.relTol = residNorm/origResidNorm;
            ++info.numIts;
            if( info.relTol < relTol )
            {
                if( progress && commRank == 0 )
                    Output("converged with relative tolerance: ",info.relTol);
                finished = true;
                info.metRequestedTol = true;
                break;
            }
            else
            {
                if( progress && commRank == 0 )
                    Output
                    ("finished iteration ",info.numIts-1,
                     " with relative accuracy ",info.relTol);
            }
            if( info.numIts == maxIts )
            {
                if( progress && commRank == 0 )
                    Output("WARNING: FGMRES did not converge");
                finished = true;
                info.metRequestedTol = false;
                break;
            }
            SetIndent( innerIndent );
        }
        SetIndent( indent );
    }
    b = x;
    return info;
}

} // namespace fgmres

// TODO(poulson): Add support for an initial guess
template<typename Field,class ApplyAType,class PrecondType>
FGMRESInfo<Base<Field>>
FGMRES
( const ApplyAType& applyA,
  const PrecondType& precond,
        DistMultiVec<Field>& B,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    EL_DEBUG_CSE
    const Int height = B.Height();
    const Int width = B.Width();
    FGMRESInfo<Base<Field>> info;
    info.relTol = 0;
    info.metRequestedTol = true;

    DistMultiVec<Field> u(B.Grid());
    Zeros( u, height, 1 );
    auto& BLoc = B.Matrix();
    auto& uLoc = u.Matrix();
    for( Int j=0; j<width; ++j )
    {
        auto bLoc = BLoc( ALL, IR(j) );
        uLoc = bLoc;
        auto singleInfo =
          fgmres::Single
          ( applyA, precond, u, relTol, restart, maxIts, progress );
        bLoc = uLoc;
        info.numIts = Max( info.numIts, singleInfo.numIts );
        info.relTol = Max( info.relTol, singleInfo.relTol );
        info.metRequestedTol =
          info.metRequestedTol && singleInfo.metRequestedTol;
    }
    return info;
}

} // namespace El

#endif // ifndef EL_SOLVE_FGMRES_HPP
