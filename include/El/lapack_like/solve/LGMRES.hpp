/*
   Copyright (c) 2009-2017, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_SOLVE_LGMRES_HPP
#define EL_SOLVE_LGMRES_HPP

// TODO(poulson): Add support for promotion?

namespace El {

template<typename Real>
struct LGMRESInfo
{
    Int numIts=0;
    Real relTol=1;
    bool metRequestedTol=false;
};

namespace lgmres {

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
LGMRESInfo<Base<Field>>
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
    LGMRESInfo<Real> info;

    // x := 0
    // ======
    Matrix<Field> x;
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    Matrix<Field> w;
    w = b;
    const Real origResidNorm = Nrm2( w );
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
    Matrix<Field> x0, q, V;
    while( !finished )
    {
        if( progress )
            Output("Starting GMRES iteration ",info.numIts);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );

        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        precond( w );

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

        // Run one round of GMRES
        // ======================
        for( Int j=0; j<restart; ++j )
        {
            if( progress )
                Output("Starting inner GMRES iteration ",j);
            const Int innerIndent = PushIndent();

            // w := A v_j
            // ----------
            applyA( Field(1), V(ALL,IR(j)), Field(0), w );

            // w := inv(M) w
            // -------------
            precond( w );

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
                // ^^^^^^^^^^^^^^^^^^^^
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
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                Axpy( y(i), V( ALL, IR(i) ), x );
            }

            // w := b - A x
            // ------------
            w = b;
            applyA( Field(-1), x, Field(1), w );

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
                    Output("WARNING: LGMRES did not converge");
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

} // namespace lgmres

// TODO(poulson): Add support for an initial guess
template<typename Field,class ApplyAType,class PrecondType>
LGMRESInfo<Base<Field>>
LGMRES
( const ApplyAType& applyA,
  const PrecondType& precond,
        Matrix<Field>& B,
        Base<Field> relTol,
        Int restart,
        Int maxIts,
        bool progress )
{
    EL_DEBUG_CSE
    LGMRESInfo<Base<Field>> info;
    info.relTol = 0;
    info.metRequestedTol = true;
    const Int width = B.Width();
    for( Int j=0; j<width; ++j )
    {
        auto b = B( ALL, IR(j) );
        auto singleInfo =
          lgmres::Single
          ( applyA, precond, b, relTol, restart, maxIts, progress );
        info.numIts = Max( info.numIts, singleInfo.numIts );
        info.relTol = Max( info.relTol, singleInfo.relTol );
        info.metRequestedTol =
          info.metRequestedTol && singleInfo.metRequestedTol;
    }
    return info;
}

namespace lgmres {

// TODO(poulson): Add support for an initial guess
template<typename Field,class ApplyAType,class PrecondType>
LGMRESInfo<Base<Field>>
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
    LGMRESInfo<Real> info;

    // x := 0
    // ======
    DistMultiVec<Field> x(grid);
    Zeros( x, n, 1 );

    // w := b (= b - A x_0)
    // ====================
    DistMultiVec<Field> w(grid);
    w = b;
    const Real origResidNorm = Nrm2( w );
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
    DistMultiVec<Field> x0(grid), q(grid), V(grid);
    while( !finished )
    {
        if( progress && commRank == 0 )
            Output("Starting GMRES iteration ",info.numIts);
        const Int indent = PushIndent();

        Zeros( cs, restart, 1 );
        Zeros( sn, restart, 1 );
        Zeros( H,  restart, restart );
        Zeros( V, n, restart );
        // TODO(poulson): Extend DistMultiVec so that it can be directly
        // manipulated rather than requiring access to the local Matrix and
        // staging through the temporary vector q
        auto& VLoc = V.Matrix();
        Zeros( q, n, 1 );

        // x0 := x
        // =======
        x0 = x;

        // w := inv(M) w
        // =============
        precond( w );

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
                Output("Starting inner GMRES iteration ",j);
            const Int innerIndent = PushIndent();

            // w := A v_j
            // ----------
            q.Matrix() = VLoc( ALL, IR(j) );
            applyA( Field(1), q, Field(0), w );

            // w := inv(M) w
            // -------------
            precond( w );

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
            // x := x0 + Vj y
            // ^^^^^^^^^^^^^^
            x = x0;
            for( Int i=0; i<=j; ++i )
            {
                Axpy( y(i), VLoc( ALL, IR(i) ), x.Matrix() );
            }

            // w := b - A x
            // ------------
            w = b;
            applyA( Field(-1), x, Field(1), w );

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
                    Output("WARNING: LGMRES did not converge");
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

} // namespace lgmres

// TODO(poulson): Add support for an initial guess
template<typename Field,class ApplyAType,class PrecondType>
LGMRESInfo<Base<Field>>
LGMRES
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
    LGMRESInfo<Base<Field>> info;
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
          lgmres::Single
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

#endif // ifndef EL_SOLVE_LGMRES_HPP
