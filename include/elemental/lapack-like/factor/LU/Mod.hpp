/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LU_MOD_HPP
#define ELEM_LU_MOD_HPP

#include ELEM_SWAP_INC
#include ELEM_TRSV_INC

#include ELEM_PERMUTEROWS_INC

#include ELEM_ZEROS_INC

namespace elem {

// Begin with an LU factorization with partial pivoting, 
//     A = P^T L U,
// and turn it into a partially-pivoted LU factorization of
//     A + u v',
// say
//     (A + u v') = P^T L ( U + w v'),
// w = inv(L) P u.

// Please see subsection 2.1 from 
//     Peter Stange, Andreas Griewank, and Matthias Bollhofer,
//     "On the efficient update of rectangular LU factorizations subject to
//      low rank modifications"
// which discusses the technique of Schwetlick and Kielbasinski described in
// "Numerische Lineare Algebra".

template<typename F>
inline void
LUMod
( Matrix<F>& A, Matrix<Int>& perm, 
  const Matrix<F>& u, const Matrix<F>& v, bool conjugate=true, 
  Base<F> tau=0.1 )
{
    DEBUG_ONLY(CallStackEntry cse("LUMod"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    if( minDim != m )
        LogicError("It is assumed that height(A) <= width(A)");
    if( u.Height() != m || u.Width() != 1 )
        LogicError("u is expected to be a conforming column vector");
    if( v.Height() != n || v.Width() != 1 )
        LogicError("v is expected to be a conforming column vector");

    // w := inv(L) P u
    auto w( u );
    PermuteRows( w, perm );
    Trsv( LOWER, NORMAL, UNIT, A, w );

    // Maintain an external vector for the temporary subdiagonal of U
    Matrix<F> uSub;
    Zeros( uSub, minDim-1, 1 );

    // Reduce w to a multiple of e0
    for( Int i=minDim-2; i>=0; --i )
    {
        // Decide if we should pivot the i'th and i+1'th rows of w
        const F lambdaSub = A.Get( i+1, i );
        const F ups_ii = A.Get(i,i); 
        const F omega_i = w.Get( i, 0 );
        const F omega_ip1 = w.Get( i+1, 0 );
        const Real rightTerm = Abs(lambdaSub*omega_i+omega_ip1);
        const bool pivot = ( Abs(omega_i) < tau*rightTerm );

        auto lBi   = ViewRange( A, i+2, i,   m,   i+1 );
        auto lBip1 = ViewRange( A, i+2, i+1, m,   i+2 );
        auto uiR   = ViewRange( A, i,   i+1, i+1, n   );
        auto uip1R = ViewRange( A, i+1, i+1, i+2, n   );

        if( pivot )
        {
            // P := P_i P
            RowSwap( perm, i, i+1 );

            // Simultaneously perform 
            //   U := P_i U and
            //   L := P_i L P_i^T
            //
            // Then update
            //     L := L T_{i,L}^{-1},
            //     U := T_{i,L} U, 
            //     w := T_{i,L} P_i w,
            // where T_{i,L} is the Gauss transform which zeros (P_i w)_{i+1}.
            // 
            // More succinctly,
            //     gamma    := w(i) / w(i+1),
            //     w(i)     := w(i+1), 
            //     w(i+1)   := 0,
            //     L(:,i)   += gamma L(:,i+1),
            //     U(i+1,:) -= gamma U(i,:).
            const F gamma = omega_i / omega_ip1;
            const F lambda_ii = F(1) + gamma*lambdaSub;
            A.Set( i,   i, gamma );
            A.Set( i+1, i, 0     );

            auto lBiCopy = lBi;
            Swap( NORMAL, lBi, lBip1 );
            Axpy( gamma, lBiCopy, lBi );

            auto uip1RCopy = uip1R;
            RowSwap( A, i, i+1 );
            Axpy( -gamma, uip1RCopy, uip1R );

            // Force L back to *unit* lower-triangular form via the transform
            //     L := L T_{i,U}^{-1} D^{-1}, 
            // where D is diagonal and responsible for forcing L(i,i) and 
            // L(i+1,i+1) back to 1. The effect on L is:
            //     eta       := L(i,i+1)/L(i,i),
            //     L(:,i+1)  -= eta L(:,i),
            //     delta_i   := L(i,i),
            //     delta_ip1 := L(i+1,i+1),
            //     L(:,i)   /= delta_i,
            //     L(:,i+1) /= delta_ip1,
            // while the effect on U is
            //     U(i,:)   += eta U(i+1,:)
            //     U(i,:)   *= delta_i,
            //     U(i+1,:) *= delta_{i+1},
            // and the effect on w is
            //     w(i) *= delta_i.
            const F eta = lambdaSub/lambda_ii;
            const F delta_i = lambda_ii;
            const F delta_ip1 = F(1) - eta*gamma;

            Axpy( -eta, lBi, lBip1 );
            A.Set( i+1, i, gamma/delta_i );
            Scale( F(1)/delta_i,   lBi   );
            Scale( F(1)/delta_ip1, lBip1 );

            A.Set( i, i, eta*ups_ii*delta_i );
            Axpy( eta, uip1R, uiR );
            Scale( delta_i,   uiR   );
            Scale( delta_ip1, uip1R );
            uSub.Set( i, 0, ups_ii*delta_ip1 );

            // Finally set w(i)
            w.Set( i, 0, omega_ip1*delta_i );
        }
        else
        {
            // Update
            //     L := L T_{i,L}^{-1},
            //     U := T_{i,L} U, 
            //     w := T_{i,L} w,
            // where T_{i,L} is the Gauss transform which zeros w_{i+1}.
            // 
            // More succinctly,
            //     gamma    := w(i+1) / w(i),
            //     L(:,i)   += gamma L(:,i+1),
            //     U(i+1,:) -= gamma U(i,:),
            //     w(i+1)   := 0.
            const F gamma = omega_ip1 / omega_i;
            A.Update( i+1, i, gamma );
            Axpy(  gamma, lBip1, lBi );
            Axpy( -gamma, uiR, uip1R );
            uSub.Set( i, 0, -gamma*ups_ii );
        }
    }

    // Add the modified w v' into U
    {
        auto a0 = View( A, 0, 0, 1, n );
        const F omega_0 = w.Get( 0, 0 ); 
        Matrix<F> vTrans;
        Transpose( v, vTrans, conjugate );
        Axpy( omega_0, vTrans, a0 );
    }

    // Transform U from upper-Hessenberg to upper-triangular form
    for( Int i=0; i<minDim-1; ++i ) 
    {
        // Decide if we should pivot the i'th and i+1'th rows U
        const F lambdaSub = A.Get(i+1,i);
        const F ups_ii = A.Get( i, i );
        const F ups_ip1i = uSub.Get( i, 0 );
        const Real rightTerm = Abs(lambdaSub*ups_ii+ups_ip1i);
        const bool pivot = ( Abs(ups_ii) < tau*rightTerm );

        auto lBi   = ViewRange( A, i+2, i,   m,   i+1 );
        auto lBip1 = ViewRange( A, i+2, i+1, m,   i+2 );
        auto uiR   = ViewRange( A, i,   i+1, i+1, n   );
        auto uip1R = ViewRange( A, i+1, i+1, i+2, n   );

        if( pivot )
        {
            // P := P_i P
            RowSwap( perm, i, i+1 );

            // Simultaneously perform 
            //   U := P_i U and
            //   L := P_i L P_i^T
            //
            // Then update
            //     L := L T_{i,L}^{-1},
            //     U := T_{i,L} U, 
            // where T_{i,L} is the Gauss transform which zeros U(i+1,i).
            // 
            // More succinctly,
            //     gamma    := U(i+1,i) / U(i,i),
            //     L(:,i)   += gamma L(:,i+1),
            //     U(i+1,:) -= gamma U(i,:).
            const F gamma = ups_ii / ups_ip1i;
            const F lambda_ii = F(1) + gamma*lambdaSub;
            A.Set( i+1, i, ups_ip1i );
            A.Set( i, i, gamma );

            auto lBiCopy = lBi;
            Swap( NORMAL, lBi, lBip1 );
            Axpy( gamma, lBiCopy, lBi );

            auto uip1RCopy = uip1R;
            RowSwap( A, i, i+1 );
            Axpy( -gamma, uip1RCopy, uip1R );

            // Force L back to *unit* lower-triangular form via the transform
            //     L := L T_{i,U}^{-1} D^{-1}, 
            // where D is diagonal and responsible for forcing L(i,i) and 
            // L(i+1,i+1) back to 1. The effect on L is:
            //     eta       := L(i,i+1)/L(i,i),
            //     L(:,i+1)  -= eta L(:,i),
            //     delta_i   := L(i,i),
            //     delta_ip1 := L(i+1,i+1),
            //     L(:,i)   /= delta_i,
            //     L(:,i+1) /= delta_ip1,
            // while the effect on U is
            //     U(i,:)   += eta U(i+1,:)
            //     U(i,:)   *= delta_i,
            //     U(i+1,:) *= delta_{i+1}.
            const F eta = lambdaSub/lambda_ii;
            const F delta_i = lambda_ii;
            const F delta_ip1 = F(1) - eta*gamma;

            Axpy( -eta, lBi, lBip1 );
            A.Set( i+1, i, gamma/delta_i );
            Scale( F(1)/delta_i,   lBi   );
            Scale( F(1)/delta_ip1, lBip1 );

            A.Set( i, i, ups_ip1i*delta_i );
            Axpy( eta, uip1R, uiR );
            Scale( delta_i,   uiR   );
            Scale( delta_ip1, uip1R );
        }
        else
        {
            // Update
            //     L := L T_{i,L}^{-1},
            //     U := T_{i,L} U, 
            // where T_{i,L} is the Gauss transform which zeros U(i+1,i).
            // 
            // More succinctly,
            //     gamma    := U(i+1,i)/ U(i,i),
            //     L(:,i)   += gamma L(:,i+1),
            //     U(i+1,:) -= gamma U(i,:).
            const F gamma = ups_ip1i / ups_ii;
            A.Update( i+1, i, gamma );
            Axpy(  gamma, lBip1, lBi );
            Axpy( -gamma, uiR, uip1R );
        }
    }
}

template<typename F,Dist UPerm>
inline void
LUMod
( DistMatrix<F>& A, DistMatrix<Int,UPerm,STAR>& perm, 
  const DistMatrix<F>& u, const DistMatrix<F>& v, bool conjugate=true, 
  Base<F> tau=0.1 )
{
    DEBUG_ONLY(CallStackEntry cse("LUMod"))
    typedef Base<F> Real;
    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);
    const Grid& g = A.Grid();
    if( minDim != m )
        LogicError("It is assumed that height(A) <= width(A)");
    if( u.Height() != m || u.Width() != 1 )
        LogicError("u is expected to be a conforming column vector");
    if( v.Height() != n || v.Width() != 1 )
        LogicError("v is expected to be a conforming column vector");
    if( A.Grid() != perm.Grid() || perm.Grid() != u.Grid() ||
        u.Grid() != v.Grid() )
        LogicError("Grids must match");

    // w := inv(L) P u
    auto w( u );
    PermuteRows( w, perm );
    Trsv( LOWER, NORMAL, UNIT, A, w );

    // Maintain an external vector for the temporary subdiagonal of U
    DistMatrix<F,MD,STAR> uSub(g);
    uSub.SetRoot( A.DiagonalRoot(-1) );
    uSub.AlignCols( A.DiagonalAlign(-1) );
    Zeros( uSub, minDim-1, 1 );

    // TODO: Consider locally maintaining all of w to avoid unnecessarily 
    //       broadcasting at every iteration.

    // Reduce w to a multiple of e0
    for( Int i=minDim-2; i>=0; --i )
    {
        // Decide if we should pivot the i'th and i+1'th rows of w
        const F lambdaSub = A.Get(i+1,i);
        const F ups_ii = A.Get(i,i); 
        const F omega_i = w.Get( i, 0 );
        const F omega_ip1 = w.Get( i+1, 0 );
        const Real rightTerm = Abs(lambdaSub*omega_i+omega_ip1);
        const bool pivot = ( Abs(omega_i) < tau*rightTerm );

        auto lBi   = ViewRange( A, i+2, i,   m,   i+1 );
        auto lBip1 = ViewRange( A, i+2, i+1, m,   i+2 );
        auto uiR   = ViewRange( A, i,   i+1, i+1, n   );
        auto uip1R = ViewRange( A, i+1, i+1, i+2, n   );

        if( pivot )
        {
            // P := P_i P
            RowSwap( perm, i, i+1 );

            // Simultaneously perform 
            //   U := P_i U and
            //   L := P_i L P_i^T
            //
            // Then update
            //     L := L T_{i,L}^{-1},
            //     U := T_{i,L} U, 
            //     w := T_{i,L} P_i w,
            // where T_{i,L} is the Gauss transform which zeros (P_i w)_{i+1}.
            // 
            // More succinctly,
            //     gamma    := w(i) / w(i+1),
            //     w(i)     := w(i+1), 
            //     w(i+1)   := 0,
            //     L(:,i)   += gamma L(:,i+1),
            //     U(i+1,:) -= gamma U(i,:).
            const F gamma = omega_i / omega_ip1;
            const F lambda_ii = F(1) + gamma*lambdaSub;
            A.Set( i,   i, gamma );
            A.Set( i+1, i, 0     );

            auto lBiCopy = lBi;
            Swap( NORMAL, lBi, lBip1 );
            Axpy( gamma, lBiCopy, lBi );

            auto uip1RCopy = uip1R;
            RowSwap( A, i, i+1 );
            Axpy( -gamma, uip1RCopy, uip1R );

            // Force L back to *unit* lower-triangular form via the transform
            //     L := L T_{i,U}^{-1} D^{-1}, 
            // where D is diagonal and responsible for forcing L(i,i) and 
            // L(i+1,i+1) back to 1. The effect on L is:
            //     eta       := L(i,i+1)/L(i,i),
            //     L(:,i+1)  -= eta L(:,i),
            //     delta_i   := L(i,i),
            //     delta_ip1 := L(i+1,i+1),
            //     L(:,i)   /= delta_i,
            //     L(:,i+1) /= delta_ip1,
            // while the effect on U is
            //     U(i,:)   += eta U(i+1,:)
            //     U(i,:)   *= delta_i,
            //     U(i+1,:) *= delta_{i+1},
            // and the effect on w is
            //     w(i) *= delta_i.
            const F eta = lambdaSub/lambda_ii;
            const F delta_i = lambda_ii;
            const F delta_ip1 = F(1) - eta*gamma;

            Axpy( -eta, lBi, lBip1 );
            A.Set( i+1, i, gamma/delta_i );
            Scale( F(1)/delta_i,   lBi   );
            Scale( F(1)/delta_ip1, lBip1 );

            A.Set( i, i, eta*ups_ii*delta_i );
            Axpy( eta, uip1R, uiR );
            Scale( delta_i,   uiR   );
            Scale( delta_ip1, uip1R );
            uSub.Set( i, 0, ups_ii*delta_ip1 );

            // Finally set w(i)
            w.Set( i, 0, omega_ip1*delta_i );
        }
        else
        {
            // Update
            //     L := L T_{i,L}^{-1},
            //     U := T_{i,L} U, 
            //     w := T_{i,L} w,
            // where T_{i,L} is the Gauss transform which zeros w_{i+1}.
            // 
            // More succinctly,
            //     gamma    := w(i+1) / w(i),
            //     L(:,i)   += gamma L(:,i+1),
            //     U(i+1,:) -= gamma U(i,:),
            //     w(i+1)   := 0.
            const F gamma = omega_ip1 / omega_i;
            A.Update( i+1, i, gamma );
            Axpy(  gamma, lBip1, lBi );
            Axpy( -gamma, uiR, uip1R );
            uSub.Set( i, 0, -gamma*ups_ii );
        }
    }

    // Add the modified w v' into U
    {
        auto a0 = View( A, 0, 0, 1, n );
        const F omega_0 = w.Get( 0, 0 ); 
        DistMatrix<F> vTrans(g);
        vTrans.AlignWith( a0 );
        Transpose( v, vTrans, conjugate );
        Axpy( omega_0, vTrans, a0 );
    }

    // Transform U from upper-Hessenberg to upper-triangular form
    for( Int i=0; i<minDim-1; ++i ) 
    {
        // Decide if we should pivot the i'th and i+1'th rows U
        const F lambdaSub = A.Get( i+1, i );
        const F ups_ii = A.Get( i, i );
        const F ups_ip1i = uSub.Get( i, 0 );
        const Real rightTerm = Abs(lambdaSub*ups_ii+ups_ip1i);
        const bool pivot = ( Abs(ups_ii) < tau*rightTerm );

        auto lBi   = ViewRange( A, i+2, i,   m,   i+1 );
        auto lBip1 = ViewRange( A, i+2, i+1, m,   i+2 );
        auto uiR   = ViewRange( A, i,   i+1, i+1, n   );
        auto uip1R = ViewRange( A, i+1, i+1, i+2, n   );

        if( pivot )
        {
            // P := P_i P
            RowSwap( perm, i, i+1 );

            // Simultaneously perform 
            //   U := P_i U and
            //   L := P_i L P_i^T
            //
            // Then update
            //     L := L T_{i,L}^{-1},
            //     U := T_{i,L} U, 
            // where T_{i,L} is the Gauss transform which zeros U(i+1,i).
            // 
            // More succinctly,
            //     gamma    := U(i+1,i) / U(i,i),
            //     L(:,i)   += gamma L(:,i+1),
            //     U(i+1,:) -= gamma U(i,:).
            const F gamma = ups_ii / ups_ip1i;
            const F lambda_ii = F(1) + gamma*lambdaSub;
            A.Set( i+1, i, ups_ip1i );
            A.Set( i, i, gamma );

            auto lBiCopy = lBi;
            Swap( NORMAL, lBi, lBip1 );
            Axpy( gamma, lBiCopy, lBi );

            auto uip1RCopy = uip1R;
            RowSwap( A, i, i+1 );
            Axpy( -gamma, uip1RCopy, uip1R );

            // Force L back to *unit* lower-triangular form via the transform
            //     L := L T_{i,U}^{-1} D^{-1}, 
            // where D is diagonal and responsible for forcing L(i,i) and 
            // L(i+1,i+1) back to 1. The effect on L is:
            //     eta       := L(i,i+1)/L(i,i),
            //     L(:,i+1)  -= eta L(:,i),
            //     delta_i   := L(i,i),
            //     delta_ip1 := L(i+1,i+1),
            //     L(:,i)   /= delta_i,
            //     L(:,i+1) /= delta_ip1,
            // while the effect on U is
            //     U(i,:)   += eta U(i+1,:)
            //     U(i,:)   *= delta_i,
            //     U(i+1,:) *= delta_{i+1}.
            const F eta = lambdaSub/lambda_ii;
            const F delta_i = lambda_ii;
            const F delta_ip1 = F(1) - eta*gamma;
            Axpy( -eta, lBi, lBip1 );
            A.Set( i+1, i, gamma/delta_i );
            Scale( F(1)/delta_i,   lBi   );
            Scale( F(1)/delta_ip1, lBip1 );

            A.Set( i, i, ups_ip1i*delta_i );
            Axpy( eta, uip1R, uiR );
            Scale( delta_i,   uiR   );
            Scale( delta_ip1, uip1R );
        }
        else
        {
            // Update
            //     L := L T_{i,L}^{-1},
            //     U := T_{i,L} U, 
            // where T_{i,L} is the Gauss transform which zeros U(i+1,i).
            // 
            // More succinctly,
            //     gamma    := U(i+1,i)/ U(i,i),
            //     L(:,i)   += gamma L(:,i+1),
            //     U(i+1,:) -= gamma U(i,:).
            const F gamma = ups_ip1i / ups_ii;
            A.Update( i+1, i, gamma );
            Axpy(  gamma, lBip1, lBi );
            Axpy( -gamma, uiR, uip1R );
        }
    }
}

} // namespace elem

#endif // ifndef ELEM_LU_MOD_HPP
