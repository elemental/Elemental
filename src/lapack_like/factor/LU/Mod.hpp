/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License,
   which can be found in the LICENSE file in the root directory, or at
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_LU_MOD_HPP
#define EL_LU_MOD_HPP

namespace El {

// Begin with an LU factorization with partial pivoting,
//     A = P^T L U,
// and turn it into a partially-pivoted LU factorization of
//     A + u v',
// say
//     (A + u v') = (PNew)^T L ( U + w v'),
// w = inv(L) P u.

// Please see subsection 2.1 from
//     Peter Stange, Andreas Griewank, and Matthias Bollhofer,
//     "On the efficient update of rectangular LU factorizations subject to
//      low rank modifications"
// which discusses the technique of Schwetlick and Kielbasinski described in
// "Numerische Lineare Algebra".

// NOTE: We will make use of (at most) 2*minDim-1 swaps. Before originally
//       forming P, reserve a large number of swaps so that running out
//       does not force an arbitrary permutation via, e.g.,
//
//       Permutation P;
//       P.ReserveSwaps( n + 2*Min(m,n)-1 );
//       LU( A, P );
//       LUMod( A, P, u, v, conjugate, tau );
//

namespace lu {

template<typename F>
void RankOneMod
( Matrix<F>& A,
        Permutation& P,
  const Matrix<F>& u,
  const Matrix<F>& v,
  bool conjugate,
  Base<F> tau )
{
    EL_DEBUG_CSE
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
    P.PermuteRows( w );
    Trsv( LOWER, NORMAL, UNIT, A, w );

    // Maintain an external vector for the temporary subdiagonal of U
    Matrix<F> uSub;
    Zeros( uSub, minDim-1, 1 );

    // Reduce w to a multiple of e0
    for( Int i=minDim-2; i>=0; --i )
    {
        // Decide if we should pivot the i'th and i+1'th rows of w
        const F lambdaSub = A(i+1,i);
        const F ups_ii = A(i,i);
        const F omega_i = w(i);
        const F omega_ip1 = w(i+1);
        const Real rightTerm = Abs(lambdaSub*omega_i+omega_ip1);
        const bool pivot = ( Abs(omega_i) < tau*rightTerm );

        const Range<Int> indi( i, i+1 ),
                         indip1( i+1, i+2 ),
                         indB( i+2, m ),
                         indR( i+1, n );

        auto lBi   = A( indB,   indi   );
        auto lBip1 = A( indB,   indip1 );
        auto uiR   = A( indi,   indR   );
        auto uip1R = A( indip1, indR   );

        if( pivot )
        {
            // P := P_i P
            P.Swap( i, i+1 );

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
            A(i,  i) = gamma;
            A(i+1,i) = 0;

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
            A(i+1,i) = gamma/delta_i;
            lBi   *= F(1)/delta_i;
            lBip1 *= F(1)/delta_ip1;

            A(i,i) = eta*ups_ii*delta_i;
            Axpy( eta, uip1R, uiR );
            uiR   *= delta_i;
            uip1R *= delta_ip1;
            uSub(i) = ups_ii*delta_ip1;

            // Finally set w(i)
            w(i) = omega_ip1*delta_i;
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
            A(i+1,i) += gamma;
            Axpy(  gamma, lBip1, lBi );
            Axpy( -gamma, uiR, uip1R );
            uSub(i) = -gamma*ups_ii;
        }
    }

    // Add the modified w v' into U
    {
        auto a0 = A( IR(0), ALL );
        const F omega_0 = w(0);
        Matrix<F> vTrans;
        Transpose( v, vTrans, conjugate );
        Axpy( omega_0, vTrans, a0 );
    }

    // Transform U from upper-Hessenberg to upper-triangular form
    for( Int i=0; i<minDim-1; ++i )
    {
        // Decide if we should pivot the i'th and i+1'th rows U
        const F lambdaSub = A(i+1,i);
        const F ups_ii = A(i,i);
        const F ups_ip1i = uSub(i);
        const Real rightTerm = Abs(lambdaSub*ups_ii+ups_ip1i);
        const bool pivot = ( Abs(ups_ii) < tau*rightTerm );

        const Range<Int> indi( i, i+1 ),
                         indip1( i+1, i+2 ),
                         indB( i+2, m ),
                         indR( i+1, n );

        auto lBi   = A( indB,   indi   );
        auto lBip1 = A( indB,   indip1 );
        auto uiR   = A( indi,   indR   );
        auto uip1R = A( indip1, indR   );

        if( pivot )
        {
            // P := P_i P
            P.Swap( i, i+1 );

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
            A(i+1,i) = ups_ip1i;
            A(i,  i) = gamma;

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
            A(i+1,i) = gamma/delta_i;
            lBi   *= F(1)/delta_i;
            lBip1 *= F(1)/delta_ip1;

            A(i,i) = ups_ip1i*delta_i;
            Axpy( eta, uip1R, uiR );
            uiR   *= delta_i;
            uip1R *= delta_ip1;
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
            A(i+1,i) += gamma;
            Axpy(  gamma, lBip1, lBi );
            Axpy( -gamma, uiR, uip1R );
        }
    }
}

template<typename F>
void RankOneMod
(       AbstractDistMatrix<F>& APre,
        DistPermutation& P,
  const AbstractDistMatrix<F>& u,
  const AbstractDistMatrix<F>& v,
  bool conjugate,
  Base<F> tau )
{
    EL_DEBUG_CSE
    const Grid& g = APre.Grid();
    typedef Base<F> Real;

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int m = A.Height();
    const Int n = A.Width();
    const Int minDim = Min(m,n);

    if( minDim != m )
        LogicError("It is assumed that height(A) <= width(A)");
    if( u.Height() != m || u.Width() != 1 )
        LogicError("u is expected to be a conforming column vector");
    if( v.Height() != n || v.Width() != 1 )
        LogicError("v is expected to be a conforming column vector");
    AssertSameGrids( A, u, v );

    // w := inv(L) P u
    // TODO: Consider locally maintaining all of w to avoid unnecessarily
    //       broadcasting at every iteration.
    DistMatrix<F> w( u );
    P.PermuteRows( w );
    Trsv( LOWER, NORMAL, UNIT, A, w );

    // Maintain an external vector for the temporary subdiagonal of U
    DistMatrix<F,MD,STAR> uSub(g);
    uSub.SetRoot( A.DiagonalRoot(-1) );
    uSub.AlignCols( A.DiagonalAlign(-1) );
    Zeros( uSub, minDim-1, 1 );

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

        const Range<Int> indB( i+2, m ),
                         indR( i+1, n ),
                         indi( i, i+1 ),
                         indip1( i+1, i+2 );

        auto lBi   = A( indB,   indi   );
        auto lBip1 = A( indB,   indip1 );
        auto uiR   = A( indi,   indR   );
        auto uip1R = A( indip1, indR   );

        if( pivot )
        {
            // P := P_i P
            P.Swap( i, i+1 );

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
            lBi   *= F(1)/delta_i;
            lBip1 *= F(1)/delta_ip1;

            A.Set( i, i, eta*ups_ii*delta_i );
            Axpy( eta, uip1R, uiR );
            uiR   *= delta_i;
            uip1R *= delta_ip1;
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
        auto a0 = A( IR(0), ALL );
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

        const Range<Int> indB( i+2, m ),
                         indR( i+1, n ),
                         indi( i, i+1 ),
                         indip1( i+1, i+2 );

        auto lBi   = A( indB,   indi   );
        auto lBip1 = A( indB,   indip1 );
        auto uiR   = A( indi,   indR   );
        auto uip1R = A( indip1, indR   );

        if( pivot )
        {
            // P := P_i P
            P.Swap( i, i+1 );

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
            lBi   *= F(1)/delta_i;
            lBip1 *= F(1)/delta_ip1;

            A.Set( i, i, ups_ip1i*delta_i );
            Axpy( eta, uip1R, uiR );
            uiR   *= delta_i;
            uip1R *= delta_ip1;
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

} // namespace lu

template<typename Field>
void LUMod
( Matrix<Field>& A,
        Permutation& P,
  const Matrix<Field>& U,
  const Matrix<Field>& V,
  bool conjugate,
  Base<Field> tau )
{
    EL_DEBUG_CSE
    // TODO(poulson): Add a higher-rank implementation.
    const Int updateRank = U.Width();
    for( Int j=0; j<updateRank; ++j)
        lu::RankOneMod( A, P, U(ALL,IR(j)), V(ALL,IR(j)), conjugate, tau );
}

template<typename Field>
void LUMod
(       AbstractDistMatrix<Field>& APre,
        DistPermutation& P,
  const AbstractDistMatrix<Field>& UPre,
  const AbstractDistMatrix<Field>& VPre,
  bool conjugate,
  Base<Field> tau )
{
    EL_DEBUG_CSE

    DistMatrixReadWriteProxy<Field,Field,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    DistMatrixReadProxy<Field,Field,MC,MR> UProx( UPre );
    auto& U = UProx.GetLocked();

    DistMatrixReadProxy<Field,Field,MC,MR> VProx( VPre );
    auto& V = VProx.GetLocked();

    // TODO(poulson): Add a higher-rank implementation.
    const Int updateRank = U.Width();
    for( Int j=0; j<updateRank; ++j)
        lu::RankOneMod( A, P, U(ALL,IR(j)), V(ALL,IR(j)), conjugate, tau );
}

} // namespace El

#endif // ifndef EL_LU_MOD_HPP
