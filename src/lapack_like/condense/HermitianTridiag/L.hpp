/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_L_HPP
#define EL_HERMITIANTRIDIAG_L_HPP

#include "./LPan.hpp"

namespace El {
namespace herm_tridiag {

// TODO(poulson): Sequential blocked implementation
template<typename F>
void L( Matrix<F>& A, Matrix<F>& phase )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Height() != A.Width() )
          LogicError("A must be square");
    )
    const Int n = A.Height();
    if( n == 0 )
    {
        phase.Resize( 0, 1 );
        return;
    }
    phase.Resize( n-1, 1 );

    Matrix<F> w21;
    for( Int k=0; k<n-1; ++k )
    {
        const Range<Int> ind1( k,   k+1 ),
                         ind2( k+1, n   );

        auto a21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );

        auto alpha21T = A( IR(k+1),     ind1 );
        auto a21B     = A( IR(k+2,END), ind1 );

        const F tau = LeftReflector( alpha21T, a21B );
        const Base<F> epsilon1 = RealPart(alpha21T(0));
        phase(k) = tau;
        alpha21T(0) = F(1);

        Zeros( w21, a21.Height(), 1 );
        Hemv( LOWER, Conj(tau), A22, a21, F(0), w21 );
        const F alpha = -Conj(tau)*Dot( w21, a21 )/F(2);
        Axpy( alpha, a21, w21 );
        Her2( LOWER, F(-1), a21, w21, A22 );
        alpha21T(0) = epsilon1;
    }
}

// TODO(poulson):
// If there is only a single MPI process, fall down to the sequential
// implementation.
template<typename F> 
void L
( ElementalMatrix<F>& APre,
  ElementalMatrix<F>& phasePre, 
  const SymvCtrl<F>& ctrl )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( APre, phasePre );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR> phaseProx( phasePre );
    auto& A = AProx.Get();
    auto& phase = phaseProx.Get();

    const Int n = A.Height();
    if( n == 0 )
    {
        phase.Resize( 0, 1 );
        return;
    }
    const Grid& g = A.Grid();
    DistMatrix<F,MD,STAR> phaseDiag(g);
    phaseDiag.SetRoot( A.DiagonalRoot(-1) );
    phaseDiag.AlignCols( A.DiagonalAlign(-1) );
    phaseDiag.Resize( n-1, 1 );

    DistMatrix<F> WPan(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), phase1_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> APan_MC_STAR(g), WPan_MC_STAR(g);
    DistMatrix<F,MR,  STAR> APan_MR_STAR(g), WPan_MR_STAR(g);

    const Int bsize = Blocksize();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k); 

        const Range<Int> ind1( k,    k+nb ),
                         indB( k,    n    ), indR( k, n ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto A22 = A( ind2, ind2 );
        auto ABR = A( indB, indR );

        const Int nbt = Min(bsize,(n-1)-k);
        auto phase1 = phaseDiag( IR(k,k+nbt), ALL );

        if( A22.Height() > 0 )
        {
            WPan.AlignWith( A11 );
            WPan.Resize( n-k, nb );
            APan_MC_STAR.AlignWith( A11 );
            APan_MC_STAR.Resize( n-k, nb );
            WPan_MC_STAR.AlignWith( A11 );
            WPan_MC_STAR.Resize( n-k, nb );
            APan_MR_STAR.AlignWith( A11 );
            APan_MR_STAR.Resize( n-k, nb );
            WPan_MR_STAR.AlignWith( A11 );
            WPan_MR_STAR.Resize( n-k, nb );

            LPan
            ( ABR, WPan, phase1,
              APan_MC_STAR, APan_MR_STAR, 
              WPan_MC_STAR, WPan_MR_STAR, ctrl );

            auto A21_MC_STAR = APan_MC_STAR( ind2-k, ind1-k );
            auto A21_MR_STAR = APan_MR_STAR( ind2-k, ind1-k );
            auto W21_MC_STAR = WPan_MC_STAR( ind2-k, ind1-k );
            auto W21_MR_STAR = WPan_MR_STAR( ind2-k, ind1-k );

            LocalTrr2k
            ( LOWER, NORMAL, ADJOINT, NORMAL, ADJOINT,
              F(-1), A21_MC_STAR, W21_MR_STAR,
              F(-1), W21_MC_STAR, A21_MR_STAR,
              F(1),  A22 );
        }
        else
        {
            A11_STAR_STAR = A11;
            phase1_STAR_STAR.Resize( nbt, 1 );
            HermitianTridiag
            ( LOWER, A11_STAR_STAR.Matrix(), phase1_STAR_STAR.Matrix() );
            A11 = A11_STAR_STAR;
            phase1 = phase1_STAR_STAR;
        }
    }
    // Redistribute from matrix-diagonal form to fully replicated
    phase = phaseDiag;
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_L_HPP
