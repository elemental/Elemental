/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_LOWER_BLOCKED_SQUARE_HPP
#define EL_HERMITIANTRIDIAG_LOWER_BLOCKED_SQUARE_HPP

#include "./LowerPanelSquare.hpp"

namespace El {
namespace herm_tridiag {

template<typename F> 
void LowerBlockedSquare
( AbstractDistMatrix<F>& APre,
  AbstractDistMatrix<F>& tPre,
  const SymvCtrl<F>& ctrl )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( APre, tPre );
      if( APre.Height() != APre.Width() )
          LogicError("A must be square");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    DistMatrixWriteProxy<F,F,STAR,STAR> tProx( tPre );
    auto& A = AProx.Get();
    auto& t = tProx.Get();

    const Grid& g = A.Grid();
    EL_DEBUG_ONLY(
      if( g.Height() != g.Width() )
          LogicError("The process grid must be square");
    )
    const Int n = A.Height();
    if( n == 0 )
    {
        t.Resize( 0, 1 );
        return;
    }
    DistMatrix<F,MD,STAR> tDiag(g);
    tDiag.SetRoot( A.DiagonalRoot(-1) );
    tDiag.AlignCols( A.DiagonalAlign(-1) );
    tDiag.Resize( A.Height()-1, 1 );

    DistMatrix<F> WPan(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), t1_STAR_STAR(g);
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
        auto t1 = tDiag( IR(k,k+nbt), ALL );
            
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

            herm_tridiag::LowerPanelSquare
            ( ABR, WPan, t1,
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
              F(1), A22 );
        }
        else
        {
            A11_STAR_STAR = A11;
            t1_STAR_STAR.Resize( nbt, 1 );
            HermitianTridiag
            ( LOWER, A11_STAR_STAR.Matrix(), t1_STAR_STAR.Matrix() );
            A11 = A11_STAR_STAR;
            t1 = t1_STAR_STAR;
        }
    }
    // Redistribute from matrix-diagonal form to fully replicated
    t = tDiag;
}

} // namespace herm_tridiag
} // namespace El

#endif // ifndef EL_HERMITIANTRIDIAG_LOWER_BLOCKED_SQUARE_HPP
