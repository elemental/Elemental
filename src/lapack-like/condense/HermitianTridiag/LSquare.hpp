/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_HERMITIANTRIDIAG_LSQUARE_HPP
#define EL_HERMITIANTRIDIAG_LSQUARE_HPP

#include "./LPanSquare.hpp"

namespace El {
namespace herm_tridiag {

template<typename F> 
void LSquare( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("herm_tridiag::LSquare");
        if( A.Grid() != t.Grid() )
            LogicError("{A,t} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Grid& g = A.Grid();
    DEBUG_ONLY(
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
        auto A11 = ViewRange( A, k,    k,    k+nb, k+nb );
        auto A21 = ViewRange( A, k+nb, k,    n,    k+nb );
        auto A22 = ViewRange( A, k+nb, k+nb, n,    n    );
        auto ABR = ViewRange( A, k,    k,    n,    n    );
        const Int nbt = Min(bsize,(n-1)-k);
        auto t1 = View( tDiag, k, 0, nbt, 1 );
            
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

            herm_tridiag::LPanSquare
            ( ABR, WPan, t1,
              APan_MC_STAR, APan_MR_STAR, 
              WPan_MC_STAR, WPan_MR_STAR );

            auto A21_MC_STAR = LockedViewRange( APan_MC_STAR, nb, 0, n-k, nb );
            auto A21_MR_STAR = LockedViewRange( APan_MR_STAR, nb, 0, n-k, nb );
            auto W21_MC_STAR = LockedViewRange( WPan_MC_STAR, nb, 0, n-k, nb );
            auto W21_MR_STAR = LockedViewRange( WPan_MR_STAR, nb, 0, n-k, nb );

            LocalTrr2k
            ( LOWER, ADJOINT, ADJOINT,
              F(-1), A21_MC_STAR, W21_MR_STAR,
                     W21_MC_STAR, A21_MR_STAR,
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

#endif // ifndef EL_HERMITIANTRIDIAG_LSQUARE_HPP
