/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef ELEM_HERMITIANTRIDIAG_USQUARE_HPP
#define ELEM_HERMITIANTRIDIAG_USQUARE_HPP

#include "./UPanSquare.hpp"

namespace elem {
namespace herm_tridiag {

template<typename F> 
void USquare( DistMatrix<F>& A, DistMatrix<F,STAR,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("herm_tridiag::USquare");
        if( A.Grid() != t.Grid() )
            LogicError("{A,t} must be distributed over the same grid");
        if( A.Height() != A.Width() )
            LogicError("A must be square");
    )
    const Grid& g = A.Grid();
    DEBUG_ONLY(
        if( g.Height() != g.Width() )
            LogicError("g must be square");
    )
    const Int n = A.Height();
    if( n == 0 )
    {
        t.Resize( 0, 1 );
        return;
    }
    DistMatrix<F,MD,STAR> tDiag(g);
    tDiag.SetRoot( A.DiagonalRoot(1) );
    tDiag.AlignCols( A.DiagonalAlign(1) );
    tDiag.Resize( A.Height()-1, 1 );

    DistMatrix<F> WPan(g);
    DistMatrix<F,STAR,STAR> A11_STAR_STAR(g), t1_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> APan_MC_STAR(g), WPan_MC_STAR(g);
    DistMatrix<F,MR,  STAR> APan_MR_STAR(g), WPan_MR_STAR(g);

    const Int bsize = Blocksize();
    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);
        auto A00 = ViewRange( A, 0, 0, k,    k    );
        auto A01 = ViewRange( A, 0, k, k,    k+nb );
        auto A11 = ViewRange( A, k, k, k+nb, k+nb );
        auto ATL = ViewRange( A, 0, 0, k+nb, k+nb );

        if( k > 0 )
        {
            auto t1 = View( tDiag, k-1, 0, nb, 1 );
            WPan.AlignWith( A01 );
            WPan.Resize( k+nb, nb );
            APan_MC_STAR.AlignWith( A00 );
            APan_MC_STAR.Resize( k+nb, nb );
            WPan_MC_STAR.AlignWith( A00 );
            WPan_MC_STAR.Resize( k+nb, nb );
            APan_MR_STAR.AlignWith( A00 );
            APan_MR_STAR.Resize( k+nb, nb );
            WPan_MR_STAR.AlignWith( A00 );
            WPan_MR_STAR.Resize( k+nb, nb );

            UPanSquare
            ( ATL, WPan, t1,
              APan_MC_STAR, APan_MR_STAR, 
              WPan_MC_STAR, WPan_MR_STAR );

            auto A01_MC_STAR = LockedViewRange( APan_MC_STAR, 0, 0, k, nb );
            auto A01_MR_STAR = LockedViewRange( APan_MR_STAR, 0, 0, k, nb );
            auto W01_MC_STAR = LockedViewRange( WPan_MC_STAR, 0, 0, k, nb );
            auto W01_MR_STAR = LockedViewRange( WPan_MR_STAR, 0, 0, k, nb );

            LocalTrr2k
            ( UPPER, ADJOINT, ADJOINT,
              F(-1), A01_MC_STAR, W01_MR_STAR,
                     W01_MC_STAR, A01_MR_STAR,
              F(1),  A00 );
        }
        else
        {
            auto t1 = View( tDiag, 0, 0, nb-1, 1 );
            A11_STAR_STAR = A11;
            t1_STAR_STAR.Resize( nb-1, 1 );
            HermitianTridiag
            ( UPPER, A11_STAR_STAR.Matrix(), t1_STAR_STAR.Matrix() );
            A11 = A11_STAR_STAR;
            t1 = t1_STAR_STAR;
        }
    }
    // Redistribute from matrix-diagonal form to fully replicated
    t = tDiag;
}

} // namespace herm_tridiag
} // namespace elem

#endif // ifndef ELEM_HERMITIANTRIDIAG_USQUARE_HPP
