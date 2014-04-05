/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRSV_LN_HPP
#define ELEM_TRSV_LN_HPP

#include ELEM_AXPY_INC
#include ELEM_ZEROS_INC
#include ELEM_GEMV_INC

namespace elem {
namespace internal {

template<typename F>
inline void
TrsvLN( UnitOrNonUnit diag, const DistMatrix<F>& L, DistMatrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrsvLN");
        if( L.Grid() != x.Grid() )
            LogicError("{L,x} must be distributed over the same grid");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( L.Width() != xLength )
            LogicError("Nonconformal TrsvLN");
    )
    const Int m = L.Height();
    const Int bsize = Blocksize();
    const Grid& g = L.Grid();

    // Matrix views 
    DistMatrix<F> L11(g), L21(g), x1(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), x1_STAR_STAR(g);

    if( x.Width() == 1 )
    {
        DistMatrix<F,MR,STAR> x1_MR_STAR(g);
        DistMatrix<F,MC,STAR> z_MC_STAR(g);

        // Views of z[MC,* ], which will store updates to x
        DistMatrix<F,MC,STAR> z1_MC_STAR(g), z2_MC_STAR(g);

        z_MC_STAR.AlignWith( L );
        Zeros( z_MC_STAR, m, 1 );

        for( Int k=0; k<m; k+=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedViewRange( L11, L, k,    k, k+nb, k+nb );
            LockedViewRange( L21, L, k+nb, k, m,    k+nb );

            ViewRange( x1, x, k, 0, k+nb, 1 );

            ViewRange( z1_MC_STAR, z_MC_STAR, k,    0, k+nb, 1 );
            ViewRange( z2_MC_STAR, z_MC_STAR, k+nb, 0, m,    1 );

            if( k != 0 )
                x1.RowSumScatterUpdate( F(1), z1_MC_STAR );

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            Trsv
            ( LOWER, NORMAL, diag,
              L11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_MR_STAR.AlignWith( L21 );
            x1_MR_STAR = x1_STAR_STAR;
            LocalGemv( NORMAL, F(-1), L21, x1_MR_STAR, F(1), z2_MC_STAR );
        }
    }
    else
    {
        DistMatrix<F,STAR,MR> x1_STAR_MR(g);
        DistMatrix<F,MC,  MR> z1(g);
        DistMatrix<F,MR,  MC> z1_MR_MC(g);
        DistMatrix<F,STAR,MC> z_STAR_MC(g);

        // Views of z[* ,MC]
        DistMatrix<F,STAR,MC> z1_STAR_MC(g), z2_STAR_MC(g);

        z_STAR_MC.AlignWith( L );
        Zeros( z_STAR_MC, 1, m );

        for( Int k=0; k<m; k+=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedViewRange( L11, L, k,    k, k+nb, k+nb );
            LockedViewRange( L21, L, k+nb, k, m,    k+nb );

            ViewRange( x1, x, 0, k, 1, k+nb );

            ViewRange( z1_STAR_MC, z_STAR_MC, 0, k,    1, k+nb );
            ViewRange( z2_STAR_MC, z_STAR_MC, 0, k+nb, 1, m    );

            if( k != 0 )
            {
                z1_MR_MC.ColSumScatterFrom( z1_STAR_MC );
                z1.AlignWith( x1 );
                z1 = z1_MR_MC;
                Axpy( F(1), z1, x1 );
            }

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            Trsv
            ( LOWER, NORMAL, diag,
              L11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MR.AlignWith( L21 );
            x1_STAR_MR = x1_STAR_STAR;
            LocalGemv( NORMAL, F(-1), L21, x1_STAR_MR, F(1), z2_STAR_MC );
        }
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRSV_LN_HPP
