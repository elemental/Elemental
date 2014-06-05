/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

#include EL_ZEROS_INC

namespace El {
namespace trsv {

template<typename F>
inline void
LT
( Orientation orientation, UnitOrNonUnit diag, 
  const DistMatrix<F>& L, DistMatrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("trsv::LT");
        if( L.Grid() != x.Grid() )
            LogicError("{L,x} must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected a (conjugate-)transpose option");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( L.Width() != xLength )
            LogicError("Nonconformal");
    )
    const Int m = L.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    const Grid& g = L.Grid();

    // Matrix views 
    DistMatrix<F> L10(g), L11(g), x1(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g), x1_STAR_STAR(g);

    if( x.Width() == 1 )
    {
        DistMatrix<F,MC,STAR> x1_MC_STAR(g);
        DistMatrix<F,MR,MC  > z1_MR_MC(g);
        DistMatrix<F,MR,STAR> z_MR_STAR(g);

        // Views of z[MR,* ]
        DistMatrix<F,MR,STAR> z0_MR_STAR(g), z1_MR_STAR(g);

        z_MR_STAR.AlignWith( L );
        Zeros( z_MR_STAR, m, 1 );

        for( Int k=kLast; k>=0; k-=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedViewRange( L10, L, k, 0, k+nb, k    );
            LockedViewRange( L11, L, k, k, k+nb, k+nb );

            ViewRange( x1, x, k, 0, k+nb, 1 );

            ViewRange( z0_MR_STAR, z_MR_STAR, 0, 0, k,    1 );
            ViewRange( z1_MR_STAR, z_MR_STAR, k, 0, k+nb, 1 );

            if( k+nb != m )
            {
                z1_MR_MC.RowSumScatterFrom( z1_MR_STAR );
                Axpy( F(1), z1_MR_MC, x1 );
            }

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            Trsv
            ( LOWER, orientation, diag,
              L11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_MC_STAR.AlignWith( L10 );
            x1_MC_STAR = x1_STAR_STAR;
            LocalGemv( orientation, F(-1), L10, x1_MC_STAR, F(1), z0_MR_STAR );
        }
    }
    else
    {
        DistMatrix<F,STAR,MC> x1_STAR_MC(g);
        DistMatrix<F,STAR,MR> z_STAR_MR(g);

        // Views of z[* ,MR], which will store updates to x
        DistMatrix<F,STAR,MR> z0_STAR_MR(g), z1_STAR_MR(g);

        z_STAR_MR.AlignWith( L );
        Zeros( z_STAR_MR, 1, m );

        for( Int k=kLast; k>=0; k-=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedViewRange( L10, L, k, 0, k+nb, k    );
            LockedViewRange( L11, L, k, k, k+nb, k+nb );

            ViewRange( x1, x, 0, k, 1, k+nb );

            ViewRange( z0_STAR_MR, z_STAR_MR, 0, 0, 1, k    );
            ViewRange( z1_STAR_MR, z_STAR_MR, 0, k, 1, k+nb );

            if( k+nb != m )
                x1.ColSumScatterUpdate( F(1), z1_STAR_MR );

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            Trsv
            ( LOWER, orientation, diag,
              L11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MC.AlignWith( L10 );
            x1_STAR_MC = x1_STAR_STAR;
            LocalGemv( orientation, F(-1), L10, x1_STAR_MC, F(1), z0_STAR_MR );
        }
    }
}

} // namespace trsv
} // namespace El
