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
UT
( Orientation orientation, UnitOrNonUnit diag, 
  const DistMatrix<F>& U, DistMatrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("trsv::UT");
        if( U.Grid() != x.Grid() )
            LogicError("{U,x} must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("Expected a (conjugate-)transpose option");
        if( U.Height() != U.Width() )
            LogicError("U must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( U.Width() != xLength )
            LogicError("Nonconformal");
    )
    const Int m = U.Height();
    const Int bsize = Blocksize();
    const Grid& g = U.Grid();

    // Matrix views 
    DistMatrix<F> U11(g), U12(g), x1(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), x1_STAR_STAR(g);

    if( x.Width() == 1 )
    {
        DistMatrix<F,MC,  STAR> x1_MC_STAR(g);
        DistMatrix<F,MR,  MC  > z1_MR_MC(g);
        DistMatrix<F,MR,  STAR> z_MR_STAR(g);

        // Views of z[MR,* ]
        DistMatrix<F,MR,STAR> z1_MR_STAR(g), z2_MR_STAR(g);

        z_MR_STAR.AlignWith( U );
        Zeros( z_MR_STAR, m, 1 );

        for( Int k=0; k<m; k+=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedViewRange( U11, U, k, k,    k+nb, k+nb );
            LockedViewRange( U12, U, k, k+nb, k+nb, m    );

            ViewRange( x1, x, k, 0, k+nb, 1 );

            ViewRange( z1_MR_STAR, z_MR_STAR, k,    0, k+nb, 1 );
            ViewRange( z2_MR_STAR, z_MR_STAR, k+nb, 0, m,    1 );

            if( k != 0 )
            {
                z1_MR_MC.RowSumScatterFrom( z1_MR_STAR );
                Axpy( F(1), z1_MR_MC, x1 );
            }

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, orientation, diag,
              U11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_MC_STAR.AlignWith( U12 );
            x1_MC_STAR = x1_STAR_STAR;
            LocalGemv( orientation, F(-1), U12, x1_MC_STAR, F(1), z2_MR_STAR );
        }
    }
    else
    {
        DistMatrix<F,STAR,MC> x1_STAR_MC(g);
        DistMatrix<F,STAR,MR> z_STAR_MR(g);

        // Views of z[* ,MR]
        DistMatrix<F,STAR,MR> z1_STAR_MR(g), z2_STAR_MR(g);

        z_STAR_MR.AlignWith( U );
        Zeros( z_STAR_MR, 1, m );

        for( Int k=0; k<m; k+=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedViewRange( U11, U, k, k,    k+nb, k+nb );
            LockedViewRange( U12, U, k, k+nb, k+nb, m    );

            ViewRange( x1, x, 0, k, 1, k+nb );

            ViewRange( z1_STAR_MR, z_STAR_MR, 0, k,    1, k+nb );
            ViewRange( z2_STAR_MR, z_STAR_MR, 0, k+nb, 1, m    );

            if( k != 0 )
                x1.ColSumScatterUpdate( F(1), z1_STAR_MR );

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, orientation, diag,
              U11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MC.AlignWith( U12 );
            x1_STAR_MC = x1_STAR_STAR;
            LocalGemv( orientation, F(-1), U12, x1_STAR_MC, F(1), z2_STAR_MR );
        }
    }
}

} // namespace trsv
} // namespace El
