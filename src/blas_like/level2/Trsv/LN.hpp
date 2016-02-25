/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace trsv {

template<typename F>
inline void
LN
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& xPre )
{
    DEBUG_ONLY(
      CSE cse("trsv::LN");
      AssertSameGrids( LPre, xPre );
      if( LPre.Height() != LPre.Width() )
          LogicError("L must be square");
      if( xPre.Width() != 1 && xPre.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength = 
          ( xPre.Width() == 1 ? xPre.Height() : xPre.Width() );
      if( LPre.Width() != xLength )
          LogicError("Nonconformal");
    )
    const Int m = LPre.Height();
    const Int bsize = Blocksize();
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> xProx( xPre );
    auto& L = LProx.GetLocked();
    auto& x = xProx.Get();

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

            LockedView( L11, L, IR(k,k+nb), IR(k,k+nb) );
            LockedView( L21, L, IR(k+nb,m), IR(k,k+nb) );

            View( x1, x, IR(k,k+nb), ALL );

            View( z1_MC_STAR, z_MC_STAR, IR(k,k+nb), ALL );
            View( z2_MC_STAR, z_MC_STAR, IR(k+nb,m), ALL );

            if( k != 0 )
                AxpyContract( F(1), z1_MC_STAR, x1 );

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
        DistMatrix<F,MR,  MC> z1_MR_MC(g);
        DistMatrix<F,STAR,MC> z_STAR_MC(g);

        // Views of z[* ,MC]
        DistMatrix<F,STAR,MC> z1_STAR_MC(g), z2_STAR_MC(g);

        z_STAR_MC.AlignWith( L );
        Zeros( z_STAR_MC, 1, m );

        for( Int k=0; k<m; k+=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedView( L11, L, IR(k,k+nb), IR(k,k+nb) );
            LockedView( L21, L, IR(k+nb,m), IR(k,k+nb) );

            View( x1, x, ALL, IR(k,k+nb) );

            View( z1_STAR_MC, z_STAR_MC, ALL, IR(k,k+nb) );
            View( z2_STAR_MC, z_STAR_MC, ALL, IR(k+nb,m) );

            if( k != 0 )
            {
                Contract( z1_STAR_MC, z1_MR_MC );
                x1 += z1_MR_MC;
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

} // namespace trsv
} // namespace El
