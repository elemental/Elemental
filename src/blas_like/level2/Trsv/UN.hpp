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
void UN
( UnitOrNonUnit diag, 
  const AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<F>& xPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( UPre, xPre );
      if( UPre.Height() != UPre.Width() )
          LogicError("U must be square");
      if( xPre.Width() != 1 && xPre.Height() != 1 )
          LogicError("x must be a vector");
      const Int xLength =
          ( xPre.Width() == 1 ? xPre.Height() : xPre.Width() );
      if( UPre.Width() != xLength )
          LogicError("Nonconformal");
    )
    const Int m = UPre.Height();
    const Int bsize = Blocksize();
    const Int kLast = LastOffset( m, bsize );
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> xProx( xPre );
    auto& U = UProx.GetLocked();
    auto& x = xProx.Get();

    // Matrix views
    DistMatrix<F> U01(g), U11(g), x1(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g), x1_STAR_STAR(g);

    if( x.Width() == 1 )
    {
        DistMatrix<F,MR,STAR> x1_MR_STAR(g);
        DistMatrix<F,MC,STAR> z_MC_STAR(g);

        // Views of z[MC,* ], which will store updates to x
        DistMatrix<F,MC,STAR> z0_MC_STAR(g), z1_MC_STAR(g);

        z_MC_STAR.AlignWith( U );
        z_MC_STAR.Resize( m, 1 );
        Zero( z_MC_STAR );

        for( Int k=kLast; k>=0; k-=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedView( U01, U, IR(0,k), IR(k,k+nb) );
            LockedView( U11, U, IR(k,k+nb), IR(k,k+nb) );

            View( x1, x, IR(k,k+nb), ALL );

            View( z0_MC_STAR, z_MC_STAR, IR(0,k), ALL );
            View( z1_MC_STAR, z_MC_STAR, IR(k,k+nb), ALL );

            if( k+nb != m )
                AxpyContract( F(1), z1_MC_STAR, x1 );

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, NORMAL, diag,
              U11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_MR_STAR.AlignWith( U01 );
            x1_MR_STAR = x1_STAR_STAR;
            LocalGemv( NORMAL, F(-1), U01, x1_MR_STAR, F(1), z0_MC_STAR );
        }
    }
    else
    {
        DistMatrix<F,STAR,MR> x1_STAR_MR(g);
        DistMatrix<F,MR,  MC> z1_MR_MC(g);
        DistMatrix<F,STAR,MC> z_STAR_MC(g);

        // Views of z[* ,MC]
        DistMatrix<F,STAR,MC> z0_STAR_MC(g), z1_STAR_MC(g);

        z_STAR_MC.AlignWith( U );
        z_STAR_MC.Resize( 1, m );
        Zero( z_STAR_MC );

        for( Int k=kLast; k>=0; k-=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedView( U01, U, IR(0,k), IR(k,k+nb) );
            LockedView( U11, U, IR(k,k+nb), IR(k,k+nb) );

            View( x1, x, ALL, IR(k,k+nb) );

            View( z0_STAR_MC, z_STAR_MC, ALL, IR(0,k) );
            View( z1_STAR_MC, z_STAR_MC, ALL, IR(k,k+nb) );

            if( k+nb != m )
            {
                Contract( z1_STAR_MC, z1_MR_MC );
                x1 += z1_MR_MC;
            }

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, NORMAL, diag,
              U11_STAR_STAR.LockedMatrix(), x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MR.AlignWith( U01 );
            x1_STAR_MR = x1_STAR_STAR;
            LocalGemv( NORMAL, F(-1), U01, x1_STAR_MR, F(1), z0_STAR_MC );
        }
    }
}

} // namespace trsv
} // namespace El
