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
void LT
( Orientation orientation,
  UnitOrNonUnit diag, 
  const AbstractDistMatrix<F>& LPre,
        AbstractDistMatrix<F>& xPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( orientation == NORMAL )
          LogicError("Expected a (conjugate-)transpose option");
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
    const Int kLast = LastOffset( m, bsize );
    const Grid& g = LPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> LProx( LPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> xProx( xPre );
    auto& L = LProx.GetLocked();
    auto& x = xProx.Get();

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
        z_MR_STAR.Resize( m, 1 );
        Zero( z_MR_STAR );

        for( Int k=kLast; k>=0; k-=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedView( L10, L, IR(k,k+nb), IR(0,k) );
            LockedView( L11, L, IR(k,k+nb), IR(k,k+nb) );

            View( x1, x, IR(k,k+nb), ALL );

            View( z0_MR_STAR, z_MR_STAR, IR(0,k), ALL );
            View( z1_MR_STAR, z_MR_STAR, IR(k,k+nb), ALL );

            if( k+nb != m )
            {
                Contract( z1_MR_STAR, z1_MR_MC );
                x1 += z1_MR_MC;
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
        z_STAR_MR.Resize( 1, m );
        Zero( z_STAR_MR );

        for( Int k=kLast; k>=0; k-=bsize )
        {
            const Int nb = Min(bsize,m-k);

            LockedView( L10, L, IR(k,k+nb), IR(0,k) );
            LockedView( L11, L, IR(k,k+nb), IR(k,k+nb) );

            View( x1, x, ALL, IR(k,k+nb) );

            View( z0_STAR_MR, z_STAR_MR, ALL, IR(0,k) );
            View( z1_STAR_MR, z_STAR_MR, ALL, IR(k,k+nb) );

            if( k+nb != m )
                AxpyContract( F(1), z1_STAR_MR, x1 );

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
