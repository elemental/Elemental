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
UT
( Orientation orientation,
  UnitOrNonUnit diag, 
  const AbstractDistMatrix<F>& UPre,
        AbstractDistMatrix<F>& xPre )
{
    DEBUG_ONLY(
      CSE cse("trsv::UT");
      if( orientation == NORMAL )
          LogicError("Expected a (conjugate-)transpose option");
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
    const Grid& g = UPre.Grid();

    DistMatrixReadProxy<F,F,MC,MR> UProx( UPre );
    DistMatrixReadWriteProxy<F,F,MC,MR> xProx( xPre );
    auto& U = UProx.GetLocked();
    auto& x = xProx.Get();

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

            LockedView( U11, U, IR(k,k+nb), IR(k,k+nb) );
            LockedView( U12, U, IR(k,k+nb), IR(k+nb,m) );

            View( x1, x, IR(k,k+nb), ALL );

            View( z1_MR_STAR, z_MR_STAR, IR(k,k+nb), ALL );
            View( z2_MR_STAR, z_MR_STAR, IR(k+nb,m), ALL );

            if( k != 0 )
            {
                Contract( z1_MR_STAR, z1_MR_MC );
                x1 += z1_MR_MC;
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

            LockedView( U11, U, IR(k,k+nb), IR(k,k+nb) );
            LockedView( U12, U, IR(k,k+nb), IR(k+nb,m) );

            View( x1, x, ALL, IR(k,k+nb) );

            View( z1_STAR_MR, z_STAR_MR, ALL, IR(k,k+nb) );
            View( z2_STAR_MR, z_STAR_MR, ALL, IR(k+nb,m) );

            if( k != 0 )
                AxpyContract( F(1), z1_STAR_MR, x1 );

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
