/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace symv {

template<typename T>
inline void LocalColAccumulateL
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("symv::LocalColAccumulateL");
        AssertSameGrids( A, x_MC_STAR, x_MR_STAR, z_MC_STAR, z_MR_STAR );
        if( x_MC_STAR.Width() != 1 || x_MR_STAR.Width() != 1 ||
            z_MC_STAR.Width() != 1 || z_MR_STAR.Width() != 1 )
            LogicError("Expected x and z to be column vectors");
        if( A.Height() != A.Width() || 
            A.Height() != x_MC_STAR.Height() ||
            A.Height() != x_MR_STAR.Height() ||
            A.Height() != z_MC_STAR.Height() ||
            A.Height() != z_MR_STAR.Height() )
            LogicError
            ("Nonconformal: \n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  x[MC,* ] ~ ",x_MC_STAR.Height()," x ",x_MC_STAR.Width(),"\n", 
             "  x[MR,* ] ~ ",x_MR_STAR.Height()," x ",x_MR_STAR.Width(),"\n", 
             "  z[MC,* ] ~ ",z_MC_STAR.Height()," x ",z_MC_STAR.Width(),"\n", 
             "  z[MR,* ] ~ ",z_MR_STAR.Height()," x ",z_MR_STAR.Width(),"\n"); 
        if( x_MC_STAR.ColAlign() != A.ColAlign() ||
            x_MR_STAR.ColAlign() != A.RowAlign() ||
            z_MC_STAR.ColAlign() != A.ColAlign() ||
            z_MR_STAR.ColAlign() != A.RowAlign() )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<T> D11(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const Int bsize = Max(g.Height(),g.Width())*LocalSymvBlocksize<T>();
    const Int n = A.Height();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const Range<Int> ind1( k, k+nb ), ind2( k+nb, n );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto x1_MC_STAR = x_MC_STAR( ind1, IR(0,1) );
        auto x2_MC_STAR = x_MC_STAR( ind2, IR(0,1) );
        auto x1_MR_STAR = x_MR_STAR( ind1, IR(0,1) );
        auto z1_MC_STAR = z_MC_STAR( ind1, IR(0,1) );
        auto z2_MC_STAR = z_MC_STAR( ind2, IR(0,1) );
        auto z1_MR_STAR = z_MR_STAR( ind1, IR(0,1) );
 
        D11.AlignWith( A11 );
        // TODO: These diagonal block updates can be greatly improved
        D11 = A11;
        MakeTrapezoidal( LOWER, D11 );
        LocalGemv( NORMAL,      alpha, D11, x1_MR_STAR, T(1), z1_MC_STAR );
        SetDiagonal( D11, T(0) );
        LocalGemv( orientation, alpha, D11, x1_MC_STAR, T(1), z1_MR_STAR );

        LocalGemv( NORMAL,      alpha, A21, x1_MR_STAR, T(1), z2_MC_STAR );
        LocalGemv( orientation, alpha, A21, x2_MC_STAR, T(1), z1_MR_STAR );
    }
}

template<typename T>
inline void LocalRowAccumulateL
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("symv::LocalRowAccumulateL");
        AssertSameGrids( A, x_STAR_MC, x_STAR_MR, z_STAR_MC, z_STAR_MR );
        if( x_STAR_MC.Height() != 1 || x_STAR_MR.Height() != 1 ||
            z_STAR_MC.Height() != 1 || z_STAR_MR.Height() != 1    )
            LogicError("Expected x and z to be row vectors");
        if( A.Height() != A.Width() || 
            A.Height() != x_STAR_MC.Width() ||
            A.Height() != x_STAR_MR.Width() ||
            A.Height() != z_STAR_MC.Width() ||
            A.Height() != z_STAR_MR.Width()   )
            LogicError
            ("Nonconformal: \n"
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  x[* ,MC] ~ ",x_STAR_MC.Height()," x ",x_STAR_MC.Width(),"\n", 
             "  x[* ,MR] ~ ",x_STAR_MR.Height()," x ",x_STAR_MR.Width(),"\n", 
             "  z[* ,MC] ~ ",z_STAR_MC.Height()," x ",z_STAR_MC.Width(),"\n", 
             "  z[* ,MR] ~ ",z_STAR_MR.Height()," x ",z_STAR_MR.Width(),"\n"); 
        if( x_STAR_MC.RowAlign() != A.ColAlign() ||
            x_STAR_MR.RowAlign() != A.RowAlign() ||
            z_STAR_MC.RowAlign() != A.ColAlign() ||
            z_STAR_MR.RowAlign() != A.RowAlign()   )
            LogicError("Partial matrix distributions are misaligned");
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrix<T> D11(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const Int bsize = Max(g.Height(),g.Width())*LocalSymvBlocksize<T>();
    const Int n = A.Height();
    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);
        const Range<Int> ind1( k, k+nb ), ind2( k+nb, n );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );
        auto x1_STAR_MC = x_STAR_MC( IR(0,1), ind1 );
        auto x2_STAR_MC = x_STAR_MC( IR(0,1), ind2 );
        auto x1_STAR_MR = x_STAR_MR( IR(0,1), ind1 );
        auto z1_STAR_MC = z_STAR_MC( IR(0,1), ind1 );
        auto z2_STAR_MC = z_STAR_MC( IR(0,1), ind2 );
        auto z1_STAR_MR = z_STAR_MR( IR(0,1), ind1 );

        D11.AlignWith( A11 );
        // TODO: These diagonal block updates can be greatly improved
        D11 = A11;
        MakeTrapezoidal( LOWER, D11 );
        LocalGemv( NORMAL,      alpha, D11, x1_STAR_MR, T(1), z1_STAR_MC );
        SetDiagonal( D11, T(0) );
        LocalGemv( orientation, alpha, D11, x1_STAR_MC, T(1), z1_STAR_MR );

        LocalGemv( NORMAL,      alpha, A21, x1_STAR_MR, T(1), z2_STAR_MC );
        LocalGemv( orientation, alpha, A21, x2_STAR_MC, T(1), z1_STAR_MR );
    }
}

} // namespace symv
} // namespace El
