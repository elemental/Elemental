/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SYMV_L_HPP
#define ELEM_SYMV_L_HPP

#include ELEM_MAKETRIANGULAR_INC
#include ELEM_SETDIAGONAL_INC
#include ELEM_GEMV_INC

namespace elem {
namespace internal {

template<typename T>
inline void
LocalSymvColAccumulateL
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::LocalSymvColAccumulateL");
        if( A.Grid() != x_MC_STAR.Grid() ||
            x_MC_STAR.Grid() != x_MR_STAR.Grid() ||
            x_MR_STAR.Grid() != z_MC_STAR.Grid() ||
            z_MC_STAR.Grid() != z_MR_STAR.Grid() )
            LogicError("{A,x,z} must be distributed over the same grid");
        if( x_MC_STAR.Width() != 1 || x_MR_STAR.Width() != 1 ||
            z_MC_STAR.Width() != 1 || z_MR_STAR.Width() != 1 )
            LogicError("Expected x and z to be column vectors");
        if( A.Height() != A.Width() || 
            A.Height() != x_MC_STAR.Height() ||
            A.Height() != x_MR_STAR.Height() ||
            A.Height() != z_MC_STAR.Height() ||
            A.Height() != z_MR_STAR.Height() )
            LogicError
            ("Nonconformal LocalSymvColAccumulateL: \n",
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

    // Matrix views
    DistMatrix<T> A11(g),
                  A21(g);
    DistMatrix<T> D11(g);

    DistMatrix<T,MC,STAR> 
        xT_MC_STAR(g),  x0_MC_STAR(g),
        xB_MC_STAR(g),  x1_MC_STAR(g),
                        x2_MC_STAR(g);
    DistMatrix<T,MR,STAR>  x1_MR_STAR(g);
    DistMatrix<T,MC,STAR> z1_MC_STAR(g),
                          z2_MC_STAR(g);
    DistMatrix<T,MR,STAR> z1_MR_STAR(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const Int ratio = Max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*LocalSymvBlocksize<T>() );
    LockedPartitionDown
    ( x_MC_STAR, xT_MC_STAR,
                 xB_MC_STAR, 0 );
    while( xT_MC_STAR.Height() < x_MC_STAR.Height() )
    {
        LockedRepartitionDown
        ( xT_MC_STAR,  x0_MC_STAR,
         /**********/ /**********/
                       x1_MC_STAR,
          xB_MC_STAR,  x2_MC_STAR );

        const Int n0 = x0_MC_STAR.Height();
        const Int n1 = x1_MC_STAR.Height();
        const Int n2 = x2_MC_STAR.Height();
        LockedView( A11, A, n0,    n0, n1, n1 );
        LockedView( A21, A, n0+n1, n0, n2, n1 );
        LockedView( x1_MR_STAR, x_MR_STAR, n0, 0, n1, 1 );
        View( z1_MC_STAR, z_MC_STAR, n0,    0, n1, 1 );
        View( z2_MC_STAR, z_MC_STAR, n0+n1, 0, n2, 1 );
        View( z1_MR_STAR, z_MR_STAR, n0,    0, n1, 1 );
 
        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        // TODO: These diagonal block updates can be greatly improved
        D11 = A11;
        MakeTriangular( LOWER, D11 );
        LocalGemv( NORMAL, alpha, D11, x1_MR_STAR, T(1), z1_MC_STAR );
        SetDiagonal( D11, T(0) );
        LocalGemv( orientation, alpha, D11, x1_MC_STAR, T(1), z1_MR_STAR );

        LocalGemv( NORMAL, alpha, A21, x1_MR_STAR, T(1), z2_MC_STAR );
        LocalGemv( orientation, alpha, A21, x2_MC_STAR, T(1), z1_MR_STAR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( xT_MC_STAR,  x0_MC_STAR,
                       x1_MC_STAR,
         /**********/ /**********/
          xB_MC_STAR,  x2_MC_STAR );
    }
    PopBlocksizeStack();
}

template<typename T>
inline void
LocalSymvRowAccumulateL
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::LocalSymvRowAccumulateL");
        if( A.Grid() != x_STAR_MC.Grid() ||
            x_STAR_MC.Grid() != x_STAR_MR.Grid() ||
            x_STAR_MR.Grid() != z_STAR_MC.Grid() ||
            z_STAR_MC.Grid() != z_STAR_MR.Grid()   )
            LogicError("{A,x,z} must be distributed over the same grid");
        if( x_STAR_MC.Height() != 1 || x_STAR_MR.Height() != 1 ||
            z_STAR_MC.Height() != 1 || z_STAR_MR.Height() != 1    )
            LogicError("Expected x and z to be row vectors");
        if( A.Height() != A.Width() || 
            A.Height() != x_STAR_MC.Width() ||
            A.Height() != x_STAR_MR.Width() ||
            A.Height() != z_STAR_MC.Width() ||
            A.Height() != z_STAR_MR.Width()   )
            LogicError
            ("Nonconformal LocalSymvRowAccumulateL: \n"
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

    // Matrix views
    DistMatrix<T> A11(g),
                  A21(g);
    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MC> 
        xL_STAR_MC(g), xR_STAR_MC(g),
        x0_STAR_MC(g), x1_STAR_MC(g), x2_STAR_MC(g);
    DistMatrix<T,STAR,MR> x1_STAR_MR(g);
    DistMatrix<T,STAR,MC> z1_STAR_MC(g), z2_STAR_MC(g);
    DistMatrix<T,STAR,MR> z1_STAR_MR(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const Int ratio = Max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*LocalSymvBlocksize<T>() );
                 
    LockedPartitionRight( x_STAR_MC,  xL_STAR_MC, xR_STAR_MC, 0 );
    while( xL_STAR_MC.Width() < x_STAR_MC.Width() )
    {
        LockedRepartitionRight
        ( xL_STAR_MC, /**/ xR_STAR_MC, 
          x0_STAR_MC, /**/ x1_STAR_MC, x2_STAR_MC );

        const Int n0 = x0_STAR_MC.Width();
        const Int n1 = x1_STAR_MC.Width();
        const Int n2 = x2_STAR_MC.Width();
        LockedView( A11, A, n0,    n0, n1, n1 );
        LockedView( A21, A, n0+n1, n0, n2, n1 );
        LockedView( x1_STAR_MR, x_STAR_MR, 0, n0, 1, n1 );
        View( z1_STAR_MC, z_STAR_MC, 0, n0,    1, n1 );
        View( z2_STAR_MC, z_STAR_MC, 0, n0+n1, 1, n2 );
        View( z1_STAR_MR, z_STAR_MR, 0, n0,    1, n1 );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        // TODO: These diagonal block updates can be greatly improved
        D11 = A11;
        MakeTriangular( LOWER, D11 );
        LocalGemv( NORMAL, alpha, D11, x1_STAR_MR, T(1), z1_STAR_MC );
        SetDiagonal( D11, T(0) );
        LocalGemv( orientation, alpha, D11, x1_STAR_MC, T(1), z1_STAR_MR );

        LocalGemv( NORMAL, alpha, A21, x1_STAR_MR, T(1), z2_STAR_MC );
        LocalGemv( orientation, alpha, A21, x2_STAR_MC, T(1), z1_STAR_MR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( xL_STAR_MC,             /**/ xR_STAR_MC,
          x0_STAR_MC, x1_STAR_MC, /**/ x2_STAR_MC );
    }
    PopBlocksizeStack();
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_SYMV_L_HPP
