/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_TRSV_LT_HPP
#define ELEM_TRSV_LT_HPP

#include ELEM_AXPY_INC
#include ELEM_ZEROS_INC
#include ELEM_GEMV_INC

namespace elem {
namespace internal {

template<typename F>
inline void
TrsvLT
( Orientation orientation, UnitOrNonUnit diag, 
  const DistMatrix<F>& L, DistMatrix<F>& x )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrsvLT");
        if( L.Grid() != x.Grid() )
            LogicError("{L,x} must be distributed over the same grid");
        if( orientation == NORMAL )
            LogicError("TrsvLT expects a (conjugate-)transpose option");
        if( L.Height() != L.Width() )
            LogicError("L must be square");
        if( x.Width() != 1 && x.Height() != 1 )
            LogicError("x must be a vector");
        const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
        if( L.Width() != xLength )
            LogicError("Nonconformal TrsvLT");
    )
    const Grid& g = L.Grid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<F> L10(g), L11(g);
        DistMatrix<F> 
            xT(g),  x0(g),
            xB(g),  x1(g),
                    x2(g);

        // Temporary distributions
        DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
        DistMatrix<F,STAR,STAR> x1_STAR_STAR(g);
        DistMatrix<F,MC,  STAR> x1_MC_STAR(g);
        DistMatrix<F,MC,  MR  > z1(g);
        DistMatrix<F,MR,  MC  > z1_MR_MC(g);
        DistMatrix<F,MR,  STAR> z_MR_STAR(g);

        // Views of z[MR,* ]
        DistMatrix<F,MR,STAR> z0_MR_STAR(g),
                              z1_MR_STAR(g);

        z_MR_STAR.AlignWith( L );
        Zeros( z_MR_STAR, x.Height(), 1 );

        // Start the algorithm
        PartitionUp
        ( x, xT,
             xB, 0 );
        while( xT.Height() > 0 )
        {
            RepartitionUp
            ( xT,  x0,
                   x1,
             /**/ /**/
              xB,  x2 );

            const Int n0 = x0.Height();
            const Int n1 = x1.Height();
            LockedView( L10, L, n0, 0,  n1, n0 );
            LockedView( L11, L, n0, n0, n1, n1 );
            View( z0_MR_STAR, z_MR_STAR, 0,  0, n0, 1 );
            View( z1_MR_STAR, z_MR_STAR, n0, 0, n1, 1 );

            x1_MC_STAR.AlignWith( L10 );
            z1.AlignWith( x1 );
            //----------------------------------------------------------------//
            if( x2.Height() != 0 )
            {
                z1_MR_MC.RowSumScatterFrom( z1_MR_STAR );
                z1 = z1_MR_MC;
                Axpy( F(1), z1, x1 );
            }

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            Trsv
            ( LOWER, orientation, diag,
              L11_STAR_STAR.LockedMatrix(),
              x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_MC_STAR = x1_STAR_STAR;
            LocalGemv( orientation, F(-1), L10, x1_MC_STAR, F(1), z0_MR_STAR );
            //----------------------------------------------------------------//

            SlidePartitionUp
            ( xT,  x0,
             /**/ /**/
                   x1,
              xB,  x2 );
        }
    }
    else
    {
        // Matrix views 
        DistMatrix<F> L10(g), L11(g);
        DistMatrix<F> 
            xL(g), xR(g),
            x0(g), x1(g), x2(g);

        // Temporary distributions
        DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
        DistMatrix<F,STAR,STAR> x1_STAR_STAR(g);
        DistMatrix<F,STAR,MC  > x1_STAR_MC(g);
        DistMatrix<F,STAR,MR  > z_STAR_MR(g);

        // Views of z[* ,MR], which will store updates to x
        DistMatrix<F,STAR,MR> z0_STAR_MR(g),
                              z1_STAR_MR(g);

        z_STAR_MR.AlignWith( L );
        Zeros( z_STAR_MR, 1, x.Width() );

        // Start the algorithm
        PartitionLeft( x,  xL, xR, 0 );
        while( xL.Width() > 0 )
        {
            RepartitionLeft
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );

            const Int n0 = x0.Width();
            const Int n1 = x1.Width();
            LockedView( L10, L, n0, 0,  n1, n0 );
            LockedView( L11, L, n0, n0, n1, n1 );
            View( z0_STAR_MR, z_STAR_MR, 0, 0,  1, n0 );
            View( z1_STAR_MR, z_STAR_MR, 0, n0, 1, n1 );

            x1_STAR_MC.AlignWith( L10 );
            //----------------------------------------------------------------//
            if( x2.Width() != 0 )
                x1.ColSumScatterUpdate( F(1), z1_STAR_MR );

            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            Trsv
            ( LOWER, orientation, diag,
              L11_STAR_STAR.LockedMatrix(),
              x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MC = x1_STAR_STAR;
            LocalGemv( orientation, F(-1), L10, x1_STAR_MC, F(1), z0_STAR_MR );
            //----------------------------------------------------------------//

            SlidePartitionLeft
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );
        }
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRSV_LT_HPP
