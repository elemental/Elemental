/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_BLAS_TRSV_UT_HPP
#define ELEM_BLAS_TRSV_UT_HPP

#include "elemental/blas-like/level1/Axpy.hpp"
#include "elemental/blas-like/level2/Gemv.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace internal {

template<typename F>
inline void
TrsvUT
( Orientation orientation, UnitOrNonUnit diag, 
  const DistMatrix<F>& U, DistMatrix<F>& x )
{
#ifndef RELEASE
    CallStackEntry cse("internal::TrsvUT");
    if( U.Grid() != x.Grid() )
        LogicError("{U,x} must be distributed over the same grid");
    if( orientation == NORMAL )
        LogicError("TrsvUT expects a (conjugate-)transpose option");
    if( U.Height() != U.Width() )
        LogicError("U must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        LogicError("x must be a vector");
    const Int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( U.Width() != xLength )
        LogicError("Nonconformal TrsvUT");
#endif
    const Grid& g = U.Grid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<F> U11(g), U12(g);
        DistMatrix<F> 
            xT(g),  x0(g),
            xB(g),  x1(g),
                    x2(g);

        // Temporary distributions
        DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
        DistMatrix<F,STAR,STAR> x1_STAR_STAR(g);
        DistMatrix<F,MC,  STAR> x1_MC_STAR(g);
        DistMatrix<F,MC,  MR  > z1(g);
        DistMatrix<F,MR,  MC  > z1_MR_MC(g);
        DistMatrix<F,MR,  STAR> z_MR_STAR(g);

        // Views of z[MR,* ]
        DistMatrix<F,MR,STAR> z1_MR_STAR(g),
                              z2_MR_STAR(g);

        z_MR_STAR.AlignWith( U );
        Zeros( z_MR_STAR, x.Height(), 1 );

        // Start the algorithm
        PartitionDown
        ( x, xT,
             xB, 0 );
        while( xB.Height() > 0 )
        {
            RepartitionDown
            ( xT,  x0,
             /**/ /**/
                   x1,
              xB,  x2 );

            const Int n0 = x0.Height();
            const Int n1 = x1.Height();
            const Int n2 = x2.Height();
            LockedView( U11, U, n0, n0,    n1, n1 );
            LockedView( U12, U, n0, n0+n1, n1, n2 );
            View( z1_MR_STAR, z_MR_STAR, n0,    0, n1, 1 );
            View( z2_MR_STAR, z_MR_STAR, n0+n1, 0, n2, 1 );

            x1_MC_STAR.AlignWith( U12 );
            z1.AlignWith( x1 );
            //----------------------------------------------------------------//
            if( x0.Height() != 0 )
            {
                z1_MR_MC.SumScatterFrom( z1_MR_STAR );
                z1 = z1_MR_MC;
                Axpy( F(1), z1, x1 );
            }

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, orientation, diag,
              U11_STAR_STAR.LockedMatrix(),
              x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_MC_STAR = x1_STAR_STAR;
            LocalGemv( orientation, F(-1), U12, x1_MC_STAR, F(1), z2_MR_STAR );
            //----------------------------------------------------------------//

            SlidePartitionDown
            ( xT,  x0,
                   x1,
             /**/ /**/
              xB,  x2 );
        }
    }
    else
    {
        // Matrix views 
        DistMatrix<F> U11(g), U12(g);
        DistMatrix<F> 
            xL(g), xR(g),
            x0(g), x1(g), x2(g);

        // Temporary distributions
        DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
        DistMatrix<F,STAR,STAR> x1_STAR_STAR(g);
        DistMatrix<F,STAR,MC  > x1_STAR_MC(g);
        DistMatrix<F,STAR,MR  > z_STAR_MR(g);

        // Views of z[* ,MR]
        DistMatrix<F,STAR,MR> z1_STAR_MR(g),
                              z2_STAR_MR(g);

        z_STAR_MR.AlignWith( U );
        Zeros( z_STAR_MR, 1, x.Width() );

        // Start the algorithm
        PartitionRight( x,  xL, xR, 0 );
        while( xR.Width() > 0 )
        {
            RepartitionRight
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );

            const Int n0 = x0.Width();
            const Int n1 = x1.Width();
            const Int n2 = x2.Width();
            LockedView( U11, U, n0, n0,    n1, n1 );
            LockedView( U12, U, n0, n0+n1, n1, n2 );
            View( z1_STAR_MR, z_STAR_MR, 0, n0,    1, n1 );
            View( z2_STAR_MR, z_STAR_MR, 0, n0+n1, 1, n2 );

            x1_STAR_MC.AlignWith( U12 );
            //----------------------------------------------------------------//
            if( x2.Width() != 0 )
                x1.SumScatterUpdate( F(1), z1_STAR_MR );

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, orientation, diag,
              U11_STAR_STAR.LockedMatrix(),
              x1_STAR_STAR.Matrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MC = x1_STAR_STAR;
            LocalGemv( orientation, F(-1), U12, x1_STAR_MC, F(1), z2_STAR_MR );
            //----------------------------------------------------------------//

            SlidePartitionRight
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );
        }
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_BLAS_TRSV_UT_HPP
