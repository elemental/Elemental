/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename F>
inline void
TrsvUN( UnitOrNonUnit diag, const DistMatrix<F>& U, DistMatrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("internal::TrsvUN");
    if( U.Grid() != x.Grid() )
        throw std::logic_error("{U,x} must be distributed over the same grid");
    if( U.Height() != U.Width() )
        throw std::logic_error("U must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( U.Width() != xLength )
        throw std::logic_error("Nonconformal TrsvUN");
#endif
    const Grid& g = U.Grid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<F> U01(g),
                      U11(g);
        DistMatrix<F> 
            xT(g),  x0(g),
            xB(g),  x1(g),
                    x2(g);

        // Temporary distributions
        DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
        DistMatrix<F,STAR,STAR> x1_STAR_STAR(g);
        DistMatrix<F,MR,  STAR> x1_MR_STAR(g);
        DistMatrix<F,MC,  STAR> z_MC_STAR(g);

        // Views of z[MC,* ], which will store updates to x
        DistMatrix<F,MC,STAR> z0_MC_STAR(g),
                              z1_MC_STAR(g);

        z_MC_STAR.AlignWith( U );
        Zeros( x.Height(), 1, z_MC_STAR );

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

            const int n0 = x0.Height();
            const int n1 = x1.Height();
            U01.LockedView( U, 0,  n0, n0, n1 );
            U11.LockedView( U, n0, n0, n1, n1 );
            z0_MC_STAR.View( z_MC_STAR, 0,  0, n0, 1 );
            z1_MC_STAR.View( z_MC_STAR, n0, 0, n1, 1 );

            x1_MR_STAR.AlignWith( U01 );
            //----------------------------------------------------------------//
            if( x2.Height() != 0 )
                x1.SumScatterUpdate( F(1), z1_MC_STAR );

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, NORMAL, diag,
              U11_STAR_STAR.LockedLocalMatrix(),
              x1_STAR_STAR.LocalMatrix() );
            x1 = x1_STAR_STAR;

            x1_MR_STAR = x1_STAR_STAR;
            Gemv
            ( NORMAL, F(-1), 
              U01.LockedLocalMatrix(), 
              x1_MR_STAR.LockedLocalMatrix(),
              F(1), z0_MC_STAR.LocalMatrix() );
            //----------------------------------------------------------------//
            x1_MR_STAR.FreeAlignments();

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
        DistMatrix<F> U01(g),
                      U11(g);
        DistMatrix<F> 
            xL(g), xR(g),
            x0(g), x1(g), x2(g);

        // Temporary distributions
        DistMatrix<F,STAR,STAR> U11_STAR_STAR(g);
        DistMatrix<F,STAR,STAR> x1_STAR_STAR(g);
        DistMatrix<F,STAR,MR  > x1_STAR_MR(g);
        DistMatrix<F,MC,  MR  > z1(g);
        DistMatrix<F,MR,  MC  > z1_MR_MC(g);
        DistMatrix<F,STAR,MC  > z_STAR_MC(g);

        // Views of z[* ,MC]
        DistMatrix<F,STAR,MC>  z0_STAR_MC(g),
                               z1_STAR_MC(g);

        z_STAR_MC.AlignWith( U );
        Zeros( 1, x.Width(), z_STAR_MC );

        // Start the algorithm
        PartitionLeft( x,  xL, xR, 0 );
        while( xL.Width() > 0 )
        {
            RepartitionLeft
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );

            const int n0 = x0.Width();
            const int n1 = x1.Width();
            U01.LockedView( U, 0,  n0, n0, n1 );
            U11.LockedView( U, n0, n0, n1, n1 );
            z0_STAR_MC.View( z_STAR_MC, 0, 0,  1, n0 );
            z1_STAR_MC.View( z_STAR_MC, 0, n0, 1, n1 );

            x1_STAR_MR.AlignWith( U01 );
            z1.AlignWith( x1 );
            //----------------------------------------------------------------//
            if( x2.Width() != 0 )
            {
                z1_MR_MC.SumScatterFrom( z1_STAR_MC );
                z1 = z1_MR_MC;
                Axpy( F(1), z1, x1 );
            }

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, NORMAL, diag,
              U11_STAR_STAR.LockedLocalMatrix(),
              x1_STAR_STAR.LocalMatrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MR = x1_STAR_STAR;
            Gemv
            ( NORMAL, F(-1), 
              U01.LockedLocalMatrix(), 
              x1_STAR_MR.LockedLocalMatrix(),
              F(1), z0_STAR_MC.LocalMatrix() );
            //----------------------------------------------------------------//
            x1_STAR_MR.FreeAlignments();
            z1.FreeAlignments(); 

            SlidePartitionLeft
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
