/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename T>
inline void
LocalSymvColAccumulateU
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalSymvColAccumulateU");
    if( A.Grid() != x_MC_STAR.Grid() ||
        x_MC_STAR.Grid() != x_MR_STAR.Grid() ||
        x_MR_STAR.Grid() != z_MC_STAR.Grid() ||
        z_MC_STAR.Grid() != z_MR_STAR.Grid() )
        throw std::logic_error
        ("{A,x,z} must be distributed over the same grid");
    if( x_MC_STAR.Width() != 1 || x_MR_STAR.Width() != 1 ||
        z_MC_STAR.Width() != 1 || z_MR_STAR.Width() != 1 )
        throw std::logic_error("Expected x and z to be column vectors");
    if( A.Height() != A.Width() || 
        A.Height() != x_MC_STAR.Height() ||
        A.Height() != x_MR_STAR.Height() ||
        A.Height() != z_MC_STAR.Height() ||
        A.Height() != z_MR_STAR.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalSymvColAccumulateU: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x[MC,* ] ~ " << x_MC_STAR.Height() << " x " 
                               << x_MC_STAR.Width() << "\n"
            << "  x[MR,* ] ~ " << x_MR_STAR.Height() << " x " 
                               << x_MR_STAR.Width() << "\n"
            << "  z[MC,* ] ~ " << z_MC_STAR.Height() << " x " 
                               << z_MC_STAR.Width() << "\n"
            << "  z[MR,* ] ~ " << z_MR_STAR.Height() << " x " 
                               << z_MR_STAR.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
    if( x_MC_STAR.ColAlignment() != A.ColAlignment() ||
        x_MR_STAR.ColAlignment() != A.RowAlignment() ||
        z_MC_STAR.ColAlignment() != A.ColAlignment() ||
        z_MR_STAR.ColAlignment() != A.RowAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> A11(g), A12(g);
    DistMatrix<T> D11(g);

    DistMatrix<T,MC,STAR> x1_MC_STAR(g);
    DistMatrix<T,MR,STAR> 
        xT_MR_STAR(g),  x0_MR_STAR(g),
        xB_MR_STAR(g),  x1_MR_STAR(g),
                        x2_MR_STAR(g);
    DistMatrix<T,MC,STAR> z1_MC_STAR(g);
    DistMatrix<T,MR,STAR> z1_MR_STAR(g),
                          z2_MR_STAR(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*LocalSymvBlocksize<T>() );
                 
    LockedPartitionDown
    ( x_MR_STAR, xT_MR_STAR,
                 xB_MR_STAR, 0 );
    while( xT_MR_STAR.Height() < x_MR_STAR.Height() )
    {
        LockedRepartitionDown
        ( xT_MR_STAR,  x0_MR_STAR,
         /**********/ /**********/
                       x1_MR_STAR,
          xB_MR_STAR,  x2_MR_STAR );

        const int n0 = x0_MR_STAR.Height();
        const int n1 = x1_MR_STAR.Height();
        const int n2 = x2_MR_STAR.Height();
        LockedView( A11, A, n0, n0,    n1, n1 );
        LockedView( A12, A, n0, n0+n1, n1, n2 );
        LockedView( x1_MC_STAR, x_MC_STAR, n0, 0, n1, 1 );
        View( z1_MC_STAR, z_MC_STAR, n0,    0, n1, 1 );
        View( z1_MR_STAR, z_MR_STAR, n0,    0, n1, 1 );
        View( z2_MR_STAR, z_MR_STAR, n0+n1, 0, n2, 1 );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        // TODO: These diagonal block updates can be greatly improved
        D11 = A11;
        MakeTrapezoidal( LEFT, UPPER, 0, D11 );
        Gemv
        ( NORMAL, 
          alpha, D11.LockedLocalMatrix(), 
                 x1_MR_STAR.LockedLocalMatrix(),
          T(1),  z1_MC_STAR.LocalMatrix() );
        MakeTrapezoidal( LEFT, UPPER, 1, D11 );
        Gemv
        ( TRANSPOSE,
          alpha, D11.LockedLocalMatrix(),
                 x1_MC_STAR.LockedLocalMatrix(),
          T(1),  z1_MR_STAR.LocalMatrix() );

        Gemv
        ( NORMAL,
          alpha, A12.LockedLocalMatrix(),
                 x2_MR_STAR.LockedLocalMatrix(),
          T(1),  z1_MC_STAR.LocalMatrix() );
        Gemv
        ( TRANSPOSE,
          alpha, A12.LockedLocalMatrix(),
                 x1_MC_STAR.LockedLocalMatrix(),
          T(1),  z2_MR_STAR.LocalMatrix() );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDown
        ( xT_MR_STAR,  x0_MR_STAR,
                       x1_MR_STAR,
         /**********/ /**********/
          xB_MR_STAR,  x2_MR_STAR );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
LocalSymvRowAccumulateU
( T alpha, 
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR )
{
#ifndef RELEASE
    PushCallStack("internal::LocalSymvRowAccumulateU");
    if( A.Grid() != x_STAR_MC.Grid() ||
        x_STAR_MC.Grid() != x_STAR_MR.Grid() ||
        x_STAR_MR.Grid() != z_STAR_MC.Grid() ||
        z_STAR_MC.Grid() != z_STAR_MR.Grid() )
        throw std::logic_error
        ("{A,x,z} must be distributed over the same grid");
    if( x_STAR_MC.Height() != 1 || x_STAR_MR.Height() != 1 ||
        z_STAR_MC.Height() != 1 || z_STAR_MR.Height() != 1 )
        throw std::logic_error("Expected x and z to be row vectors");
    if( A.Height() != A.Width() || 
        A.Height() != x_STAR_MC.Width() ||
        A.Height() != x_STAR_MR.Width() ||
        A.Height() != z_STAR_MC.Width() ||
        A.Height() != z_STAR_MR.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalSymvRowAccumulateU: \n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  x[* ,MC] ~ " << x_STAR_MC.Height() << " x " 
                               << x_STAR_MC.Width() << "\n"
            << "  x[* ,MR] ~ " << x_STAR_MR.Height() << " x " 
                               << x_STAR_MR.Width() << "\n"
            << "  z[* ,MC] ~ " << z_STAR_MC.Height() << " x " 
                               << z_STAR_MC.Width() << "\n"
            << "  z[* ,MR] ~ " << z_STAR_MR.Height() << " x " 
                               << z_STAR_MR.Width() << "\n";
        throw std::logic_error( msg.str() );
    }
    if( x_STAR_MC.RowAlignment() != A.ColAlignment() ||
        x_STAR_MR.RowAlignment() != A.RowAlignment() ||
        z_STAR_MC.RowAlignment() != A.ColAlignment() ||
        z_STAR_MR.RowAlignment() != A.RowAlignment() )
        throw std::logic_error("Partial matrix distributions are misaligned");
#endif
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> A11(g), A12(g);
    DistMatrix<T> D11(g);

    DistMatrix<T,STAR,MC> x1_STAR_MC(g);
    DistMatrix<T,STAR,MR> 
        xL_STAR_MR(g), xR_STAR_MR(g),
        x0_STAR_MR(g), x1_STAR_MR(g), x2_STAR_MR(g);
    DistMatrix<T,STAR,MC> z1_STAR_MC(g);
    DistMatrix<T,STAR,MR> z1_STAR_MR(g), z2_STAR_MR(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const int ratio = std::max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*LocalSymvBlocksize<T>() );
                 
    LockedPartitionRight( x_STAR_MR,  xL_STAR_MR, xR_STAR_MR, 0 );
    while( xL_STAR_MR.Width() < x_STAR_MR.Width() )
    {
        LockedRepartitionRight
        ( xL_STAR_MR, /**/ xR_STAR_MR, 
          x0_STAR_MR, /**/ x1_STAR_MR, x2_STAR_MR );

        const int n0 = x0_STAR_MR.Width();
        const int n1 = x1_STAR_MR.Width();
        const int n2 = x2_STAR_MR.Width();
        LockedView( A11, A, n0, n0,    n1, n1 );
        LockedView( A12, A, n0, n0+n1, n1, n2 );
        LockedView( x1_STAR_MC, x_STAR_MC, 0, n0, 1, n1 );
        View( z1_STAR_MC, z_STAR_MC, 0, n0,    1, n1 );
        View( z1_STAR_MR, z_STAR_MR, 0, n0,    1, n1 );
        View( z2_STAR_MR, z_STAR_MR, 0, n0+n1, 1, n2 );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        // TODO: These diagonal block updates can be greatly improved
        D11 = A11;
        MakeTrapezoidal( LEFT, UPPER, 0, D11 );
        Gemv
        ( NORMAL, 
          alpha, D11.LockedLocalMatrix(), 
                 x1_STAR_MR.LockedLocalMatrix(),
          T(1),  z1_STAR_MC.LocalMatrix() );
        MakeTrapezoidal( LEFT, UPPER, 1, D11 );
        Gemv
        ( TRANSPOSE,
          alpha, D11.LockedLocalMatrix(),
                 x1_STAR_MC.LockedLocalMatrix(),
          T(1),  z1_STAR_MR.LocalMatrix() );

        Gemv
        ( NORMAL,
          alpha, A12.LockedLocalMatrix(),
                 x2_STAR_MR.LockedLocalMatrix(),
          T(1),  z1_STAR_MC.LocalMatrix() );
        Gemv
        ( TRANSPOSE,
          alpha, A12.LockedLocalMatrix(),
                 x1_STAR_MC.LockedLocalMatrix(),
          T(1),  z2_STAR_MR.LocalMatrix() );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionRight
        ( xL_STAR_MR,             /**/ xR_STAR_MR,
          x0_STAR_MR, x1_STAR_MR, /**/ x2_STAR_MR );
    }
    PopBlocksizeStack();
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
