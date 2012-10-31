/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental.

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    - Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright notice,
      this list of conditions and the following disclaimer in the documentation
      and/or other materials provided with the distribution.

    - Neither the name of the owner nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
*/

namespace elem {
namespace internal {

template<typename F>
inline void
TrsvUT
( Orientation orientation, UnitOrNonUnit diag, 
  const DistMatrix<F>& U, DistMatrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("internal::TrsvUT");
    if( U.Grid() != x.Grid() )
        throw std::logic_error("{U,x} must be distributed over the same grid");
    if( orientation == NORMAL )
        throw std::logic_error("TrsvUT expects a (conjugate-)transpose option");
    if( U.Height() != U.Width() )
        throw std::logic_error("U must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( U.Width() != xLength )
        throw std::logic_error("Nonconformal TrsvUT");
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
        Zeros( x.Height(), 1, z_MR_STAR );

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

            const int n0 = x0.Height();
            const int n1 = x1.Height();
            const int n2 = x2.Height();
            U11.LockedView( U, n0, n0,    n1, n1 );
            U12.LockedView( U, n0, n0+n1, n1, n2 );
            z1_MR_STAR.View( z_MR_STAR, n0,    0, n1, 1 );
            z2_MR_STAR.View( z_MR_STAR, n0+n1, 0, n2, 1 );

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
              U11_STAR_STAR.LockedLocalMatrix(),
              x1_STAR_STAR.LocalMatrix() );
            x1 = x1_STAR_STAR;

            x1_MC_STAR = x1_STAR_STAR;
            Gemv
            ( orientation, F(-1), 
              U12.LockedLocalMatrix(), 
              x1_MC_STAR.LockedLocalMatrix(),
              F(1), z2_MR_STAR.LocalMatrix() );
            //----------------------------------------------------------------//
            x1_MC_STAR.FreeAlignments();
            z1.FreeAlignments();

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
        Zeros( 1, x.Width(), z_STAR_MR );

        // Start the algorithm
        PartitionRight( x,  xL, xR, 0 );
        while( xR.Width() > 0 )
        {
            RepartitionRight
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );

            const int n0 = x0.Width();
            const int n1 = x1.Width();
            const int n2 = x2.Width();
            U11.LockedView( U, n0, n0,    n1, n1 );
            U12.LockedView( U, n0, n0+n1, n1, n2 );
            z1_STAR_MR.View( z_STAR_MR, 0, n0,    1, n1 );
            z2_STAR_MR.View( z_STAR_MR, 0, n0+n1, 1, n2 );

            x1_STAR_MC.AlignWith( U12 );
            //----------------------------------------------------------------//
            if( x2.Width() != 0 )
                x1.SumScatterUpdate( F(1), z1_STAR_MR );

            x1_STAR_STAR = x1;
            U11_STAR_STAR = U11;
            Trsv
            ( UPPER, orientation, diag,
              U11_STAR_STAR.LockedLocalMatrix(),
              x1_STAR_STAR.LocalMatrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MC = x1_STAR_STAR;
            Gemv
            ( orientation, F(-1), 
              U12.LockedLocalMatrix(), 
              x1_STAR_MC.LockedLocalMatrix(),
              F(1), z2_STAR_MR.LocalMatrix() );
            //----------------------------------------------------------------//
            x1_STAR_MC.FreeAlignments();

            SlidePartitionRight
            ( xL,     /**/ xR,
              x0, x1, /**/ x2 );
        }
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
