/*
   Copyright (c) 2009-2011, Jack Poulson
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

template<typename T>
inline void
elemental::basic::internal::LocalHemvColAccumulateL
( T alpha, 
  const DistMatrix<T,MC,MR  >& A,
  const DistMatrix<T,MC,STAR>& x_MC_STAR,
  const DistMatrix<T,MR,STAR>& x_MR_STAR,
        DistMatrix<T,MC,STAR>& z_MC_STAR,
        DistMatrix<T,MR,STAR>& z_MR_STAR )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalHemvColAccumulateL");
    if( A.Grid() != x_MC_STAR.Grid() ||
        x_MC_STAR.Grid() != x_MR_STAR.Grid() ||
        x_MR_STAR.Grid() != z_MC_STAR.Grid() ||
        z_MC_STAR.Grid() != z_MR_STAR.Grid() )
        throw std::logic_error
        ("{A,x,z} must be distributed over the same grid");
    if( x_MC_STAR.Width() != 1 || x_MR_STAR.Width() != 1 ||
        z_MC_STAR.Width() != 1 || z_MR_STAR.Width() != 1  )
        throw std::logic_error("Expected x and z to be column vectors");
    if( A.Height() != A.Width() || 
        A.Height() != x_MC_STAR.Height() ||
        A.Height() != x_MR_STAR.Height() ||
        A.Height() != z_MC_STAR.Height() ||
        A.Height() != z_MR_STAR.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalHemvColAccumulateL: \n"
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
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,MC,STAR> 
        xT_MC_STAR(g),  x0_MC_STAR(g),
        xB_MC_STAR(g),  x1_MC_STAR(g),
                        x2_MC_STAR(g);

    DistMatrix<T,MR,STAR> 
        xT_MR_STAR(g),  x0_MR_STAR(g),
        xB_MR_STAR(g),  x1_MR_STAR(g),
                        x2_MR_STAR(g);

    DistMatrix<T,MC,STAR> 
        zT_MC_STAR(g),  z0_MC_STAR(g),
        zB_MC_STAR(g),  z1_MC_STAR(g),
                        z2_MC_STAR(g);

    DistMatrix<T,MR,STAR> 
        zT_MR_STAR(g),  z0_MR_STAR(g),
        zB_MR_STAR(g),  z1_MR_STAR(g),
                        z2_MR_STAR(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*LocalHemvBlocksize<T>() );
                 
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionDown
    ( x_MC_STAR, xT_MC_STAR,
                 xB_MC_STAR, 0 );
    LockedPartitionDown
    ( x_MR_STAR, xT_MR_STAR,
                 xB_MR_STAR, 0 );
    PartitionDown
    ( z_MC_STAR, zT_MC_STAR,
                 zB_MC_STAR, 0 );
    PartitionDown
    ( z_MR_STAR, zT_MR_STAR,
                 zB_MR_STAR, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionDown
        ( xT_MC_STAR,  x0_MC_STAR,
         /**********/ /**********/
                       x1_MC_STAR,
          xB_MC_STAR,  x2_MC_STAR );

        LockedRepartitionDown
        ( xT_MR_STAR,  x0_MR_STAR,
         /**********/ /**********/
                       x1_MR_STAR,
          xB_MR_STAR,  x2_MR_STAR );

        RepartitionDown
        ( zT_MC_STAR,  z0_MC_STAR,
         /**********/ /**********/
                       z1_MC_STAR,
          zB_MC_STAR,  z2_MC_STAR );

        RepartitionDown
        ( zT_MR_STAR,  z0_MR_STAR,
         /**********/ /**********/
                       z1_MR_STAR,
          zB_MR_STAR,  z2_MR_STAR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( LEFT, LOWER );
        basic::Gemv
        ( NORMAL, 
          alpha, D11.LockedLocalMatrix(), 
                 x1_MR_STAR.LockedLocalMatrix(),
          (T)1,  z1_MC_STAR.LocalMatrix() );
        D11.MakeTrapezoidal( LEFT, LOWER, -1 );

        basic::Gemv
        ( ADJOINT,
          alpha, D11.LockedLocalMatrix(),
                 x1_MC_STAR.LockedLocalMatrix(),
          (T)1,  z1_MR_STAR.LocalMatrix() );

        basic::Gemv
        ( NORMAL,
          alpha, A21.LockedLocalMatrix(),
                 x1_MR_STAR.LockedLocalMatrix(),
          (T)1,  z2_MC_STAR.LocalMatrix() );

        basic::Gemv
        ( ADJOINT,
          alpha, A21.LockedLocalMatrix(),
                 x2_MC_STAR.LockedLocalMatrix(),
          (T)1,  z1_MR_STAR.LocalMatrix() );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionDown
        ( xT_MC_STAR,  x0_MC_STAR,
                       x1_MC_STAR,
         /**********/ /**********/
          xB_MC_STAR,  x2_MC_STAR );

        SlideLockedPartitionDown
        ( xT_MR_STAR,  x0_MR_STAR,
                       x1_MR_STAR,
         /**********/ /**********/
          xB_MR_STAR,  x2_MR_STAR );
        
        SlidePartitionDown
        ( zT_MC_STAR,  z0_MC_STAR,
                       z1_MC_STAR,
         /**********/ /**********/
          zB_MC_STAR,  z2_MC_STAR );

        SlidePartitionDown
        ( zT_MR_STAR,  z0_MR_STAR,
                       z1_MR_STAR,
         /**********/ /**********/
          zB_MR_STAR,  z2_MR_STAR );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename T>
inline void
elemental::basic::internal::LocalHemvRowAccumulateL
( T alpha, 
  const DistMatrix<T,MC,  MR>& A,
  const DistMatrix<T,STAR,MC>& x_STAR_MC,
  const DistMatrix<T,STAR,MR>& x_STAR_MR,
        DistMatrix<T,STAR,MC>& z_STAR_MC,
        DistMatrix<T,STAR,MR>& z_STAR_MR )
{
#ifndef RELEASE
    PushCallStack("basic::internal::LocalHemvRowAccumulateL");
    if( A.Grid() != x_STAR_MC.Grid() ||
        x_STAR_MC.Grid() != x_STAR_MR.Grid() ||
        x_STAR_MR.Grid() != z_STAR_MC.Grid() ||
        z_STAR_MC.Grid() != z_STAR_MR.Grid() )
        throw std::logic_error
        ("{A,x,z} must be distributed over the same grid");
    if( x_STAR_MC.Height() != 1 || x_STAR_MR.Height() != 1 ||
        z_STAR_MC.Height() != 1 || z_STAR_MR.Height() != 1  )
        throw std::logic_error("Expected x and z to be row vectors");
    if( A.Height() != A.Width() || 
        A.Height() != x_STAR_MC.Width() ||
        A.Height() != x_STAR_MR.Width() ||
        A.Height() != z_STAR_MC.Width() ||
        A.Height() != z_STAR_MR.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal LocalHemvRowAccumulateL: \n"
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
    DistMatrix<T,MC,MR> 
        ATL(g), ATR(g),  A00(g), A01(g), A02(g),
        ABL(g), ABR(g),  A10(g), A11(g), A12(g),
                         A20(g), A21(g), A22(g);

    DistMatrix<T,MC,MR> D11(g);

    DistMatrix<T,STAR,MC> 
        xL_STAR_MC(g), xR_STAR_MC(g),
        x0_STAR_MC(g), x1_STAR_MC(g), x2_STAR_MC(g);

    DistMatrix<T,STAR,MR> 
        xL_STAR_MR(g), xR_STAR_MR(g),
        x0_STAR_MR(g), x1_STAR_MR(g), x2_STAR_MR(g);

    DistMatrix<T,STAR,MC> 
        zL_STAR_MC(g), zR_STAR_MC(g),
        z0_STAR_MC(g), z1_STAR_MC(g), z2_STAR_MC(g);

    DistMatrix<T,STAR,MR> 
        zL_STAR_MR(g), zR_STAR_MR(g),
        z0_STAR_MR(g), z1_STAR_MR(g), z2_STAR_MR(g);

    // We want our local gemvs to be of width blocksize, so we will 
    // temporarily change to max(r,c) times the current blocksize
    const int ratio = max( g.Height(), g.Width() );
    PushBlocksizeStack( ratio*LocalHemvBlocksize<T>() );
                 
    LockedPartitionDownDiagonal
    ( A, ATL, ATR,
         ABL, ABR, 0 );
    LockedPartitionRight( x_STAR_MC, xL_STAR_MC, xR_STAR_MC, 0 );
    LockedPartitionRight( x_STAR_MR, xL_STAR_MR, xR_STAR_MR, 0 );
    PartitionRight( z_STAR_MC, zL_STAR_MC, zR_STAR_MC, 0 );
    PartitionRight( z_STAR_MR, zL_STAR_MR, zR_STAR_MR, 0 );
    while( ATL.Height() < A.Height() )
    {
        LockedRepartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, /**/ A01, A02,
         /*************/ /******************/
               /**/       A10, /**/ A11, A12,
          ABL, /**/ ABR,  A20, /**/ A21, A22 );

        LockedRepartitionRight
        ( xL_STAR_MC, /**/ xR_STAR_MC, 
          x0_STAR_MC, /**/ x1_STAR_MC, x2_STAR_MC );

        LockedRepartitionRight
        ( xL_STAR_MR, /**/ xR_STAR_MR, 
          x0_STAR_MR, /**/ x1_STAR_MR, x2_STAR_MR );

        RepartitionRight
        ( zL_STAR_MC, /**/ zR_STAR_MC,
          z0_STAR_MC, /**/ z1_STAR_MC, z2_STAR_MC );

        RepartitionRight
        ( zL_STAR_MR, /**/ zR_STAR_MR,
          z0_STAR_MR, /**/ z1_STAR_MR, z2_STAR_MR );

        D11.AlignWith( A11 );
        //--------------------------------------------------------------------//
        D11 = A11;
        D11.MakeTrapezoidal( LEFT, LOWER );
        basic::Gemv
        ( NORMAL, 
          alpha, D11.LockedLocalMatrix(), 
                 x1_STAR_MR.LockedLocalMatrix(),
          (T)1,  z1_STAR_MC.LocalMatrix() );
        D11.MakeTrapezoidal( LEFT, LOWER, -1 );

        basic::Gemv
        ( ADJOINT,
          alpha, D11.LockedLocalMatrix(),
                 x1_STAR_MC.LockedLocalMatrix(),
          (T)1,  z1_STAR_MR.LocalMatrix() );

        basic::Gemv
        ( NORMAL,
          alpha, A21.LockedLocalMatrix(),
                 x1_STAR_MR.LockedLocalMatrix(),
          (T)1,  z2_STAR_MC.LocalMatrix() );

        basic::Gemv
        ( ADJOINT,
          alpha, A21.LockedLocalMatrix(),
                 x2_STAR_MC.LockedLocalMatrix(),
          (T)1,  z1_STAR_MR.LocalMatrix() );
        //--------------------------------------------------------------------//
        D11.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( ATL, /**/ ATR,  A00, A01, /**/ A02,
               /**/       A10, A11, /**/ A12,
         /*************/ /******************/
          ABL, /**/ ABR,  A20, A21, /**/ A22 );

        SlideLockedPartitionRight
        ( xL_STAR_MC,             /**/ xR_STAR_MC,
          x0_STAR_MC, x1_STAR_MC, /**/ x2_STAR_MC );
        
        SlideLockedPartitionRight
        ( xL_STAR_MR,             /**/ xR_STAR_MR,
          x0_STAR_MR, x1_STAR_MR, /**/ x2_STAR_MR );

        SlidePartitionRight
        ( zL_STAR_MC,             /**/ zR_STAR_MC,
          z0_STAR_MC, z1_STAR_MC, /**/ z2_STAR_MC );

        SlidePartitionRight
        ( zL_STAR_MR,             /**/ zR_STAR_MR,
          z0_STAR_MR, z1_STAR_MR, /**/ z2_STAR_MR );
    }

    PopBlocksizeStack();

#ifndef RELEASE
    PopCallStack();
#endif
}
