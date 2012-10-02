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
TrsvLN
( UnitOrNonUnit diag, 
  const DistMatrix<F>& L, 
        DistMatrix<F>& x )
{
#ifndef RELEASE
    PushCallStack("internal::TrsvLN");
    if( L.Grid() != x.Grid() )
        throw std::logic_error("{L,x} must be distributed over the same grid");
    if( L.Height() != L.Width() )
        throw std::logic_error("L must be square");
    if( x.Width() != 1 && x.Height() != 1 )
        throw std::logic_error("x must be a vector");
    const int xLength = ( x.Width() == 1 ? x.Height() : x.Width() );
    if( L.Width() != xLength )
        throw std::logic_error("Nonconformal TrsvLN");
#endif
    const Grid& g = L.Grid();

    if( x.Width() == 1 )
    {
        // Matrix views 
        DistMatrix<F> 
            LTL(g), LTR(g),  L00(g), L01(g), L02(g),
            LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                             L20(g), L21(g), L22(g);

        DistMatrix<F> 
            xT(g),  x0(g),
            xB(g),  x1(g),
                    x2(g);

        // Temporary distributions
        DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
        DistMatrix<F,STAR,STAR> x1_STAR_STAR(g);
        DistMatrix<F,MR,  STAR> x1_MR_STAR(g);
        DistMatrix<F,MC,  STAR> z2_MC_STAR(g);

        // Start the algorithm
        LockedPartitionDownDiagonal
        ( L, LTL, LTR,
             LBL, LBR, 0 );
        PartitionDown
        ( x, xT,
             xB, 0 );
        while( xB.Height() > 0 )
        {
            LockedRepartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            RepartitionDown
            ( xT,  x0,
             /**/ /**/
                   x1,
              xB,  x2 );

            x1_MR_STAR.AlignWith( L21 );
            z2_MC_STAR.AlignWith( L21 );
            z2_MC_STAR.ResizeTo( x2.Height(), 1 );
            //----------------------------------------------------------------//
            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            Trsv
            ( LOWER, NORMAL, diag,
              L11_STAR_STAR.LockedLocalMatrix(),
              x1_STAR_STAR.LocalMatrix() );
            x1 = x1_STAR_STAR;

            x1_MR_STAR = x1_STAR_STAR;
            Gemv
            ( NORMAL, F(-1), 
              L21.LockedLocalMatrix(), 
              x1_MR_STAR.LockedLocalMatrix(),
              F(0), z2_MC_STAR.LocalMatrix() );
            x2.SumScatterUpdate( F(1), z2_MC_STAR );
            //----------------------------------------------------------------//
            x1_MR_STAR.FreeAlignments();
            z2_MC_STAR.FreeAlignments();

            SlideLockedPartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

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
        DistMatrix<F> 
            LTL(g), LTR(g),  L00(g), L01(g), L02(g),
            LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                             L20(g), L21(g), L22(g);

        DistMatrix<F> 
            xL(g), xR(g),
            x0(g), x1(g), x2(g);

        // Temporary distributions
        DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
        DistMatrix<F,STAR,STAR> x1_STAR_STAR(g);
        DistMatrix<F,STAR,MR  > x1_STAR_MR(g);
        DistMatrix<F,STAR,MC  > z2_STAR_MC(g);
        DistMatrix<F,MR,  MC  > z2_MR_MC(g);
        DistMatrix<F> z2(g);

        // Start the algorithm
        LockedPartitionDownDiagonal
        ( L, LTL, LTR,
             LBL, LBR, 0 );
        PartitionRight( x,  xL, xR, 0 );
        while( xR.Width() > 0 )
        {
            LockedRepartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, /**/ L01, L02,
             /*************/ /******************/
                   /**/       L10, /**/ L11, L12,
              LBL, /**/ LBR,  L20, /**/ L21, L22 );

            RepartitionRight
            ( xL, /**/ xR,
              x0, /**/ x1, x2 );

            x1_STAR_MR.AlignWith( L21 );
            z2_STAR_MC.AlignWith( L21 );
            z2.AlignWith( x2 );
            z2_STAR_MC.ResizeTo( 1, x2.Width() );
            //----------------------------------------------------------------//
            x1_STAR_STAR = x1;
            L11_STAR_STAR = L11;
            Trsv
            ( LOWER, NORMAL, diag,
              L11_STAR_STAR.LockedLocalMatrix(),
              x1_STAR_STAR.LocalMatrix() );
            x1 = x1_STAR_STAR;

            x1_STAR_MR = x1_STAR_STAR;
            Gemv
            ( NORMAL, F(-1), 
              L21.LockedLocalMatrix(), 
              x1_STAR_MR.LockedLocalMatrix(),
              F(0), z2_STAR_MC.LocalMatrix() );
            z2_MR_MC.SumScatterFrom( z2_STAR_MC );
            z2 = z2_MR_MC;
            Axpy( F(1), z2, x2 );
            //----------------------------------------------------------------//
            x1_STAR_MR.FreeAlignments();
            z2_STAR_MC.FreeAlignments();
            z2.FreeAlignments(); 

            SlideLockedPartitionDownDiagonal
            ( LTL, /**/ LTR,  L00, L01, /**/ L02,
                   /**/       L10, L11, /**/ L12,
             /*************/ /******************/
              LBL, /**/ LBR,  L20, L21, /**/ L22 );

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
