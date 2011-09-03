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
#include "elemental/basic_internal.hpp"
using namespace std;
using namespace elemental;

// Template conventions:
//   G: general datatype
//
//   T: any ring, e.g., the (Gaussian) integers and the real/complex numbers
//   Z: representation of a real ring, e.g., the integers or real numbers
//   std::complex<Z>: representation of a complex ring, e.g. Gaussian integers
//                    or complex numbers
//
//   F: representation of real or complex number
//   R: representation of real number
//   std::complex<R>: representation of complex number

// Left Lower NORMAL (Non)Unit Trsm 
//   X := tril(L)^-1  X, or
//   X := trilu(L)^-1 X

// For large numbers of RHS's, e.g., width(X) >> p
template<typename F>
void
elemental::basic::internal::TrsmLLNLarge
( Diagonal diagonal,
  F alpha, const DistMatrix<F,MC,MR>& L,
                 DistMatrix<F,MC,MR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrsmLLNLarge");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmLLN: " << endl
            << "  L ~ " << L.Height() << " x " << L.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,MC,MR> XT(g),  X0(g), 
                        XB(g),  X1(g),
                                X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    // Start the algorithm
    basic::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        L21_MC_STAR.AlignWith( X2 );
        X1_STAR_MR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]

        // X1[* ,VR] := L11^-1[* ,* ] X1[* ,VR]
        basic::internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, diagonal, (F)1, L11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR] <- X1[* ,MR]
        L21_MC_STAR = L21;        // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, L21_MC_STAR, X1_STAR_MR, (F)1, X2 );
        //--------------------------------------------------------------------//
        L21_MC_STAR.FreeAlignments();
        X1_STAR_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// For medium numbers of RHS's, e.g., width(X) ~= p
template<typename F>
void
elemental::basic::internal::TrsmLLNMedium
( Diagonal diagonal,
  F alpha, const DistMatrix<F,MC,MR>& L,
                 DistMatrix<F,MC,MR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrsmLLNMedium");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmLLN: " << endl
            << "  L ~ " << L.Height() << " x " << L.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,MC,MR> XT(g),  X0(g), 
                        XB(g),  X1(g),
                                X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,MC,  STAR> L21_MC_STAR(g);
    DistMatrix<F,MR,  STAR> X1Trans_MR_STAR(g);

    // Start the algorithm
    basic::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        L21_MC_STAR.AlignWith( X2 );
        X1Trans_MR_STAR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11;                 // L11[* ,* ] <- L11[MC,MR]
        X1Trans_MR_STAR.TransposeFrom( X1 ); // X1[* ,MR] <- X1[MC,MR]

        // X1[* ,MR] := L11^-1[* ,* ] X1[* ,MR]
        // X1^T[MR,* ] := X1^T[MR,* ] L11^-T[* ,* ]
        basic::internal::LocalTrsm
        ( RIGHT, LOWER, TRANSPOSE, diagonal, 
          (F)1, L11_STAR_STAR, X1Trans_MR_STAR, checkIfSingular );

        X1.TransposeFrom( X1Trans_MR_STAR ); // X1[MC,MR]  <- X1[* ,MR]
        L21_MC_STAR = L21;                   // L21[MC,* ] <- L21[MC,MR]
        
        // X2[MC,MR] -= L21[MC,* ] X1[* ,MR]
        basic::internal::LocalGemm
        ( NORMAL, TRANSPOSE, (F)-1, L21_MC_STAR, X1Trans_MR_STAR, (F)1, X2 );
        //--------------------------------------------------------------------//
        L21_MC_STAR.FreeAlignments();
        X1Trans_MR_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

// For small numbers of RHS's, e.g., width(X) < p
template<typename F>
void
elemental::basic::internal::TrsmLLNSmall
( Diagonal diagonal,
  F alpha, const DistMatrix<F,VC,STAR>& L,
                 DistMatrix<F,VC,STAR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("basic::internal::TrsmLLNSmall");
    if( L.Grid() != X.Grid() )
        throw logic_error( "L and X must be distributed over the same grid." );
    if( L.Height() != L.Width() || L.Width() != X.Height() )
    {
        ostringstream msg;
        msg << "Nonconformal TrsmLLN: " << endl
            << "  L ~ " << L.Height() << " x " << L.Width() << endl
            << "  X ~ " << X.Height() << " x " << X.Width() << endl;
        throw logic_error( msg.str() );
    }
    if( L.ColAlignment() != X.ColAlignment() )
        throw std::logic_error("L and X are assumed to be aligned");
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,VC,STAR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,VC,STAR> XT(g),  X0(g), 
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    // Start the algorithm
    basic::Scal( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; // L11[* ,* ] <- L11[VC,* ]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VC,* ]

        // X1[* ,* ] := (L11[* ,* ])^-1 X1[* ,* ]
        basic::internal::LocalTrsm
        ( LEFT, LOWER, NORMAL, diagonal, 
          (F)1, L11_STAR_STAR, X1_STAR_STAR, checkIfSingular );

        // X2[VC,* ] -= L21[VC,* ] X1[* ,* ]
        basic::internal::LocalGemm
        ( NORMAL, NORMAL, (F)-1, L21, X1_STAR_STAR, (F)1, X2 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12, 
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionDown
        ( XT,  X0,
               X1,
         /**/ /**/
          XB,  X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template void elemental::basic::internal::TrsmLLNLarge
( Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X,
  bool checkIfSingular );
template void elemental::basic::internal::TrsmLLNMedium
( Diagonal diagonal,
  float alpha, 
  const DistMatrix<float,MC,MR>& L,
        DistMatrix<float,MC,MR>& X,
  bool checkIfSingular );
template void elemental::basic::internal::TrsmLLNSmall
( Diagonal diagonal,
  float alpha,
  const DistMatrix<float,VC,STAR>& L,
        DistMatrix<float,VC,STAR>& X,
  bool checkIfSingular );

template void elemental::basic::internal::TrsmLLNLarge
( Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X,
  bool checkIfSingular );
template void elemental::basic::internal::TrsmLLNMedium
( Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,MC,MR>& L,
        DistMatrix<double,MC,MR>& X,
  bool checkIfSingular );
template void elemental::basic::internal::TrsmLLNSmall
( Diagonal diagonal,
  double alpha, 
  const DistMatrix<double,VC,STAR>& L,
        DistMatrix<double,VC,STAR>& X,
  bool checkIfSingular );

#ifndef WITHOUT_COMPLEX
template void elemental::basic::internal::TrsmLLNLarge
( Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X,
  bool checkIfSingular );
template void elemental::basic::internal::TrsmLLNMedium
( Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,MC,MR>& L,
        DistMatrix<scomplex,MC,MR>& X,
  bool checkIfSingular );
template void elemental::basic::internal::TrsmLLNSmall
( Diagonal diagonal,
  scomplex alpha, 
  const DistMatrix<scomplex,VC,STAR>& L,
        DistMatrix<scomplex,VC,STAR>& X,
  bool checkIfSingular );

template void elemental::basic::internal::TrsmLLNLarge
( Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X,
  bool checkIfSingular );
template void elemental::basic::internal::TrsmLLNMedium
( Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,MC,MR>& L,
        DistMatrix<dcomplex,MC,MR>& X,
  bool checkIfSingular );
template void elemental::basic::internal::TrsmLLNSmall
( Diagonal diagonal,
  dcomplex alpha, 
  const DistMatrix<dcomplex,VC,STAR>& L,
        DistMatrix<dcomplex,VC,STAR>& X,
  bool checkIfSingular );
#endif

