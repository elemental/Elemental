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

// Right Lower (Conjugate)Transpose (Non)Unit Trsm
//   X := X tril(L)^-T,
//   X := X tril(L)^-H,
//   X := X trilu(L)^-T, or
//   X := X trilu(L)^-H
template<typename F>
inline void
TrsmRLT
( Orientation orientation, 
  UnitOrNonUnit diag,
  F alpha, 
  const DistMatrix<F,MC,MR>& L,
        DistMatrix<F,MC,MR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmRLT");
    if( L.Grid() != X.Grid() )
        throw std::logic_error
        ("L and X must be distributed over the same grid");
    if( orientation == NORMAL )
        throw std::logic_error("TrsmRLT expects a (Conjugate)Transpose option");
    if( L.Height() != L.Width() || X.Width() != L.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrsmRLT: \n"
            << "  L ~ " << L.Height() << " x " << L.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = L.Grid();

    // Matrix views
    DistMatrix<F,MC,MR> 
        LTL(g), LTR(g),  L00(g), L01(g), L02(g),
        LBL(g), LBR(g),  L10(g), L11(g), L12(g),
                         L20(g), L21(g), L22(g);

    DistMatrix<F,MC,MR> XL(g), XR(g),
                        X0(g), X1(g), X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> L11_STAR_STAR(g);
    DistMatrix<F,VR,  STAR> L21_VR_STAR(g);
    DistMatrix<F,STAR,MR  > L21AdjOrTrans_STAR_MR(g);
    DistMatrix<F,VC,  STAR> X1_VC_STAR(g);
    DistMatrix<F,STAR,MC  > X1Trans_STAR_MC(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( L, LTL, LTR,
         LBL, LBR, 0 );
    PartitionRight( X, XL, XR, 0 );
    while( XR.Width() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, /**/ L01, L02,
         /*************/ /******************/
               /**/       L10, /**/ L11, L12,
          LBL, /**/ LBR,  L20, /**/ L21, L22 );

        RepartitionRight
        ( XL, /**/     XR,
          X0, /**/ X1, X2 );

        X1_VC_STAR.AlignWith( X2 );
        X1Trans_STAR_MC.AlignWith( X2 );
        L21_VR_STAR.AlignWith( X2 );
        L21AdjOrTrans_STAR_MR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        L11_STAR_STAR = L11; 
        X1_VC_STAR = X1;  
        
        LocalTrsm
        ( RIGHT, LOWER, orientation, diag, 
          (F)1, L11_STAR_STAR, X1_VC_STAR, checkIfSingular );

        X1Trans_STAR_MC.TransposeFrom( X1_VC_STAR );
        X1.TransposeFrom( X1Trans_STAR_MC );
        L21_VR_STAR = L21;
        if( orientation == ADJOINT )
            L21AdjOrTrans_STAR_MR.AdjointFrom( L21_VR_STAR );
        else
            L21AdjOrTrans_STAR_MR.TransposeFrom( L21_VR_STAR );

        // X2[MC,MR] -= X1[MC,*] (L21[MR,*])^(T/H)
        //            = X1^T[* ,MC] (L21^(T/H))[*,MR]
        LocalGemm
        ( TRANSPOSE, NORMAL, 
          (F)-1, X1Trans_STAR_MC, L21AdjOrTrans_STAR_MR, (F)1, X2 );
        //--------------------------------------------------------------------//
        X1_VC_STAR.FreeAlignments();
        X1Trans_STAR_MC.FreeAlignments();
        L21_VR_STAR.FreeAlignments();
        L21AdjOrTrans_STAR_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( LTL, /**/ LTR,  L00, L01, /**/ L02,
               /**/       L10, L11, /**/ L12,
         /*************/ /******************/
          LBL, /**/ LBR,  L20, L21, /**/ L22 );

        SlidePartitionRight
        ( XL,     /**/ XR,
          X0, X1, /**/ X2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
