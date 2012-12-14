/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Left Upper (Conjugate)Transpose (Non)Unit Trsm
//   X := triu(U)^-T  X, 
//   X := triu(U)^-H  X,
//   X := triuu(U)^-T X, or
//   X := triuu(U)^-H X

// width(X) >> p
template<typename F>
inline void
TrsmLUTLarge
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& U, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmLUTLarge");
    if( orientation == NORMAL )
        throw std::logic_error("TrsmLUT expects a (Conjugate)Transpose option");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<F> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,STAR,MR  > X1_STAR_MR(g);
    DistMatrix<F,STAR,VR  > X1_STAR_VR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        U12_STAR_MC.AlignWith( X2 );
        X1_STAR_MR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        X1_STAR_VR    = X1;  // X1[* ,VR] <- X1[MC,MR]
        
        // X1[* ,VR] := U11^-[T/H][*,*] X1[* ,VR]
        LocalTrsm
        ( LEFT, UPPER, orientation, diag, F(1), U11_STAR_STAR, X1_STAR_VR,
          checkIfSingular );

        X1_STAR_MR  = X1_STAR_VR; // X1[* ,MR]  <- X1[* ,VR]
        X1          = X1_STAR_MR; // X1[MC,MR]  <- X1[* ,MR]
        U12_STAR_MC = U12;        // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])^(T/H) X1[* ,MR]
        //            = U12^(T/H)[MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, NORMAL, F(-1), U12_STAR_MC, X1_STAR_MR, F(1), X2 );
        //--------------------------------------------------------------------//
        U12_STAR_MC.FreeAlignments();
        X1_STAR_MR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

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

// width(X) ~= p
template<typename F>
inline void
TrsmLUTMedium
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F>& U, DistMatrix<F>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmLUTMedium");
    if( orientation == NORMAL )
        throw std::logic_error("TrsmLUT expects a (Conjugate)Transpose option");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<F> XT(g),  X0(g),
                  XB(g),  X1(g),
                          X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,MC  > U12_STAR_MC(g);
    DistMatrix<F,MR,  STAR> X1AdjOrTrans_MR_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        U12_STAR_MC.AlignWith( X2 );
        X1AdjOrTrans_MR_STAR.AlignWith( X2 );
        //--------------------------------------------------------------------//
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[MC,MR]
        // X1[* ,VR] <- X1[MC,MR]
        if( orientation == TRANSPOSE )
            X1AdjOrTrans_MR_STAR.TransposeFrom( X1 );
        else
            X1AdjOrTrans_MR_STAR.AdjointFrom( X1 );
        
        // X1[* ,MR] := U11^-[T/H][*,*] X1[* ,MR]
        //
        // X1^[T/H][MR,* ] := X1^[T/H][MR,* ] U11^-1[* ,* ]
        LocalTrsm
        ( RIGHT, UPPER, NORMAL, diag, 
          F(1), U11_STAR_STAR, X1AdjOrTrans_MR_STAR, checkIfSingular );

        if( orientation == TRANSPOSE )
            X1.TransposeFrom( X1AdjOrTrans_MR_STAR );
        else
            X1.AdjointFrom( X1AdjOrTrans_MR_STAR );

        U12_STAR_MC = U12;        // U12[* ,MC] <- U12[MC,MR]

        // X2[MC,MR] -= (U12[* ,MC])^[T/H] X1[* ,MR]
        //            = U12^[T/H][MC,*] X1[* ,MR]
        LocalGemm
        ( orientation, orientation, 
          F(-1), U12_STAR_MC, X1AdjOrTrans_MR_STAR, F(1), X2 );
        //--------------------------------------------------------------------//
        U12_STAR_MC.FreeAlignments();
        X1AdjOrTrans_MR_STAR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

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

// width(X) << p
template<typename F>
inline void
TrsmLUTSmall
( Orientation orientation, UnitOrNonUnit diag,
  F alpha, const DistMatrix<F,STAR,VR>& U, DistMatrix<F,VR,STAR>& X,
  bool checkIfSingular )
{
#ifndef RELEASE
    PushCallStack("internal::TrsmLUTSmall");
    if( U.Grid() != X.Grid() )
        throw std::logic_error
        ("U and X must be distributed over the same grid");
    if( orientation == NORMAL )
        throw std::logic_error("TrsmLUT expects a (Conjugate)Transpose option");
    if( U.Height() != U.Width() || U.Height() != X.Height() )
    {
        std::ostringstream msg;
        msg << "Nonconformal TrsmLUT: \n"
            << "  U ~ " << U.Height() << " x " << U.Width() << "\n"
            << "  X ~ " << X.Height() << " x " << X.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
    if( U.RowAlignment() != X.ColAlignment() )
        throw std::logic_error("U and X are assumed to be aligned");
#endif
    const Grid& g = U.Grid();

    // Matrix views
    DistMatrix<F,STAR,VR> 
        UTL(g), UTR(g),  U00(g), U01(g), U02(g),
        UBL(g), UBR(g),  U10(g), U11(g), U12(g),
                         U20(g), U21(g), U22(g);

    DistMatrix<F,VR,STAR> XT(g),  X0(g),
                          XB(g),  X1(g),
                                  X2(g);

    // Temporary distributions
    DistMatrix<F,STAR,STAR> U11_STAR_STAR(g); 
    DistMatrix<F,STAR,STAR> X1_STAR_STAR(g);

    // Start the algorithm
    Scale( alpha, X );
    LockedPartitionDownDiagonal
    ( U, UTL, UTR,
         UBL, UBR, 0 );
    PartitionDown
    ( X, XT,
         XB, 0 );
    while( XB.Height() > 0 )
    {
        LockedRepartitionDownDiagonal
        ( UTL, /**/ UTR,  U00, /**/ U01, U02,
         /*************/ /******************/
               /**/       U10, /**/ U11, U12,
          UBL, /**/ UBR,  U20, /**/ U21, U22 );

        RepartitionDown
        ( XT,  X0,
         /**/ /**/
               X1,
          XB,  X2 );

        //--------------------------------------------------------------------//
        U11_STAR_STAR = U11; // U11[* ,* ] <- U11[* ,VR]
        X1_STAR_STAR = X1;   // X1[* ,* ] <- X1[VR,* ]
        
        // X1[* ,* ] := U11^-[T/H][* ,* ] X1[* ,* ]
        LocalTrsm
        ( LEFT, UPPER, orientation, diag,
          F(1), U11_STAR_STAR, X1_STAR_STAR, checkIfSingular );

        X1 = X1_STAR_STAR;

        // X2[VR,* ] -= U12[* ,VR]^[T/H] X1[* ,* ]
        LocalGemm( orientation, NORMAL, F(-1), U12, X1_STAR_STAR, F(1), X2 );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( UTL, /**/ UTR,   U00, U01, /**/ U02,
               /**/        U10, U11, /**/ U12,
         /*************/  /******************/
          UBL, /**/ UBR,   U20, U21, /**/ U22 );

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

} // namespace internal
} // namespace elem
