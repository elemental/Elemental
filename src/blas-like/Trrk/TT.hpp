/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef ELEM_TRRK_TT_HPP
#define ELEM_TRRK_TT_HPP

namespace elem {
namespace internal {

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::TrrkTN");
        if( C.Height() != C.Width() || A.Width() != C.Height() || 
            B.Height() != C.Width() || A.Height() != B.Width() )
            LogicError("Nonconformal TrrkTN");
        if( orientationOfA == NORMAL || orientationOfB == NORMAL )
            LogicError("Orientations must be TRANSPOSE or ADJOINT");
    )
    const Grid& g = C.Grid();

    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);

    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);

    A1_STAR_MC.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );

    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    while( AT.Width() < A.Height() )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );
        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        //--------------------------------------------------------------------//
        A1_STAR_MC = A1;
        B1_VR_STAR = B1;
        B1_VR_STAR.TransposePartialColAllGather
        ( B1Trans_STAR_MR, (orientationOfB==ADJOINT) );
        LocalTrrk
        ( uplo, orientationOfA,
          alpha, A1_STAR_MC, B1Trans_STAR_MR, beta, C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_TRRK_TT_HPP
