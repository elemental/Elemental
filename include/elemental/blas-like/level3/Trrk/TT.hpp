/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRRK_TT_HPP
#define BLAS_TRRK_TT_HPP

namespace elem {
namespace internal {

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
inline void
TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::TrrkTN");
    if( C.Height() != C.Width() ||
        A.Width() != C.Height() || 
        B.Height() != C.Width() ||
        A.Height() != B.Width() )
        throw std::logic_error("Nonconformal TrrkTN");
    if( orientationOfA == NORMAL || orientationOfB == NORMAL )
        throw std::logic_error("Orientations must be TRANSPOSE or ADJOINT");
#endif
    const Grid& g = C.Grid();

    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);

    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1AdjOrTrans_STAR_MR(g);

    A1_STAR_MC.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1AdjOrTrans_STAR_MR.AlignWith( C );

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
        if( orientationOfB == ADJOINT )
            B1AdjOrTrans_STAR_MR.AdjointFrom( B1_VR_STAR );
        else
            B1AdjOrTrans_STAR_MR.TransposeFrom( B1_VR_STAR );
        LocalTrrk
        ( uplo, orientationOfA,
          alpha, A1_STAR_MC, B1AdjOrTrans_STAR_MR, beta, C );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TRRK_TT_HPP
