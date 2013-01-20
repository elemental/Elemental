/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRR2K_TTTT_HPP
#define BLAS_TRR2K_TTTT_HPP

namespace elem {
namespace internal {

// Distributed E := alpha (A^{T/H} B^{T/H} + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
Trr2kTTTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  Orientation orientationOfC, Orientation orientationOfD, 
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kTTTT");
    if( E.Height() != E.Width()  || A.Height() != C.Height() ||
        A.Width()  != E.Height() || C.Width()  != E.Height() ||
        B.Height() != E.Width()  || D.Height() != E.Width()  ||
        A.Height() != B.Width()  || C.Height() != D.Width() )
        throw std::logic_error("Nonconformal Trr2kTTTT");
#endif
    const Grid& g = E.Grid();

    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);

    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);
    DistMatrix<T> DL(g), DR(g),
                  D0(g), D1(g), D2(g);

    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1AdjOrTrans_STAR_MR(g);
    DistMatrix<T,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<T,VR,  STAR> D1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > D1AdjOrTrans_STAR_MR(g);

    A1_STAR_MC.AlignWith( E );
    B1_VR_STAR.AlignWith( E );
    B1AdjOrTrans_STAR_MR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1_VR_STAR.AlignWith( E );
    D1AdjOrTrans_STAR_MR.AlignWith( E );

    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    LockedPartitionDown
    ( C, CT,
         CB, 0 );
    LockedPartitionRight( D, DL, DR, 0 );
    while( AT.Height() < A.Height() )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );
        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );
        LockedRepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );
        LockedRepartitionRight
        ( DL, /**/ DR,
          D0, /**/ D1, D2 );

        //--------------------------------------------------------------------//
        A1_STAR_MC = A1;
        C1_STAR_MC = C1;
        B1_VR_STAR = B1;
        D1_VR_STAR = D1;
        if( orientationOfB == ADJOINT )
            B1AdjOrTrans_STAR_MR.AdjointFrom( B1_VR_STAR );
        else
            B1AdjOrTrans_STAR_MR.TransposeFrom( B1_VR_STAR );
        if( orientationOfD == ADJOINT )
            D1AdjOrTrans_STAR_MR.AdjointFrom( D1_VR_STAR );
        else
            D1AdjOrTrans_STAR_MR.TransposeFrom( D1_VR_STAR );
        LocalTrr2k
        ( uplo, orientationOfA, orientationOfC,
          alpha, A1_STAR_MC, B1AdjOrTrans_STAR_MR,
                 C1_STAR_MC, D1AdjOrTrans_STAR_MR,
          beta,  E );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( DL,     /**/ DR,
          D0, D1, /**/ D2 );
        SlideLockedPartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
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

#endif // ifndef BLAS_TRR2K_TTTT_HPP
