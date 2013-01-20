/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRR2K_NTNT_HPP
#define BLAS_TRR2K_NTNT_HPP

namespace elem {
namespace internal {

// Distributed E := alpha (A B^{T/H} + C D^{T/H}) + beta E
template<typename T>
inline void
Trr2kNTNT
( UpperOrLower uplo,
  Orientation orientationOfB, Orientation orientationOfD,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kNTNT");
    if( E.Height() != E.Width()  || A.Width()  != C.Width()  ||
        A.Height() != E.Height() || C.Height() != E.Height() ||
        B.Height() != E.Width()  || D.Height() != E.Width()  ||
        A.Width()  != B.Width()  || C.Width()  != D.Width() )
        throw std::logic_error("Nonconformal Trr2kNTNT");
#endif
    const Grid& g = E.Grid();

    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);

    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g);
    DistMatrix<T> DL(g), DR(g),
                  D0(g), D1(g), D2(g);

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1AdjOrTrans_STAR_MR(g);
    DistMatrix<T,MC,  STAR> C1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> D1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > D1AdjOrTrans_STAR_MR(g);

    A1_MC_STAR.AlignWith( E );
    B1_VR_STAR.AlignWith( E );
    B1AdjOrTrans_STAR_MR.AlignWith( E );
    C1_MC_STAR.AlignWith( E );
    D1_VR_STAR.AlignWith( E );
    D1AdjOrTrans_STAR_MR.AlignWith( E );

    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    LockedPartitionRight( C, CL, CR, 0 );
    LockedPartitionRight( D, DL, DR, 0 );
    while( AL.Width() < A.Width() )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );
        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );
        LockedRepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );
        LockedRepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1;
        C1_MC_STAR = C1;
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
        ( uplo, 
          alpha, A1_MC_STAR, B1AdjOrTrans_STAR_MR, 
                 C1_MC_STAR, D1AdjOrTrans_STAR_MR,
          beta,  E );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( DL,     /**/ DR,
          D0, D1, /**/ D2 );
        SlideLockedPartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_TRR2K_NTNT_HPP
