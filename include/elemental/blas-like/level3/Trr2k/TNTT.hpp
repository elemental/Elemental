/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Distributed E := alpha (A^{T/H} B + C^{T/H} D^{T/H}) + beta E
template<typename T>
inline void
Trr2kTNTT
( UpperOrLower uplo,
  Orientation orientationOfA,
  Orientation orientationOfC, Orientation orientationOfD,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kTNTT");
    if( E.Height() != E.Width()  || A.Height() != C.Height() ||
        A.Width()  != E.Height() || C.Width()  != E.Height() ||
        B.Width()  != E.Width()  || D.Height() != E.Width()  ||
        A.Height() != B.Height() || C.Height() != D.Width() )
        throw std::logic_error("Nonconformal Trr2kNNTT");
#endif
    const Grid& g = E.Grid();

    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);
    DistMatrix<T> DL(g), DR(g),
                  D0(g), D1(g), D2(g);

    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<T,VR,  STAR> D1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > D1AdjOrTrans_STAR_MR(g);

    A1_STAR_MC.AlignWith( E );
    B1Trans_MR_STAR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1_VR_STAR.AlignWith( E );
    D1AdjOrTrans_STAR_MR.AlignWith( E );

    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
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
        LockedRepartitionDown
        ( BT,  B0,
         /**/ /**/
               B1,
          BB,  B2 );
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
        B1Trans_MR_STAR.TransposeFrom( B1 );
        D1_VR_STAR = D1;
        if( orientationOfD == ADJOINT )
            D1AdjOrTrans_STAR_MR.AdjointFrom( D1_VR_STAR );
        else
            D1AdjOrTrans_STAR_MR.TransposeFrom( D1_VR_STAR );
        LocalTrr2k
        ( uplo, orientationOfA, TRANSPOSE, orientationOfC,
          alpha, A1_STAR_MC, B1Trans_MR_STAR, 
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
        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );
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
