/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Distributed E := alpha (A B^{T/H} + C^{T/H} D) + beta E
template<typename T>
inline void
Trr2kNTTN
( UpperOrLower uplo,
  Orientation orientationOfB, Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
           const DistMatrix<T,MC,MR>& C, const DistMatrix<T,MC,MR>& D,
  T beta,        DistMatrix<T,MC,MR>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kNTTN");
    if( E.Height() != E.Width()  || A.Width()  != C.Height() ||
        A.Height() != E.Height() || C.Width()  != E.Height() ||
        B.Height() != E.Width()  || D.Width()  != E.Width()  ||
        A.Width()  != B.Width()  || C.Height() != D.Height() )
        throw std::logic_error("Nonconformal Trr2kNTTN");
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,MR> AL(g), AR(g),
                        A0(g), A1(g), A2(g);
    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    DistMatrix<T,MC,MR> CT(g),  C0(g),
                        CB(g),  C1(g),
                                C2(g);
    DistMatrix<T,MC,MR> DT(g),  D0(g),
                        DB(g),  D1(g),
                                D2(g);

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1AdjOrTrans_STAR_MR(g);
    DistMatrix<T,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( E );
    B1_VR_STAR.AlignWith( E );
    B1AdjOrTrans_STAR_MR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1Trans_MR_STAR.AlignWith( E );

    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    LockedPartitionDown
    ( C, CT,
         CB, 0 );
    LockedPartitionDown
    ( D, DT,
         DB, 0 );
    while( AL.Width() < A.Width() )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );
        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );
        LockedRepartitionDown
        ( CT,  C0,
         /**/ /**/
               C1,
          CB,  C2 );
        LockedRepartitionDown
        ( DT,  D0,
         /**/ /**/
               D1,
          DB,  D2 );

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1;
        C1_STAR_MC = C1;
        B1_VR_STAR = B1;
        if( orientationOfB == ADJOINT )
            B1AdjOrTrans_STAR_MR.AdjointFrom( B1_VR_STAR );
        else
            B1AdjOrTrans_STAR_MR.TransposeFrom( B1_VR_STAR );
        D1Trans_MR_STAR.TransposeFrom( D1 );
        LocalTrr2k 
        ( uplo, orientationOfC, TRANSPOSE,
          alpha, A1_MC_STAR, B1AdjOrTrans_STAR_MR,
                 C1_STAR_MC, D1Trans_MR_STAR,
          beta,  E );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
        SlideLockedPartitionDown
        ( CT,  C0,
               C1,
         /**/ /**/
          CB,  C2 );
        SlideLockedPartitionDown
        ( DT,  D0,
               D1,
         /**/ /**/
          DB,  D2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
