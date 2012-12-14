/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Distributed E := alpha (A^{T/H} B + C^{T/H} D) + beta E
template<typename T>
inline void
Trr2kTNTN
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfC,
  T alpha, const DistMatrix<T,MC,MR>& A, const DistMatrix<T,MC,MR>& B,
           const DistMatrix<T,MC,MR>& C, const DistMatrix<T,MC,MR>& D,
  T beta,        DistMatrix<T,MC,MR>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kTNTN");
    if( E.Height() != E.Width()  || A.Height() != C.Height() ||
        A.Width()  != E.Height() || C.Width()  != E.Height() ||
        B.Width()  != E.Width()  || D.Width()  != E.Width()  ||
        A.Height() != B.Height() || C.Height() != D.Height() )
        throw std::logic_error("Nonconformal Trr2kTNTN");
#endif
    const Grid& g = E.Grid();

    DistMatrix<T,MC,MR> AT(g),  A0(g),
                        AB(g),  A1(g),
                                A2(g);
    DistMatrix<T,MC,MR> BT(g),  B0(g),
                        BB(g),  B1(g),
                                B2(g);

    DistMatrix<T,MC,MR> CT(g),  C0(g),
                        CB(g),  C1(g),
                                C2(g);
    DistMatrix<T,MC,MR> DT(g),  D0(g),
                        DB(g),  D1(g),
                                D2(g);

    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( E );
    B1Trans_MR_STAR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1Trans_MR_STAR.AlignWith( E );

    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    LockedPartitionDown
    ( C, CT,
         CB, 0 );
    LockedPartitionDown
    ( D, DT,
         DB, 0 );
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
        LockedRepartitionDown
        ( DT,  D0,
         /**/ /**/
               D1,
          DB,  D2 );

        //--------------------------------------------------------------------//
        A1_STAR_MC = A1;
        C1_STAR_MC = C1;
        B1Trans_MR_STAR.TransposeFrom( B1 );
        D1Trans_MR_STAR.TransposeFrom( D1 );
        LocalTrr2k
        ( uplo, orientationOfA, TRANSPOSE, orientationOfC, TRANSPOSE,
          alpha, A1_STAR_MC, B1Trans_MR_STAR, 
                 C1_STAR_MC, D1Trans_MR_STAR,
          beta,  E );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( DT,  D0,
               D1,
         /**/ /**/
          DB,  D2 );
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
