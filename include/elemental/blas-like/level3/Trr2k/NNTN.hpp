/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRR2K_NNTN_HPP
#define BLAS_TRR2K_NNTN_HPP

namespace elem {
namespace internal {

// Distributed E := alpha (A B + C^{T/H} D) + beta E
template<typename T>
inline void
Trr2kNNTN
( UpperOrLower uplo,
  Orientation orientationOfC,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kNNTN");
    if( E.Height() != E.Width()  || A.Width()  != C.Height() ||
        A.Height() != E.Height() || C.Width()  != E.Height() ||
        B.Width()  != E.Width()  || D.Width()  != E.Width()  ||
        A.Width()  != B.Height() || C.Height() != D.Height() )
        throw std::logic_error("Nonconformal Trr2kNNTN");
#endif
    const Grid& g = E.Grid();

    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    DistMatrix<T> CT(g),  C0(g),
                  CB(g),  C1(g),
                          C2(g);
    DistMatrix<T> DT(g),  D0(g),
                  DB(g),  D1(g),
                          D2(g);

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( E );
    B1Trans_MR_STAR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1Trans_MR_STAR.AlignWith( E );

    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
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
        A1_MC_STAR = A1;
        C1_STAR_MC = C1;
        B1Trans_MR_STAR.TransposeFrom( B1 );
        D1Trans_MR_STAR.TransposeFrom( D1 );
        LocalTrr2k
        ( uplo, TRANSPOSE, orientationOfC, TRANSPOSE,
          alpha, A1_MC_STAR, B1Trans_MR_STAR, 
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

#endif // ifndef BLAS_TRR2K_NNTN_HPP
