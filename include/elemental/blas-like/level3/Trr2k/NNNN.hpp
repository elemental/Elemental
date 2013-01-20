/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_TRR2K_NNNN_HPP
#define BLAS_TRR2K_NNNN_HPP

namespace elem {
namespace internal {

// Distributed E := alpha (A B + C D) + beta E
template<typename T>
inline void
Trr2kNNNN
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
           const DistMatrix<T>& C, const DistMatrix<T>& D,
  T beta,        DistMatrix<T>& E )
{
#ifndef RELEASE
    PushCallStack("internal::Trr2kNNNN");
    if( E.Height() != E.Width()  || A.Width()  != C.Width()  ||
        A.Height() != E.Height() || C.Height() != E.Height() ||
        B.Width()  != E.Width()  || D.Width()  != E.Width()  ||
        A.Width()  != B.Height() || C.Width()  != D.Height() )
        throw std::logic_error("Nonconformal Trr2kNNNN");
#endif
    const Grid& g = E.Grid();

    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    DistMatrix<T> CL(g), CR(g),
                  C0(g), C1(g), C2(g);
    DistMatrix<T> DT(g),  D0(g),
                  DB(g),  D1(g),
                          D2(g);

    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,MC,STAR> C1_MC_STAR(g);
    DistMatrix<T,MR,STAR> D1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( E );
    B1Trans_MR_STAR.AlignWith( E );
    C1_MC_STAR.AlignWith( E );
    D1Trans_MR_STAR.AlignWith( E );

    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    LockedPartitionRight( C, CL, CR, 0 );
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
        LockedRepartitionRight
        ( CL, /**/ CR,
          C0, /**/ C1, C2 );
        LockedRepartitionDown
        ( DT,  D0,
         /**/ /**/
               D1,
          DB,  D2 );

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1;
        C1_MC_STAR = C1;
        B1Trans_MR_STAR.TransposeFrom( B1 );
        D1Trans_MR_STAR.TransposeFrom( D1 );
        LocalTrr2k
        ( uplo, TRANSPOSE, TRANSPOSE, 
          alpha, A1_MC_STAR, B1Trans_MR_STAR, 
                 C1_MC_STAR, D1Trans_MR_STAR,
          beta,  E );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( DT,  D0,
               D1,
         /**/ /**/
          DB,  D2 );
        SlideLockedPartitionRight
        ( CL,     /**/ CR,
          C0, C1, /**/ C2 );
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

#endif // ifndef BLAS_TRR2K_NNNN_HPP
