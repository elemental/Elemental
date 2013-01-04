/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Distributed C := alpha A^{T/H} B + beta C
template<typename T>
inline void
TrrkTN
( UpperOrLower uplo,
  Orientation orientationOfA,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::TrrkTN");
    if( C.Height() != C.Width() ||
        A.Width() != C.Height() || 
        B.Width() != C.Width() ||
        A.Height() != B.Height() )
        throw std::logic_error("Nonconformal TrrkTN");
    if( orientationOfA == NORMAL )
        throw std::logic_error("Orientation must be ADJOINT or NORMAL");
#endif
    const Grid& g = C.Grid();

    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    LockedPartitionDown
    ( A, AT,
         AB, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
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

        //--------------------------------------------------------------------//
        A1_STAR_MC = A1;
        B1Trans_MR_STAR.TransposeFrom( B1 );
        LocalTrrk
        ( uplo, orientationOfA, TRANSPOSE, 
          alpha, A1_STAR_MC, B1Trans_MR_STAR, beta, C );
        //--------------------------------------------------------------------//

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
