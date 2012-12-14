/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

// Distributed C := alpha A B^{T/H} + beta C
template<typename T>
inline void
TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha, const DistMatrix<T,MC,MR>& A,
           const DistMatrix<T,MC,MR>& B,
  T beta,        DistMatrix<T,MC,MR>& C )
{
#ifndef RELEASE
    PushCallStack("internal::TrrkNT");
    if( C.Height() != C.Width() ||
        A.Height() != C.Height() || 
        B.Height() != C.Width() ||
        A.Width() != B.Width() )
        throw std::logic_error("Nonconformal TrrkNT");
    if( orientationOfB == NORMAL )
        throw std::logic_error("Orientation must be ADJOINT or TRANSPOSE");
#endif
    const Grid& g = C.Grid();

    DistMatrix<T,MC,MR> AL(g), AR(g),
                        A0(g), A1(g), A2(g);
    DistMatrix<T,MC,MR> BL(g), BR(g),
                        B0(g), B1(g), B2(g);

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1AdjOrTrans_STAR_MR(g);

    A1_MC_STAR.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1AdjOrTrans_STAR_MR.AlignWith( C );

    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    while( AL.Width() < A.Width() )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );
        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1;
        B1_VR_STAR = B1;
        if( orientationOfB == ADJOINT )
            B1AdjOrTrans_STAR_MR.AdjointFrom( B1_VR_STAR );
        else
            B1AdjOrTrans_STAR_MR.TransposeFrom( B1_VR_STAR );
        LocalTrrk( uplo, alpha, A1_MC_STAR, B1AdjOrTrans_STAR_MR, beta, C );
        //--------------------------------------------------------------------//

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
