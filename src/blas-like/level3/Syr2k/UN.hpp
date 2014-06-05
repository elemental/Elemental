/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace syr2k {

template<typename T>
inline void
UN
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("syr2k::UN");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
        if( A.Height() != C.Height() || A.Height() != C.Width() ||
            B.Height() != C.Height() || B.Height() != C.Width() ||
            A.Width() != B.Width() )
            LogicError
            ("Nonconformal:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = C.Grid();

    // Matrix views 
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> BL(g), BR(g),
                  B0(g), B1(g), B2(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,MC,  STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> A1_VR_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > A1Trans_STAR_MR(g);
    DistMatrix<T,STAR,MR  > B1Trans_STAR_MR(g);

    A1_MC_STAR.AlignWith( C );
    B1_MC_STAR.AlignWith( C );
    A1_VR_STAR.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    A1Trans_STAR_MR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );

    // Start the algorithm
    ScaleTrapezoid( beta, UPPER, C );
    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionRight( B, BL, BR, 0 );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        LockedRepartitionRight
        ( BL, /**/ BR,
          B0, /**/ B1, B2 );

        //--------------------------------------------------------------------//
        A1_VR_STAR = A1_MC_STAR = A1;
        A1_VR_STAR.TransposePartialColAllGather( A1Trans_STAR_MR, conjugate );

        B1_VR_STAR = B1_MC_STAR = B1;
        B1_VR_STAR.TransposePartialColAllGather( B1Trans_STAR_MR, conjugate );

        LocalTrr2k
        ( UPPER, 
          alpha, A1_MC_STAR, B1Trans_STAR_MR,
                 B1_MC_STAR, A1Trans_STAR_MR, T(1), C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
    }
}

} // namespace syr2k
} // namespace El
