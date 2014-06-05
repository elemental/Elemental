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
LT
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("syr2k::LT");
        if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
        if( A.Width() != C.Height() || A.Width() != C.Width()  ||
            B.Width() != C.Height() || B.Width() != C.Width()  ||
            A.Height() != B.Height() )
            LogicError
            ("Nonconformal:\n",
             DimsString(A,"A"),"\n",DimsString(B,"B"),"\n",DimsString(C,"C"));
    )
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Matrix views
    DistMatrix<T> AT(g),  A0(g),
                  AB(g),  A1(g),
                          A2(g);
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    // Temporary distributions
    DistMatrix<T,MR,  STAR> A1Trans_MR_STAR(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,STAR,VR  > A1_STAR_VR(g);
    DistMatrix<T,STAR,VR  > B1_STAR_VR(g);
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);

    A1Trans_MR_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );
    A1_STAR_MC.AlignWith( C );
    B1_STAR_MC.AlignWith( C );

    // Start the algorithm
    ScaleTrapezoid( beta, LOWER, C );
    LockedPartitionDown
    ( A, AT, 
         AB, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
    while( AB.Height() > 0 )
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
        A1.TransposeColAllGather( A1Trans_MR_STAR );
        A1_STAR_VR.TransposePartialRowFilterFrom( A1Trans_MR_STAR );
        A1_STAR_MC = A1_STAR_VR;

        B1.TransposeColAllGather( B1Trans_MR_STAR );
        B1_STAR_VR.TransposePartialRowFilterFrom( B1Trans_MR_STAR );
        B1_STAR_MC = B1_STAR_VR;

        LocalTrr2k
        ( LOWER, orientation, TRANSPOSE, orientation, TRANSPOSE, 
          alpha, A1_STAR_MC, B1Trans_MR_STAR,
                 B1_STAR_MC, A1Trans_MR_STAR, T(1), C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );
    }
}

} // namespace syr2k
} // namespace El
