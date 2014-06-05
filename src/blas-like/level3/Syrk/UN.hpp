/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace syrk {

template<typename T>
inline void
UN
( T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C, 
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("syrk::UN");
        if( A.Grid() != C.Grid() )
            LogicError("A and C must be distributed over the same grid");
        if( A.Height() != C.Height() || A.Height() != C.Width() )
            LogicError
            ("Nonconformal:\n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  C ~ ",C.Height()," x ",C.Width(),"\n");
    )
    const Grid& g = A.Grid();

    // Matrix views
    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);

    // Temporary distributions
    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> A1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > A1Trans_STAR_MR(g);

    A1_MC_STAR.AlignWith( C );
    A1_VR_STAR.AlignWith( C );
    A1Trans_STAR_MR.AlignWith( C );

    // Start the algorithm
    ScaleTrapezoid( beta, UPPER, C );
    LockedPartitionRight( A, AL, AR, 0 );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        //--------------------------------------------------------------------//
        A1_VR_STAR = A1_MC_STAR = A1;
        A1_VR_STAR.TransposePartialColAllGather( A1Trans_STAR_MR, conjugate );
        LocalTrrk( UPPER, alpha, A1_MC_STAR, A1Trans_STAR_MR, T(1), C ); 
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
}

} // namespace syrk
} // namespace El
