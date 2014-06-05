/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRRK_NN_HPP
#define EL_TRRK_NN_HPP

namespace El {
namespace trrk {

// Distributed C := alpha A B + beta C
template<typename T>
void TrrkNN
( UpperOrLower uplo,
  T alpha, const DistMatrix<T>& A,
           const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
    DEBUG_ONLY(
        CallStackEntry cse("trrk::TrrkNN");
        if( C.Height() != C.Width() || A.Height() != C.Height() || 
            B.Width() != C.Width() || A.Width() != B.Height() )
            LogicError("Nonconformal TrrkNN");
    )
    const Grid& g = C.Grid();

    DistMatrix<T> AL(g), AR(g),
                  A0(g), A1(g), A2(g);
    DistMatrix<T> BT(g),  B0(g),
                  BB(g),  B1(g),
                          B2(g);

    DistMatrix<T,MC,STAR> A1_MC_STAR(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    LockedPartitionRight( A, AL, AR, 0 );
    LockedPartitionDown
    ( B, BT,
         BB, 0 );
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

        //--------------------------------------------------------------------//
        A1_MC_STAR = A1;
        B1.TransposeColAllGather( B1Trans_MR_STAR );
        LocalTrrk
        ( uplo, TRANSPOSE, alpha, A1_MC_STAR, B1Trans_MR_STAR, beta, C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( BT,  B0,
               B1,
         /**/ /**/
          BB,  B2 );
        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
}

} // namespace trrk
} // namespace El

#endif // ifndef EL_TRRK_NN_HPP
