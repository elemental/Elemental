/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

template<typename T>
inline void
Syr2kLN
( T alpha, const DistMatrix<T>& A, const DistMatrix<T>& B,
  T beta,        DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::Syr2kLN");
    if( A.Grid() != B.Grid() || B.Grid() != C.Grid() )
        throw std::logic_error
        ("{A,B,C} must be distributed over the same grid");
    if( A.Height() != C.Height() || A.Height() != C.Width() ||
        B.Height() != C.Height() || B.Height() != C.Width() ||
        A.Width() != B.Width()                                 )
    {
        std::ostringstream msg;
        msg << "Nonconformal Syr2kLN:\n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  B ~ " << B.Height() << " x " << B.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();

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
    ScaleTrapezoid( beta, LEFT, LOWER, 0, C );
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
        A1Trans_STAR_MR.TransposeFrom( A1_VR_STAR );

        B1_VR_STAR = B1_MC_STAR = B1;
        B1Trans_STAR_MR.TransposeFrom( B1_VR_STAR );

        LocalTrr2k
        ( LOWER, 
          alpha, A1_MC_STAR, B1Trans_STAR_MR,
                 B1_MC_STAR, A1Trans_STAR_MR,
          T(1),  C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionRight
        ( BL,     /**/ BR,
          B0, B1, /**/ B2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
