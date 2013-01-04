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
SyrkLN( T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C )
{
#ifndef RELEASE
    PushCallStack("internal::SyrkLN");
    if( A.Grid() != C.Grid() )
        throw std::logic_error
        ("A and C must be distributed over the same grid");
    if( A.Height() != C.Height() || A.Height() != C.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal SyrkLN:\n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
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
    ScaleTrapezoid( beta, LEFT, LOWER, 0, C );
    LockedPartitionRight( A, AL, AR, 0 );
    while( AR.Width() > 0 )
    {
        LockedRepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        //--------------------------------------------------------------------//
        A1_VR_STAR = A1_MC_STAR = A1;
        A1Trans_STAR_MR.TransposeFrom( A1_VR_STAR );
        LocalTrrk( LOWER, alpha, A1_MC_STAR, A1Trans_STAR_MR, T(1), C );
        //--------------------------------------------------------------------//

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
