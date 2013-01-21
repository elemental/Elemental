/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef BLAS_SYRK_UT_HPP
#define BLAS_SYRK_UT_HPP

namespace elem {
namespace internal {

template<typename T>
inline void
SyrkUT
( T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C, bool conjugate )
{
#ifndef RELEASE
    PushCallStack("internal::SyrkUT");
    if( A.Grid() != C.Grid() )
        throw std::logic_error
        ("A and C must be distributed over the same grid");
    if( A.Width() != C.Height() || A.Width() != C.Width() )
    {
        std::ostringstream msg;
        msg << "Nonconformal SyrkUT:\n"
            << "  A ~ " << A.Height() << " x " << A.Width() << "\n"
            << "  C ~ " << C.Height() << " x " << C.Width() << "\n";
        throw std::logic_error( msg.str().c_str() );
    }
#endif
    const Grid& g = A.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Matrix views
    DistMatrix<T> AT(g),  A0(g), 
                  AB(g),  A1(g),
                          A2(g);

    // Temporary distributions
    DistMatrix<T,MR,  STAR> A1Trans_MR_STAR(g);
    DistMatrix<T,STAR,VR  > A1_STAR_VR(g);
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);

    A1Trans_MR_STAR.AlignWith( C );
    A1_STAR_MC.AlignWith( C );

    // Start the algorithm
    ScaleTrapezoid( beta, LEFT, UPPER, 0, C );
    LockedPartitionUp
    ( A, AT, 
         AB, 0 );
    while( AT.Height() > 0 )
    {
        LockedRepartitionUp
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );

        //--------------------------------------------------------------------//
        A1Trans_MR_STAR.TransposeFrom( A1 );
        A1_STAR_VR.TransposeFrom( A1Trans_MR_STAR );
        A1_STAR_MC = A1_STAR_VR;

        LocalTrrk
        ( UPPER, orientation, TRANSPOSE, 
          alpha, A1_STAR_MC, A1Trans_MR_STAR, T(1), C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionUp
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem

#endif // ifndef BLAS_SYRK_UT_HPP
