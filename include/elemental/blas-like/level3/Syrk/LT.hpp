/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_SYRK_LT_HPP
#define ELEM_SYRK_LT_HPP

#include ELEM_SCALETRAPEZOID_INC

namespace elem {
namespace internal {

template<typename T>
inline void
SyrkLT
( T alpha, const DistMatrix<T>& A, T beta, DistMatrix<T>& C, 
  bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("internal::SyrkLT");
        if( A.Grid() != C.Grid() )
            LogicError("A and C must be distributed over the same grid");
        if( A.Width() != C.Height() || A.Width() != C.Width() )
            LogicError
            ("Nonconformal SyrkLT:\n",
             "  A ~ ",A.Height()," x ",A.Width(),"\n",
             "  C ~ ",C.Height()," x ",C.Width());
    )
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
    ScaleTrapezoid( beta, LOWER, C );
    LockedPartitionDown
    ( A, AT, 
         AB, 0 );
    while( AB.Height() > 0 )
    {
        LockedRepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        //--------------------------------------------------------------------//
        A1.TransposeColAllGather( A1Trans_MR_STAR );
        A1_STAR_VR.TransposePartialRowFilterFrom( A1Trans_MR_STAR );
        A1_STAR_MC = A1_STAR_VR;

        LocalTrrk
        ( LOWER, orientation, TRANSPOSE, 
          alpha, A1_STAR_MC, A1Trans_MR_STAR, T(1), C );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    }
}

} // namespace internal
} // namespace elem

#endif // ifndef ELEM_SYRK_LT_HPP
