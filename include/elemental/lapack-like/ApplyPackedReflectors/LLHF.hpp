/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_LLHF_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_LLHF_HPP

#include "elemental/blas-like/level1/MakeTrapezoidal.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Syrk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace apply_packed_reflectors {

//
// Since applying Householder transforms from vectors stored top-to-bottom
// implies that we will be forming a generalization of
//
//   (I - tau_1 v_1^H v_1) (I - tau_0 v_0^H v_0) = 
//   I - tau_0 v_0^H v_0 - tau_1 v_1^H v_1 + (tau_0 tau_1 v_1 v_0^H) v_1^H v_0 =
//   I - [ v_0^H, v_1^H ] [  tau_0,                 0     ] [ v_0 ]
//                        [ -tau_0 tau_1 v_1 v_0^H, tau_1 ] [ v_1 ],
//
// which has a lower-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   tril(T) = tril( V V^H ),  diag(T) = 1/t or 1/conj(t),
//
// where V is the matrix of Householder vectors and t is the vector of scalars.
//

template<typename F> 
inline void
LLHF
( Conjugation conjugation, int offset, 
  const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("apply_packed_reflectors::LLHF");
    // TODO: Proper dimension checks
    if( t.Height() != H.DiagonalLength(offset) )
        throw std::logic_error("t must be the same length as H's offset diag");
#endif
    Matrix<F>
        HTL, HTR,  H00, H01, H02,  HPan,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<F>
        AT,  A0,  ATop,
        AB,  A1,
             A2;
    Matrix<F>
        tT,  t0,
        tB,  t1,
             t2;
    Matrix<F> HPanCopy;
    Matrix<F> SInv, Z;

    LockedPartitionDownOffsetDiagonal
    ( offset,
      H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionDown
    ( t, tT,
         tB, 0 );
    PartitionDown
    ( A, AT,
         AB, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        RepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2, H11.Height() );

        LockedView1x2( HPan, H10, H11 );
        View2x1( ATop, A0, A1 );

        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LOWER, HPanCopy, offset, RIGHT );
        SetDiagonal( HPanCopy, F(1), offset, RIGHT );
        Herk( LOWER, NORMAL, F(1), HPanCopy, SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, NORMAL, F(1), HPanCopy, ATop, Z );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), SInv, Z );
        Gemm( ADJOINT, NORMAL, F(-1), HPanCopy, Z, F(1), ATop );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlideLockedPartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        SlidePartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    }
}

template<typename F> 
inline void
LLHF
( Conjugation conjugation, int offset, 
  const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("apply_packed_reflectors::LLHF");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        throw std::logic_error("{H,t,A} must be distributed over same grid");
    // TODO: Proper dimension checks
    if( t.Height() != H.DiagonalLength(offset) )
        throw std::logic_error("t must be the same length as H's offset diag");
    if( !t.AlignedWithDiagonal( H, offset ) )
        throw std::logic_error("t must be aligned with H's 'offset' diagonal");
#endif
    const Grid& g = H.Grid();
    DistMatrix<F>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<F>
        AT(g),  A0(g),  ATop(g),
        AB(g),  A1(g),
                A2(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<F> HPanCopy(g);
    DistMatrix<F,STAR,VR  > HPan_STAR_VR(g); 
    DistMatrix<F,STAR,MC  > HPan_STAR_MC(g);
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > Z_STAR_MR(g);
    DistMatrix<F,STAR,VR  > Z_STAR_VR(g);

    LockedPartitionDownOffsetDiagonal
    ( offset,
      H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionDown
    ( t, tT,
         tB, 0 );
    PartitionDown
    ( A, AT,
         AB, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        RepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2, H11.Height() );

        LockedView1x2( HPan, H10, H11 );
        View2x1( ATop, A0, A1 );

        HPan_STAR_MC.AlignWith( ATop );
        Z_STAR_MR.AlignWith( ATop );
        Z_STAR_VR.AlignWith( ATop );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LOWER, HPanCopy, offset, RIGHT );
        SetDiagonal( HPanCopy, F(1), offset, RIGHT );
        HPan_STAR_VR = HPanCopy;
        Zeros( SInv_STAR_STAR, HPan.Height(), HPan.Height() );
        Herk
        ( LOWER, NORMAL,
          F(1), HPan_STAR_VR.LockedMatrix(),
          F(0), SInv_STAR_STAR.Matrix() );
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_STAR_MC = HPan_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(1), HPan_STAR_MC, ATop, Z_STAR_MR );
        Z_STAR_VR.SumScatterFrom( Z_STAR_MR );

        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, F(1), SInv_STAR_STAR, Z_STAR_VR );

        Z_STAR_MR = Z_STAR_VR;
        LocalGemm
        ( ADJOINT, NORMAL, F(-1), HPan_STAR_MC, Z_STAR_MR, F(1), ATop );
        //--------------------------------------------------------------------//
        HPan_STAR_MC.FreeAlignments();
        Z_STAR_MR.FreeAlignments();
        Z_STAR_VR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlideLockedPartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        SlidePartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    }
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_LLHF_HPP
