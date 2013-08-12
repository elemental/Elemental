/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_APPLYPACKEDREFLECTORS_RUHB_HPP
#define ELEM_LAPACK_APPLYPACKEDREFLECTORS_RUHB_HPP

#include "elemental/blas-like/level1/MakeTriangular.hpp"
#include "elemental/blas-like/level1/SetDiagonal.hpp"
#include "elemental/blas-like/level3/Gemm.hpp"
#include "elemental/blas-like/level3/Herk.hpp"
#include "elemental/blas-like/level3/Syrk.hpp"
#include "elemental/blas-like/level3/Trsm.hpp"
#include "elemental/matrices/Zeros.hpp"

namespace elem {
namespace apply_packed_reflectors {

//
// Since applying Householder transforms from vectors stored bottom-to-top
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
//   tril(T) = triu( V V^H ),  diag(T) = 1/t or 1/conj(t),
//
// where V is the matrix of Householder vectors and t is the vector of scalars.
// V is stored row-wise in the matrix.
//

template<typename F>
inline void
RUHB
( Conjugation conjugation, Int offset, 
  const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("apply_packed_reflectors::RUHB");
    // TODO: Proper dimension checks
    if( t.Height() != H.DiagonalLength(offset) )
        LogicError("t must be the same length as H's offset diag");
#endif
    Matrix<F>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<F>
        AL, AR,         ARight,
        A0, A1, A2;
    Matrix<F>
        tT,  t0,
        tB,  t1,
             t2;
    Matrix<F> SInv, Z;

    LockedPartitionUpOffsetDiagonal
    ( offset,
      H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    PartitionLeft( A, AL, AR, HTR.Width() );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        RepartitionLeft
        ( AL,     /**/ AR,
          A0, A1, /**/ A2, H11.Height() );

        LockedView1x2( HPan, H11, H12 );
        View1x2( ARight, A1, A2 );

        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTriangular( UPPER, HPanCopy );
        SetDiagonal( HPanCopy, F(1) );

        Herk( LOWER, NORMAL, F(1), HPanCopy, SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, ADJOINT, F(1), ARight, HPanCopy, Z );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), SInv, Z );
        Gemm( NORMAL, NORMAL, F(-1), Z, HPanCopy, F(1), ARight );
        //--------------------------------------------------------------------//

        SlidePartitionLeft
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        SlideLockedPartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );
    }
}

template<typename F>
inline void
RUHB
( Conjugation conjugation, Int offset, 
  const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("apply_packed_reflectors::RUHB");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        LogicError("{H,t,A} must be distributed over the same grid");
    // TODO: Proper dimension checks
    if( t.Height() != H.DiagonalLength(offset) )
        LogicError("t must be the same length as H's offset diag");
    if( !t.AlignedWithDiagonal( H, offset ) )
        LogicError("t must be aligned with H's offset diagonal");
#endif
    const Grid& g = H.Grid();
    DistMatrix<F>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<F>
        AL(g), AR(g),         ARight(g),
        A0(g), A1(g), A2(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<F,STAR,VR  > HPan_STAR_VR(g);
    DistMatrix<F,STAR,MR  > HPan_STAR_MR(g);
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > ZAdj_STAR_MC(g);
    DistMatrix<F,STAR,VC  > ZAdj_STAR_VC(g);

    LockedPartitionUpOffsetDiagonal
    ( offset,
      H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    PartitionLeft( A, AL, AR, HTR.Width() );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        RepartitionLeft
        ( AL,     /**/ AR,
          A0, A1, /**/ A2, H11.Height() );

        LockedView1x2( HPan, H11, H12 );
        View1x2( ARight, A1, A2 );

        HPan_STAR_MR.AlignWith( ARight );
        ZAdj_STAR_MC.AlignWith( ARight );
        ZAdj_STAR_VC.AlignWith( ARight );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTriangular( UPPER, HPanCopy );
        SetDiagonal( HPanCopy, F(1) );

        HPan_STAR_VR = HPanCopy;
        Zeros( SInv_STAR_STAR, HPan.Height(), HPan.Height() );
        Herk
        ( LOWER, NORMAL,
          F(1), HPan_STAR_VR.LockedMatrix(),
          F(0), SInv_STAR_STAR.Matrix() );
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_STAR_MR = HPan_STAR_VR;
        LocalGemm( NORMAL, ADJOINT, F(1), HPan_STAR_MR, ARight, ZAdj_STAR_MC );
        ZAdj_STAR_VC.SumScatterFrom( ZAdj_STAR_MC );

        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, NORMAL, F(-1), ZAdj_STAR_MC, HPan_STAR_MR, F(1), ARight );
        //--------------------------------------------------------------------//

        SlidePartitionLeft
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        SlideLockedPartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );
    }
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef ELEM_LAPACK_APPLYPACKEDREFLECTORS_RUHB_HPP
