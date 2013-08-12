/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef ELEM_LAPACK_APPLYPACKEDREFLECTORS_RUVB_HPP
#define ELEM_LAPACK_APPLYPACKEDREFLECTORS_RUVB_HPP

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
// Since applying Householder transforms from vectors stored right-to-left
// implies that we will be forming a generalization of
//
//   (I - tau_1 u_1 u_1^H) (I - tau_0 u_0 u_0^H) = 
//   I - tau_0 u_0 u_0^H - tau_1 u_1 u_1^H + (tau_0 tau_1 u_1^H u_0) u_1 u_0^H =
//   I - [ u_0, u_1 ] [  tau_0,                 0     ] [ u_0^H ]
//                    [ -tau_0 tau_1 u_1^H u_0, tau_1 ] [ u_1^H ],
//
// which has a lower-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   tril(T) = tril( U^H U ),  diag(T) = 1/t or 1/conj(t),
//
// where U is the matrix of Householder vectors and t is the vector of scalars.
//

template<typename F> 
inline void
RUVB
( Conjugation conjugation, Int offset, 
  const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("apply_packed_reflectors::RUVB");
    // TODO: Proper dimension checks
    if( t.Height() != H.DiagonalLength(offset) )
        LogicError("t must be the same length as H's offset diag");
#endif
    Matrix<F>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<F> ALeft;
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

        LockedView2x1( HPan, H01, H11 );
        View( ALeft, A, 0, 0, A.Height(), HPan.Height() );

        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( UPPER, HPanCopy, 0, RIGHT );
        SetDiagonal( HPanCopy, F(1), 0, RIGHT );

        Herk( LOWER, ADJOINT, F(1), HPanCopy, SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, NORMAL, F(1), ALeft, HPanCopy, Z );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), SInv, Z );
        Gemm( NORMAL, ADJOINT, F(-1), Z, HPanCopy, F(1), ALeft );
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        SlideLockedPartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );
    }
}

template<typename F> 
inline void
RUVB
( Conjugation conjugation, Int offset, 
  const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry cse("apply_packed_reflectors::RUVB");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        LogicError("{H,t,A} must be distributed over the same grid");
    // TODO: Proper dimension checks
    if( t.Height() != H.DiagonalLength(offset) )
        LogicError("t must be the same length as H's offset diag");
    if( !t.AlignedWithDiagonal( H, offset ) )
        LogicError("t must be aligned with H's 'offset' diagonal");
#endif
    const Grid& g = H.Grid();
    DistMatrix<F>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<F> ALeft(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<F,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<F,MR,  STAR> HPan_MR_STAR(g);
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

        LockedView2x1( HPan, H01, H11 );
        View( ALeft, A, 0, 0, A.Height(), HPan.Height() );

        HPan_MR_STAR.AlignWith( ALeft );
        ZAdj_STAR_MC.AlignWith( ALeft );
        ZAdj_STAR_VC.AlignWith( ALeft );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( UPPER, HPanCopy, 0, RIGHT );
        SetDiagonal( HPanCopy, F(1), 0, RIGHT );

        HPan_VC_STAR = HPanCopy;
        Zeros( SInv_STAR_STAR, HPan.Width(), HPan.Width() );
        Herk
        ( LOWER, ADJOINT, 
          F(1), HPan_VC_STAR.LockedMatrix(),
          F(0), SInv_STAR_STAR.Matrix() );     
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_MR_STAR = HPan_VC_STAR;
        LocalGemm( ADJOINT, ADJOINT, F(1), HPan_MR_STAR, ALeft, ZAdj_STAR_MC );
        ZAdj_STAR_VC.SumScatterFrom( ZAdj_STAR_MC );
        
        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, ADJOINT, F(-1), ZAdj_STAR_MC, HPan_MR_STAR, F(1), ALeft );
        //--------------------------------------------------------------------//

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        SlideLockedPartitionUp
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2 );
    }
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef ELEM_LAPACK_APPLYPACKEDREFLECTORS_RUVB_HPP
