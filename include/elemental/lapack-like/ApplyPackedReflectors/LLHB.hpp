/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_LLHB_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_LLHB_HPP

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
// Since applying Householder transforms from vectors stored bottom-to-top
// implies that we will be forming a generalization of 
//
//   (I - tau_0 v_0^H v_0) (I - tau_1 v_1^H v_1) = 
//   I - tau_0 v_0^H v_0 - tau_1 v_1^H v_1 + (tau_0 tau_1 v_0 v_1^H) v_0^H v_1 =
//   I - [ v_0^H, v_1^H ] [  tau_0, -tau_0 tau_1 v_0 v_1^H ] [ v_0 ]
//                        [  0,      tau_1                 ] [ v_1 ],
//
// which has a upper-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   triu(T) = triu( V V^H ),  diag(T) = 1/t or 1/conj(t),
//
// where V is the matrix of Householder vectors and t is the vector of scalars.
// V is stored row-wise in the matrix.
//

template<typename F> 
inline void
LLHB
( Conjugation conjugation, int offset, 
  const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("apply_packed_reflectors::LLHB");
    if( offset > 0 || offset < -H.Width() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Height() )
        throw std::logic_error
        ("Width of transforms must equal height of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
#endif
    Matrix<F>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<F> ATop;
    Matrix<F>
        tT,  t0,
        tB,  t1,
             t2;

    Matrix<F> SInv, Z;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
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

        const int HPanWidth = H10.Width() + H11.Width();
        const int HPanOffset = 
            std::min( H11.Height(), std::max(-offset-H00.Height(),0) );
        const int HPanHeight = H11.Height()-HPanOffset;
        LockedView
        ( HPan, H, H00.Height()+HPanOffset, 0, HPanHeight, HPanWidth );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/ 
          tB,  t2, HPanHeight );

        View( ATop, A, 0, 0, HPanWidth, A.Width() );

        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LOWER, HPanCopy, offset, RIGHT );
        SetDiagonal( HPanCopy, F(1), offset, RIGHT );

        Herk( UPPER, NORMAL, F(1), HPanCopy, SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, NORMAL, F(1), HPanCopy, ATop, Z );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), SInv, Z );
        Gemm( ADJOINT, NORMAL, F(-1), HPanCopy, Z, F(1), ATop );
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
LLHB
( Conjugation conjugation, int offset, 
  const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("apply_packed_reflectors::LLHB");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        throw std::logic_error
        ("H, t, and A must be distributed over the same grid");
    if( offset > 0 || offset < -H.Width() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Height() )
        throw std::logic_error
        ("Width of transforms must equal height of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
    if( !t.AlignedWithDiagonal( H, offset ) )
        throw std::logic_error("t must be aligned with H's offset diagonal");
#endif
    const Grid& g = H.Grid();

    DistMatrix<F>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<F> ATop(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<F,STAR,VR  > HPan_STAR_VR(g);
    DistMatrix<F,STAR,MC  > HPan_STAR_MC(g); 
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > Z_STAR_MR(g);
    DistMatrix<F,STAR,VR  > Z_STAR_VR(g);

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
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

        const int HPanWidth = H10.Width() + H11.Width();
        const int HPanOffset = 
            std::min( H11.Height(), std::max(-offset-H00.Height(),0) );
        const int HPanHeight = H11.Height()-HPanOffset;
        LockedView
        ( HPan, H, H00.Height()+HPanOffset, 0, HPanHeight, HPanWidth );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/ 
          tB,  t2, HPanHeight );

        View( ATop, A, 0, 0, HPanWidth, A.Width() );

        HPan_STAR_MC.AlignWith( ATop );
        Z_STAR_MR.AlignWith( ATop );
        Z_STAR_VR.AlignWith( ATop );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LOWER, HPanCopy, offset, RIGHT );
        SetDiagonal( HPanCopy, F(1), offset, RIGHT );

        HPan_STAR_VR = HPanCopy;
        Zeros( SInv_STAR_STAR, HPanHeight, HPanHeight );
        Herk
        ( UPPER, NORMAL,
          F(1), HPan_STAR_VR.LockedMatrix(),
          F(0), SInv_STAR_STAR.Matrix() );
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_STAR_MC = HPan_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(1), HPan_STAR_MC, ATop, Z_STAR_MR );
        Z_STAR_VR.SumScatterFrom( Z_STAR_MR );

        LocalTrsm
        ( LEFT, UPPER, NORMAL, NON_UNIT, F(1), SInv_STAR_STAR, Z_STAR_VR );

        Z_STAR_MR = Z_STAR_VR;
        LocalGemm
        ( ADJOINT, NORMAL, F(-1), HPan_STAR_MC, Z_STAR_MR, F(1), ATop );
        //--------------------------------------------------------------------//
        HPan_STAR_MC.FreeAlignments();
        Z_STAR_MR.FreeAlignments();
        Z_STAR_VR.FreeAlignments();

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

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_LLHB_HPP
