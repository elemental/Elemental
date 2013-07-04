/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_LLVF_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_LLVF_HPP

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
// Since applying Householder transforms from vectors stored left-to-right
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
LLVF
( Conjugation conjugation, int offset, 
  const Matrix<F>& H, const Matrix<F>& t, Matrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("apply_packed_reflectors::LLVF");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Height() != A.Height() )
        throw std::logic_error
        ("Height of transforms must equal height of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag.");
#endif
    Matrix<F>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<F>
        AT,  A0,
        AB,  A1,
             A2;
    Matrix<F>
        tT,  t0,
        tB,  t1,
             t2;

    Matrix<F> SInv, Z;

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
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

        int HPanHeight = H11.Height() + H21.Height();
        int HPanWidth = std::min( H11.Width(), std::max(HPanHeight+offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2, HPanWidth );

        RepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LOWER, HPanCopy, offset );
        SetDiagonal( HPanCopy, F(1), offset );

        Herk( LOWER, ADJOINT, F(1), HPanCopy, SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( ADJOINT, NORMAL, F(1), HPanCopy, AB, Z );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), SInv, Z );
        Gemm( NORMAL, NORMAL, F(-1), HPanCopy, Z, F(1), AB );
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
LLVF
( Conjugation conjugation, int offset, 
  const DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t, DistMatrix<F>& A )
{
#ifndef RELEASE
    CallStackEntry entry("apply_packed_reflectors::LLVF");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        throw std::logic_error
        ("{H,t,A} must be distributed over the same grid");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Height() != A.Height() )
        throw std::logic_error
        ("Height of transforms must equal height of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag.");
    if( !t.AlignedWithDiagonal( H, offset ) )
        throw std::logic_error("t must be aligned with H's 'offset' diagonal");
#endif
    const Grid& g = H.Grid();

    DistMatrix<F>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<F>
        AT(g),  A0(g),
        AB(g),  A1(g),
                A2(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<F,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<F,MC,  STAR> HPan_MC_STAR(g);
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > Z_STAR_MR(g);
    DistMatrix<F,STAR,VR  > Z_STAR_VR(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
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

        int HPanHeight = H11.Height() + H21.Height();
        int HPanWidth = std::min( H11.Width(), std::max(HPanHeight+offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2, HPanWidth );

        RepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        HPan_MC_STAR.AlignWith( AB );
        Z_STAR_MR.AlignWith( AB );
        Z_STAR_VR.AlignWith( AB );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LOWER, HPanCopy, offset );
        SetDiagonal( HPanCopy, F(1), offset );

        HPan_VC_STAR = HPanCopy;
        Zeros( SInv_STAR_STAR, HPan.Width(), HPan.Width() );
        Herk
        ( LOWER, ADJOINT, 
          F(1), HPan_VC_STAR.LockedMatrix(),
          F(0), SInv_STAR_STAR.Matrix() );     
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_MC_STAR = HPanCopy;
        LocalGemm( ADJOINT, NORMAL, F(1), HPan_MC_STAR, AB, Z_STAR_MR );
        Z_STAR_VR.SumScatterFrom( Z_STAR_MR );
        
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, F(1), SInv_STAR_STAR, Z_STAR_VR );

        Z_STAR_MR = Z_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), HPan_MC_STAR, Z_STAR_MR, F(1), AB );
        //--------------------------------------------------------------------//
        HPan_MC_STAR.FreeAlignments();
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

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_LLVF_HPP
