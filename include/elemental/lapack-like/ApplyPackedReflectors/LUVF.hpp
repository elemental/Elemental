/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_LUVF_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_LUVF_HPP

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

template<typename R> 
inline void
LUVF( int offset, const Matrix<R>& H, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::LUVF");
    if( offset < 0 || offset > H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Height() )
        throw std::logic_error
        ("Width of transforms must equal height of target matrix");
#endif
    Matrix<R>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<R>
        AT,  A0,  ATop,
        AB,  A1,
             A2;

    Matrix<R> SInv, Z;

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
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

        const int HPanHeight = H01.Height() + H11.Height();
        const int HPanOffset = 
            std::min( H11.Width(), std::max(offset-H00.Width(),0) );
        const int HPanWidth = H11.Width()-HPanOffset;
        LockedView( HPan, H, 0, H00.Width()+HPanOffset, HPanHeight, HPanWidth );

        RepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        View2x1( ATop, A0, 
                       A1 );

        Zeros( HPan.Width(), ATop.Width(), Z );
        Zeros( HPan.Width(), HPan.Width(), SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( RIGHT, UPPER, offset, HPanCopy );
        SetDiagonal( RIGHT, offset, HPanCopy, R(1) );
        Syrk( LOWER, TRANSPOSE, R(1), HPanCopy, R(0), SInv );
        HalveMainDiagonal( SInv );

        Gemm( TRANSPOSE, NORMAL, R(1), HPanCopy, ATop, R(0), Z );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, R(1), SInv, Z );
        Gemm( NORMAL, NORMAL, R(-1), HPanCopy, Z, R(1), ATop );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlidePartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
LUVF
( int offset, 
  const DistMatrix<R>& H,
        DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::LUVF");
    if( H.Grid() != A.Grid() )
        throw std::logic_error("{H,A} must be distributed over the same grid");
    if( offset < 0 || offset > H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Height() )
        throw std::logic_error
        ("Width of transforms must equal height of target matrix");
#endif
    const Grid& g = H.Grid();

    DistMatrix<R>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<R>
        AT(g),  A0(g),  ATop(g),
        AB(g),  A1(g),
                A2(g);

    DistMatrix<R> HPanCopy(g);
    DistMatrix<R,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<R,MC,  STAR> HPan_MC_STAR(g);
    DistMatrix<R,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<R,STAR,MR  > Z_STAR_MR(g);
    DistMatrix<R,STAR,VR  > Z_STAR_VR(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
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

        const int HPanHeight = H01.Height() + H11.Height();
        const int HPanOffset = 
            std::min( H11.Width(), std::max(offset-H00.Width(),0) );
        const int HPanWidth = H11.Width()-HPanOffset;
        LockedView( HPan, H, 0, H00.Width()+HPanOffset, HPanHeight, HPanWidth );

        RepartitionDown
        ( AT,  A0,
         /**/ /**/
               A1,
          AB,  A2 );

        View2x1( ATop, A0, 
                       A1 );

        HPan_MC_STAR.AlignWith( ATop );
        Z_STAR_MR.AlignWith( ATop );
        Z_STAR_VR.AlignWith( ATop );
        Zeros( HPan.Width(), ATop.Width(), Z_STAR_MR );
        Zeros( HPan.Width(), HPan.Width(), SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( RIGHT, UPPER, offset, HPanCopy );
        SetDiagonal( RIGHT, offset, HPanCopy, R(1) );
        HPan_VC_STAR = HPanCopy;
        Syrk
        ( LOWER, TRANSPOSE, 
          R(1), HPan_VC_STAR.LockedMatrix(),
          R(0), SInv_STAR_STAR.Matrix() ); 
        SInv_STAR_STAR.SumOverGrid();
        HalveMainDiagonal( SInv_STAR_STAR );

        HPan_MC_STAR = HPanCopy;
        LocalGemm
        ( TRANSPOSE, NORMAL, R(1), HPan_MC_STAR, ATop, R(0), Z_STAR_MR );
        Z_STAR_VR.SumScatterFrom( Z_STAR_MR );
        
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, 
          R(1), SInv_STAR_STAR, Z_STAR_VR );

        Z_STAR_MR = Z_STAR_VR;
        LocalGemm( NORMAL, NORMAL, R(-1), HPan_MC_STAR, Z_STAR_MR, R(1), ATop );
        //--------------------------------------------------------------------//
        HPan_MC_STAR.FreeAlignments();
        Z_STAR_MR.FreeAlignments();
        Z_STAR_VR.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlidePartitionDown
        ( AT,  A0,
               A1,
         /**/ /**/
          AB,  A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
LUVF
( Conjugation conjugation, int offset, 
  const Matrix<Complex<R> >& H,
  const Matrix<Complex<R> >& t,
        Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::LUVF");
    if( offset < 0 || offset > H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Height() )
        throw std::logic_error
        ("Width of transforms must equal height of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
#endif
    typedef Complex<R> C;

    Matrix<C>
        HTL, HTR,  H00, H01, H02,  HPan,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<C>
        AT,  A0,  ATop,
        AB,  A1,
             A2;
    Matrix<C>
        tT,  t0,
        tB,  t1,
             t2;

    Matrix<C> HPanCopy;
    Matrix<C> SInv, Z;

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

        const int HPanHeight = H01.Height() + H11.Height();
        const int HPanOffset = 
            std::min( H11.Width(), std::max(offset-H00.Width(),0) );
        const int HPanWidth = H11.Width()-HPanOffset;
        LockedView( HPan, H, 0, H00.Width()+HPanOffset, HPanHeight, HPanWidth );

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

        View2x1( ATop, A0,
                       A1 );

        Zeros( HPan.Width(), ATop.Width(), Z );
        Zeros( HPan.Width(), HPan.Width(), SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( RIGHT, UPPER, offset, HPanCopy );
        SetDiagonal( RIGHT, offset, HPanCopy, C(1) );
        Herk( LOWER, ADJOINT, C(1), HPanCopy, C(0), SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( ADJOINT, NORMAL, C(1), HPanCopy, ATop, C(0), Z );
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, C(1), SInv, Z );
        Gemm( NORMAL, NORMAL, C(-1), HPanCopy, Z, C(1), ATop );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
LUVF
( Conjugation conjugation, int offset, 
  const DistMatrix<Complex<R> >& H,
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::LUVF");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        throw std::logic_error
        ("{H,t,A} must be distributed over the same grid");
    if( offset < 0 || offset > H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Height() )
        throw std::logic_error
        ("Width of transforms must equal height of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
    if( !t.AlignedWithDiagonal( H, offset ) )
        throw std::logic_error("t must be aligned with H's 'offset' diagonal");
#endif
    typedef Complex<R> C;
    const Grid& g = H.Grid();

    DistMatrix<C>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<C>
        AT(g),  A0(g),  ATop(g),
        AB(g),  A1(g),
                A2(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<C> HPanCopy(g);
    DistMatrix<C,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<C,MC,  STAR> HPan_MC_STAR(g);
    DistMatrix<C,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<C,STAR,MR  > Z_STAR_MR(g);
    DistMatrix<C,STAR,VR  > Z_STAR_VR(g);

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

        const int HPanHeight = H01.Height() + H11.Height();
        const int HPanOffset = 
            std::min( H11.Width(), std::max(offset-H00.Width(),0) );
        const int HPanWidth = H11.Width()-HPanOffset;
        LockedView( HPan, H, 0, H00.Width()+HPanOffset, HPanHeight, HPanWidth );

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

        View2x1( ATop, A0,
                       A1 );

        HPan_MC_STAR.AlignWith( ATop );
        Z_STAR_MR.AlignWith( ATop );
        Z_STAR_VR.AlignWith( ATop );
        Zeros( HPan.Width(), ATop.Width(), Z_STAR_MR );
        Zeros( HPan.Width(), HPan.Width(), SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( RIGHT, UPPER, offset, HPanCopy );
        SetDiagonal( RIGHT, offset, HPanCopy, C(1) );
        HPan_VC_STAR = HPanCopy;
        Herk
        ( LOWER, ADJOINT, 
          C(1), HPan_VC_STAR.LockedMatrix(),
          C(0), SInv_STAR_STAR.Matrix() ); 
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_MC_STAR = HPanCopy;
        LocalGemm( ADJOINT, NORMAL, C(1), HPan_MC_STAR, ATop, C(0), Z_STAR_MR );
        Z_STAR_VR.SumScatterFrom( Z_STAR_MR );
        
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, C(1), SInv_STAR_STAR, Z_STAR_VR );

        Z_STAR_MR = Z_STAR_VR;
        LocalGemm( NORMAL, NORMAL, C(-1), HPan_MC_STAR, Z_STAR_MR, C(1), ATop );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_LUVF_HPP
