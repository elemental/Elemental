/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_RUHB_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_RUHB_HPP

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

template<typename R> 
inline void
RUHB( int offset, const Matrix<R>& H, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUHB");
    if( offset < 0 || offset > H.Width() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Width() )
        throw std::logic_error
        ("Length of transforms must equal width of target matrix");
#endif
    Matrix<R>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy,
        HBL, HBR,  H10, H11, H12, 
                   H20, H21, H22;
    Matrix<R>
        AL, AR,        ARight,
        A0, A1, A2;

    Matrix<R> SInv, Z;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    PartitionLeft( A, AL, AR, std::max(0,H.Width()-H.Height()) );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /**************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanWidth = H11.Width() + H12.Width();
        const int HPanHeight = 
            std::min( H11.Height(), std::max(HPanWidth-offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        RepartitionLeft
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        View1x2( ARight, A1, A2 );

        Zeros( ARight.Height(), HPanHeight, Z );
        Zeros( HPanHeight, HPanHeight, SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, UPPER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, R(1) );

        Syrk( LOWER, NORMAL, R(1), HPanCopy, R(0), SInv );
        HalveMainDiagonal( SInv );

        Gemm( NORMAL, TRANSPOSE, R(1), ARight, HPanCopy, R(0), Z );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, R(1), SInv, Z );
        Gemm( NORMAL, NORMAL, R(-1), Z, HPanCopy, R(1), ARight );
        //--------------------------------------------------------------------//

        SlidePartitionLeft
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
RUHB
( int offset, 
  const DistMatrix<R>& H,
        DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUHB");
    if( H.Grid() != A.Grid() )
        throw std::logic_error("{H,A} must be distributed over the same grid");
    if( offset < 0 || offset > H.Width() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Width() )
        throw std::logic_error
        ("Length of transforms must equal width of target matrix");
#endif
    const Grid& g = H.Grid();

    DistMatrix<R>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g), 
                         H20(g), H21(g), H22(g);
    DistMatrix<R>
        AL(g), AR(g),        ARight(g),
        A0(g), A1(g), A2(g);

    DistMatrix<R,STAR,VR  > HPan_STAR_VR(g);
    DistMatrix<R,STAR,MR  > HPan_STAR_MR(g);
    DistMatrix<R,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<R,STAR,MC  > ZTrans_STAR_MC(g);
    DistMatrix<R,STAR,VC  > ZTrans_STAR_VC(g);

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    PartitionLeft( A, AL, AR, std::max(0,H.Width()-H.Height()) );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /**************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanWidth = H11.Width() + H12.Width();
        const int HPanHeight = 
            std::min( H11.Height(), std::max(HPanWidth-offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        RepartitionLeft
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        View1x2( ARight, A1, A2 );

        HPan_STAR_MR.AlignWith( ARight );
        ZTrans_STAR_MC.AlignWith( ARight );
        ZTrans_STAR_VC.AlignWith( ARight );
        Zeros( HPanHeight, ARight.Height(), ZTrans_STAR_MC );
        Zeros( HPanHeight, HPanHeight, SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, UPPER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, R(1) );

        HPan_STAR_VR = HPanCopy;
        Syrk
        ( LOWER, NORMAL,
          R(1), HPan_STAR_VR.LockedMatrix(),
          R(0), SInv_STAR_STAR.Matrix() );
        SInv_STAR_STAR.SumOverGrid();
        HalveMainDiagonal( SInv_STAR_STAR );

        HPan_STAR_MR = HPan_STAR_VR;
        LocalGemm
        ( NORMAL, TRANSPOSE,
          R(1), HPan_STAR_MR, ARight, R(0), ZTrans_STAR_MC );
        ZTrans_STAR_VC.SumScatterFrom( ZTrans_STAR_MC );

        LocalTrsm
        ( LEFT, LOWER, TRANSPOSE, NON_UNIT,
          R(1), SInv_STAR_STAR, ZTrans_STAR_VC );

        ZTrans_STAR_MC = ZTrans_STAR_VC;
        LocalGemm
        ( TRANSPOSE, NORMAL,
          R(-1), ZTrans_STAR_MC, HPan_STAR_MR, R(1), ARight );
        //--------------------------------------------------------------------//
        ZTrans_STAR_VC.FreeAlignments();
        ZTrans_STAR_MC.FreeAlignments();
        HPan_STAR_MR.FreeAlignments();

        SlidePartitionLeft
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
RUHB
( Conjugation conjugation, int offset, 
  const Matrix<Complex<R> >& H,
  const Matrix<Complex<R> >& t,
        Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUHB");
    if( offset < 0 || offset > H.Width() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Width() )
        throw std::logic_error
        ("Length of transforms must equal width of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
#endif
    typedef Complex<R> C;

    Matrix<C>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<C>
        AL, AR,         ARight,
        A0, A1, A2;
    Matrix<C>
        tT,  t0,
        tB,  t1,
             t2;

    Matrix<C> SInv, Z;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    PartitionLeft
    ( A, AL, AR, std::max(0,H.Width()-H.Height()) );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanWidth = H11.Width() + H12.Width();
        const int HPanHeight = 
            std::min( H11.Height(), std::max(HPanWidth-offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, HPanHeight );

        RepartitionLeft
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        View1x2( ARight, A1, A2 );

        Zeros( ARight.Height(), HPanHeight, Z );
        Zeros( HPanHeight, HPanHeight, SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, UPPER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, C(1) );

        Herk( LOWER, NORMAL, C(1), HPanCopy, C(0), SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, ADJOINT, C(1), ARight, HPanCopy, C(0), Z );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, C(1), SInv, Z );
        Gemm( NORMAL, NORMAL, C(-1), Z, HPanCopy, C(1), ARight );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
RUHB
( Conjugation conjugation, int offset, 
  const DistMatrix<Complex<R> >& H,
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUHB");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        throw std::logic_error
        ("{H,t,A} must be distributed over the same grid");
    if( offset < 0 || offset > H.Width() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Width() != A.Width() )
        throw std::logic_error
        ("Length of transforms must equal width of target matrix");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
    if( !t.AlignedWithDiagonal( H, offset ) )
        throw std::logic_error("t must be aligned with H's offset diagonal");
#endif
    typedef Complex<R> C;
    const Grid& g = H.Grid();

    DistMatrix<C>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<C>
        AL(g), AR(g),         ARight(g),
        A0(g), A1(g), A2(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<C,STAR,VR  > HPan_STAR_VR(g);
    DistMatrix<C,STAR,MR  > HPan_STAR_MR(g);
    DistMatrix<C,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<C,STAR,MC  > ZAdj_STAR_MC(g);
    DistMatrix<C,STAR,VC  > ZAdj_STAR_VC(g);

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    PartitionLeft
    ( A, AL, AR, std::max(0,H.Width()-H.Height()) );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanWidth = H11.Width() + H12.Width();
        const int HPanHeight = 
            std::min( H11.Height(), std::max(HPanWidth-offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, HPanHeight );

        RepartitionLeft
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        View1x2( ARight, A1, A2 );

        HPan_STAR_MR.AlignWith( ARight );
        ZAdj_STAR_MC.AlignWith( ARight );
        ZAdj_STAR_VC.AlignWith( ARight );
        Zeros( HPanHeight, ARight.Height(), ZAdj_STAR_MC );
        Zeros( HPanHeight, HPanHeight, SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, UPPER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, C(1) );

        HPan_STAR_VR = HPanCopy;
        Herk
        ( LOWER, NORMAL,
          C(1), HPan_STAR_VR.LockedMatrix(),
          C(0), SInv_STAR_STAR.Matrix() );
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_STAR_MR = HPan_STAR_VR;
        LocalGemm
        ( NORMAL, ADJOINT, C(1), HPan_STAR_MR, ARight, C(0), ZAdj_STAR_MC );
        ZAdj_STAR_VC.SumScatterFrom( ZAdj_STAR_MC );

        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT,
          C(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, NORMAL, C(-1), ZAdj_STAR_MC, HPan_STAR_MR, C(1), ARight );
        //--------------------------------------------------------------------//
        ZAdj_STAR_VC.FreeAlignments();
        ZAdj_STAR_MC.FreeAlignments();
        HPan_STAR_MR.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_RUHB_HPP
