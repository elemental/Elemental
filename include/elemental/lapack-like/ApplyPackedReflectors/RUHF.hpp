/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_RUHF_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_RUHF_HPP

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
//   (I - tau_0 v_0^H v_0) (I - tau_1 v_1^H v_1) = 
//   I - tau_0 v_0^H v_0 - tau_1 v_1^H v_1 + (tau_0 tau_1 v_0 v_1^H) v_0^H v_1 =
//   I - [ v_0^H, v_1^H ] [ tau_0, -tau_0 tau_1 v_0 v_1^H ] [ v_0 ]
//                        [ 0,      tau_1                 ] [ v_1 ],
//
// which has an upper-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   triu(T) = triu( V V^H ),  diag(T) = 1/t or 1/conj(t),
//
// where V is the matrix of Householder vectors and t is the vector of scalars.
// V is stored row-wise in the matrix.
//

template<typename R>
inline void
RUHF( int offset, const Matrix<R>& H, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUHF");
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
        AL, AR,
        A0, A1, A2;

    Matrix<R> SInv, Z;

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    PartitionRight( A, AL, AR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL,  /**/ HTR,  H00, /**/ H01, H02,
         /**************/ /******************/
                /**/       H10, /**/ H11, H12,
          HBL,  /**/ HBR,  H20, /**/ H21, H22 );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        const int HPanWidth = H11.Width() + H12.Width();
        const int HPanHeight = 
            std::min( H11.Height(), std::max(HPanWidth-offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        Zeros( AR.Height(), HPanHeight, Z );
        Zeros( HPanHeight, HPanHeight, SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, UPPER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, R(1) );

        Syrk( UPPER, NORMAL, R(1), HPanCopy, R(0), SInv );
        HalveMainDiagonal( SInv );

        Gemm( NORMAL, TRANSPOSE, R(1), AR, HPanCopy, R(0), Z );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, R(1), SInv, Z );
        Gemm( NORMAL, NORMAL, R(-1), Z, HPanCopy, R(1), AR );
        //--------------------------------------------------------------------//

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
RUHF
( int offset, 
  const DistMatrix<R>& H,
        DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUHF");
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
        AL(g), AR(g),
        A0(g), A1(g), A2(g);

    DistMatrix<R,STAR,VR  > HPan_STAR_VR(g);
    DistMatrix<R,STAR,MR  > HPan_STAR_MR(g);
    DistMatrix<R,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<R,STAR,MC  > ZTrans_STAR_MC(g);
    DistMatrix<R,STAR,VC  > ZTrans_STAR_VC(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    PartitionRight( A, AL, AR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL,  /**/ HTR,  H00, /**/ H01, H02,
         /**************/ /******************/
                /**/       H10, /**/ H11, H12,
          HBL,  /**/ HBR,  H20, /**/ H21, H22 );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        const int HPanWidth = H11.Width() + H12.Width();
        const int HPanHeight = 
            std::min( H11.Height(), std::max(HPanWidth-offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        HPan_STAR_MR.AlignWith( AR );
        ZTrans_STAR_MC.AlignWith( AR );
        ZTrans_STAR_VC.AlignWith( AR );
        Zeros( HPanHeight, AR.Height(), ZTrans_STAR_MC );
        Zeros( HPanHeight, HPanHeight, SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, UPPER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, R(1) );

        HPan_STAR_VR = HPanCopy;
        Syrk
        ( UPPER, NORMAL,
          R(1), HPan_STAR_VR.LockedMatrix(),
          R(0), SInv_STAR_STAR.Matrix() );
        SInv_STAR_STAR.SumOverGrid();
        HalveMainDiagonal( SInv_STAR_STAR );

        HPan_STAR_MR = HPan_STAR_VR;
        LocalGemm
        ( NORMAL, TRANSPOSE, R(1), HPan_STAR_MR, AR, R(0), ZTrans_STAR_MC );
        ZTrans_STAR_VC.SumScatterFrom( ZTrans_STAR_MC );

        LocalTrsm
        ( LEFT, UPPER, TRANSPOSE, NON_UNIT,
          R(1), SInv_STAR_STAR, ZTrans_STAR_VC );

        ZTrans_STAR_MC = ZTrans_STAR_VC;
        LocalGemm
        ( TRANSPOSE, NORMAL, R(-1), ZTrans_STAR_MC, HPan_STAR_MR, R(1), AR );
        //--------------------------------------------------------------------//
        HPan_STAR_MR.FreeAlignments();
        ZTrans_STAR_MC.FreeAlignments();
        ZTrans_STAR_VC.FreeAlignments();

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
RUHF
( Conjugation conjugation, int offset, 
  const Matrix<Complex<R> >& H,
  const Matrix<Complex<R> >& t,
        Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUHF");
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
        AL, AR,
        A0, A1, A2;
    Matrix<C>
        tT,  t0,
        tB,  t1,
             t2;

    Matrix<C> SInv, Z;

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionDown
    ( t, tT,
         tB, 0 );
    PartitionRight
    ( A, AL, AR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        const int HPanWidth = H11.Width() + H12.Width();
        const int HPanHeight = 
            std::min( H11.Height(), std::max(HPanWidth-offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2, HPanHeight );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        Zeros( AR.Height(), HPanHeight, Z );
        Zeros( HPanHeight, HPanHeight, SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, UPPER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, C(1) );

        Herk( UPPER, NORMAL, C(1), HPanCopy, C(0), SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, ADJOINT, C(1), AR, HPanCopy, C(0), Z );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, C(1), SInv, Z );
        Gemm( NORMAL, NORMAL, C(-1), Z, HPanCopy, C(1), AR );
        //--------------------------------------------------------------------//

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
void
RUHF
( Conjugation conjugation, int offset, 
  const DistMatrix<Complex<R> >& H,
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUHF");
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
        throw std::logic_error("t must be aligned with H's 'offset' diagonal");
#endif
    typedef Complex<R> C;
    const Grid& g = H.Grid();

    DistMatrix<C>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<C>
        AL(g), AR(g),
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

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionDown
    ( t, tT,
         tB, 0 );
    PartitionRight
    ( A, AL, AR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        const int HPanWidth = H11.Width() + H12.Width();
        const int HPanHeight = 
            std::min( H11.Height(), std::max(HPanWidth-offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2, HPanHeight );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        HPan_STAR_MR.AlignWith( AR );
        ZAdj_STAR_MC.AlignWith( AR );
        ZAdj_STAR_VC.AlignWith( AR );
        Zeros( HPanHeight, AR.Height(), ZAdj_STAR_MC );
        Zeros( HPanHeight, HPanHeight, SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, UPPER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, C(1) );

        HPan_STAR_VR = HPanCopy;
        Herk
        ( UPPER, NORMAL,
          C(1), HPan_STAR_VR.LockedMatrix(),
          C(0), SInv_STAR_STAR.Matrix() );
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_STAR_MR = HPan_STAR_VR;
        LocalGemm
        ( NORMAL, ADJOINT, C(1), HPan_STAR_MR, AR, C(0), ZAdj_STAR_MC );
        ZAdj_STAR_VC.SumScatterFrom( ZAdj_STAR_MC );

        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, C(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, NORMAL, C(-1), ZAdj_STAR_MC, HPan_STAR_MR, C(1), AR );
        //--------------------------------------------------------------------//
        HPan_STAR_MR.FreeAlignments();
        ZAdj_STAR_MC.FreeAlignments();
        ZAdj_STAR_VC.FreeAlignments();

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );

        SlideLockedPartitionDown
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2 );

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_RUHF_HPP
