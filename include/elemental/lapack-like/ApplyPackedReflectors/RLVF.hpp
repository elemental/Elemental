/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_RLVF_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_RLVF_HPP

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
//   (I - tau_0 u_0 u_0^H) (I - tau_1 u_1 u_1^H) = 
//   I - tau_0 u_0 u_0^H - tau_1 u_1 u_1^H + (tau_0 tau_1 u_0^H u_1) u_0 u_1^H =
//   I - [ u_0, u_1 ] [ tau_0, -tau_0 tau_1 u_0^H u_1 ] [ u_0^H ]
//                    [ 0,      tau_1                 ] [ u_1^H ],
//
// which has an upper-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   triu(T) = triu( U^H U ),  diag(T) = 1/t or 1/conj(t),
//
// where U is the matrix of Householder vectors and t is the vector of scalars.
//

template<typename R> 
inline void
RLVF( int offset, const Matrix<R>& H, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RLVF");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Height() != A.Width() )
        throw std::logic_error
        ("Height of transforms must equal width of target matrix");
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
    PartitionRight
    ( A, AL, AR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int HPanWidth = 
            std::min( H11.Width(), std::max(HPanHeight+offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        Zeros( AR.Height(), HPanWidth, Z );
        Zeros( HPanWidth, HPanWidth, SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, LOWER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, R(1) );

        Syrk( UPPER, TRANSPOSE, R(1), HPanCopy, R(0), SInv );
        HalveMainDiagonal( SInv );

        Gemm( NORMAL, NORMAL, R(1), AR, HPanCopy, R(0), Z );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, R(1), SInv, Z );
        Gemm( NORMAL, TRANSPOSE, R(-1), Z, HPanCopy, R(1), AR );
        //--------------------------------------------------------------------//

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
RLVF
( int offset, 
  const DistMatrix<R>& H,
        DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RLVF");
    if( H.Grid() != A.Grid() )
        throw std::logic_error("{H,A} must be distributed over the same grid");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Height() != A.Width() )
        throw std::logic_error
        ("Height of transforms must equal width of target matrix");
#endif
    const Grid& g = H.Grid();

    DistMatrix<R>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g), HPanCopy(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),
                         H20(g), H21(g), H22(g);
    DistMatrix<R>
        AL(g), AR(g),
        A0(g), A1(g), A2(g);

    DistMatrix<R,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<R,MR,  STAR> HPan_MR_STAR(g);
    DistMatrix<R,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<R,STAR,MC  > ZTrans_STAR_MC(g);
    DistMatrix<R,STAR,VC  > ZTrans_STAR_VC(g);

    LockedPartitionDownDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    PartitionRight
    ( A, AL, AR, 0 );
    while( HTL.Height() < H.Height() && HTL.Width() < H.Width() )
    {
        LockedRepartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int HPanWidth = 
            std::min( H11.Width(), std::max(HPanHeight+offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        HPan_MR_STAR.AlignWith( AR );
        ZTrans_STAR_MC.AlignWith( AR );
        ZTrans_STAR_VC.AlignWith( AR );
        Zeros( HPanWidth, AR.Height(), ZTrans_STAR_MC );
        Zeros( HPanWidth, HPanWidth, SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, LOWER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, R(1) );

        HPan_VC_STAR = HPanCopy;
        Syrk
        ( UPPER, TRANSPOSE, 
          R(1), HPan_VC_STAR.LockedMatrix(),
          R(0), SInv_STAR_STAR.Matrix() );     
        SInv_STAR_STAR.SumOverGrid();
        HalveMainDiagonal( SInv_STAR_STAR );

        HPan_MR_STAR = HPan_VC_STAR;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, R(1), HPan_MR_STAR, AR, R(0), ZTrans_STAR_MC );
        ZTrans_STAR_VC.SumScatterFrom( ZTrans_STAR_MC );
        
        LocalTrsm
        ( LEFT, UPPER, TRANSPOSE, NON_UNIT, 
          R(1), SInv_STAR_STAR, ZTrans_STAR_VC );

        ZTrans_STAR_MC = ZTrans_STAR_VC;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, R(-1), ZTrans_STAR_MC, HPan_MR_STAR, R(1), AR );
        //--------------------------------------------------------------------//
        HPan_MR_STAR.FreeAlignments();
        ZTrans_STAR_MC.FreeAlignments();
        ZTrans_STAR_VC.FreeAlignments();

        SlideLockedPartitionDownDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
void
RLVF
( Conjugation conjugation, int offset, 
  const Matrix<Complex<R> >& H,
  const Matrix<Complex<R> >& t,
        Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RLVF");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Height() != A.Width() )
        throw std::logic_error
        ("Height of transforms must equal width of target matrix");
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

        const int HPanHeight = H11.Height() + H21.Height();
        const int HPanWidth = 
            std::min( H11.Width(), std::max(HPanHeight+offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2, HPanWidth );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        Zeros( AR.Height(), HPanWidth, Z );
        Zeros( HPanWidth, HPanWidth, SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, LOWER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, C(1) );

        Herk( UPPER, ADJOINT, C(1), HPanCopy, C(0), SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, NORMAL, C(1), AR, HPanCopy, C(0), Z );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, C(1), SInv, Z );
        Gemm( NORMAL, ADJOINT, C(-1), Z, HPanCopy, C(1), AR );
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

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
void
RLVF
( Conjugation conjugation, int offset, 
  const DistMatrix<Complex<R> >& H,
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RLVF");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        throw std::logic_error
        ("{H,t,A} must be distributed over the same grid");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Height() != A.Width() )
        throw std::logic_error
        ("Height of transforms must equal width of target matrix");
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
        AL(g), AR(g),
        A0(g), A1(g), A2(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<C,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<C,MR,  STAR> HPan_MR_STAR(g);
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

        const int HPanHeight = H11.Height() + H21.Height();
        const int HPanWidth = 
            std::min( H11.Width(), std::max(HPanHeight+offset,0) );
        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        LockedRepartitionDown
        ( tT,  t0,
         /**/ /**/
               t1,
          tB,  t2, HPanWidth );

        RepartitionRight
        ( AL, /**/ AR,
          A0, /**/ A1, A2 );

        HPan_MR_STAR.AlignWith( AR );
        ZAdj_STAR_MC.AlignWith( AR );
        ZAdj_STAR_VC.AlignWith( AR );
        Zeros( HPanWidth, AR.Height(), ZAdj_STAR_MC );
        Zeros( HPanWidth, HPanWidth, SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( LEFT, LOWER, offset, HPanCopy );
        SetDiagonal( LEFT, offset, HPanCopy, C(1) );

        HPan_VC_STAR = HPanCopy;
        Herk
        ( UPPER, ADJOINT, 
          C(1), HPan_VC_STAR.LockedMatrix(),
          C(0), SInv_STAR_STAR.Matrix() );     
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_MR_STAR = HPan_VC_STAR;
        LocalGemm
        ( ADJOINT, ADJOINT, C(1), HPan_MR_STAR, AR, C(0), ZAdj_STAR_MC );
        ZAdj_STAR_VC.SumScatterFrom( ZAdj_STAR_MC );
        
        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, C(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, ADJOINT, C(-1), ZAdj_STAR_MC, HPan_MR_STAR, C(1), AR );
        //--------------------------------------------------------------------//
        HPan_MR_STAR.FreeAlignments();
        ZAdj_STAR_MC.FreeAlignments();
        ZAdj_STAR_VC.FreeAlignments();

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

        SlidePartitionRight
        ( AL,     /**/ AR,
          A0, A1, /**/ A2 );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_RLVF_HPP
