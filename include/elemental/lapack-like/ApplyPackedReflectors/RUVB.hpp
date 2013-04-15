/*
   Copyright (c) 2009-2013, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef LAPACK_APPLYPACKEDREFLECTORS_RUVB_HPP
#define LAPACK_APPLYPACKEDREFLECTORS_RUVB_HPP

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

template<typename R>
inline void
RUVB( int offset, const Matrix<R>& H, Matrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUVB");
    if( offset < 0 || offset > H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Height() != A.Width() )
        throw std::logic_error
        ("Height of transforms must equal width of target matrix");
#endif
    Matrix<R>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy,
        HBL, HBR,  H10, H11, H12,
                   H20, H21, H22;
    Matrix<R> ALeft;

    Matrix<R> SInv, Z;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H01.Height() + H11.Height();
        const int HPanOffset = 
            std::min( H11.Width(), std::max(offset-H00.Width(),0) );
        const int HPanWidth = H11.Width()-HPanOffset;
        LockedView( HPan, H, 0, H00.Width()+HPanOffset, HPanHeight, HPanWidth );

        View( ALeft, A, 0, 0, A.Height(), HPanHeight );

        Zeros( ALeft.Height(), HPan.Width(), Z );
        Zeros( HPan.Width(), HPan.Width(), SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( RIGHT, UPPER, offset, HPanCopy );
        SetDiagonal( RIGHT, offset, HPanCopy, R(1) );

        Syrk( LOWER, TRANSPOSE, R(1), HPanCopy, R(0), SInv );
        HalveMainDiagonal( SInv );

        Gemm( NORMAL, NORMAL, R(1), ALeft, HPanCopy, R(0), Z );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, R(1), SInv, Z );
        Gemm( NORMAL, TRANSPOSE, R(-1), Z, HPanCopy, R(1), ALeft );
        //--------------------------------------------------------------------//

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
RUVB
( int offset, 
  const DistMatrix<R>& H,
        DistMatrix<R>& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUVB");
    if( H.Grid() != A.Grid() )
        throw std::logic_error("{H,A} must be distributed over the same grid");
    if( offset < 0 || offset > H.Height() )
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
    DistMatrix<R> ALeft(g);

    DistMatrix<R,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<R,MR,  STAR> HPan_MR_STAR(g);
    DistMatrix<R,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<R,STAR,MC  > ZTrans_STAR_MC(g);
    DistMatrix<R,STAR,VC  > ZTrans_STAR_VC(g);

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H01.Height() + H11.Height();
        const int HPanOffset = 
            std::min( H11.Width(), std::max(offset-H00.Width(),0) );
        const int HPanWidth = H11.Width()-HPanOffset;
        LockedView( HPan, H, 0, H00.Width()+HPanOffset, HPanHeight, HPanWidth );

        View( ALeft, A, 0, 0, A.Height(), HPanHeight );

        HPan_MR_STAR.AlignWith( ALeft );
        ZTrans_STAR_MC.AlignWith( ALeft );
        ZTrans_STAR_VC.AlignWith( ALeft );
        Zeros( HPan.Width(), ALeft.Height(), ZTrans_STAR_MC );
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

        HPan_MR_STAR = HPan_VC_STAR;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE,
          R(1), HPan_MR_STAR, ALeft, R(0), ZTrans_STAR_MC );
        ZTrans_STAR_VC.SumScatterFrom( ZTrans_STAR_MC );
        
        LocalTrsm
        ( LEFT, LOWER, TRANSPOSE, NON_UNIT, 
          R(1), SInv_STAR_STAR, ZTrans_STAR_VC );

        ZTrans_STAR_MC = ZTrans_STAR_VC;
        LocalGemm
        ( TRANSPOSE, TRANSPOSE, 
          R(-1), ZTrans_STAR_MC, HPan_MR_STAR, R(1), ALeft );
        //--------------------------------------------------------------------//
        HPan_MR_STAR.FreeAlignments();
        ZTrans_STAR_MC.FreeAlignments();
        ZTrans_STAR_VC.FreeAlignments();

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
RUVB
( Conjugation conjugation, int offset, 
  const Matrix<Complex<R> >& H,
  const Matrix<Complex<R> >& t,
        Matrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUVB");
    if( offset < 0 || offset > H.Height() )
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
    Matrix<C> ALeft;
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
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H01.Height() + H11.Height();
        const int HPanOffset = 
            std::min( H11.Width(), std::max(offset-H00.Width(),0) );
        const int HPanWidth = H11.Width()-HPanOffset;
        LockedView( HPan, H, 0, H00.Width()+HPanOffset, HPanHeight, HPanWidth );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, HPanWidth );

        View( ALeft, A, 0, 0, A.Height(), HPanHeight );

        Zeros( ALeft.Height(), HPan.Width(), Z );
        Zeros( HPan.Width(), HPan.Width(), SInv );
        //--------------------------------------------------------------------//
        HPanCopy = HPan;
        MakeTrapezoidal( RIGHT, UPPER, offset, HPanCopy );
        SetDiagonal( RIGHT, offset, HPanCopy, C(1) );

        Herk( LOWER, ADJOINT, C(1), HPanCopy, C(0), SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, NORMAL, C(1), ALeft, HPanCopy, C(0), Z );
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, C(1), SInv, Z );
        Gemm( NORMAL, ADJOINT, C(-1), Z, HPanCopy, C(1), ALeft );
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
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
RUVB
( Conjugation conjugation, int offset, 
  const DistMatrix<Complex<R> >& H,
  const DistMatrix<Complex<R>,MD,STAR>& t,
        DistMatrix<Complex<R> >& A )
{
#ifndef RELEASE
    PushCallStack("apply_packed_reflectors::RUVB");
    if( H.Grid() != t.Grid() || t.Grid() != A.Grid() )
        throw std::logic_error
        ("{H,t,A} must be distributed over the same grid");
    if( offset < 0 || offset > H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( H.Height() != A.Width() )
        throw std::logic_error
        ("Height of transforms must equal width of target matrix");
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
    DistMatrix<C> ALeft(g);
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

        const int HPanHeight = H01.Height() + H11.Height();
        const int HPanOffset = 
            std::min( H11.Width(), std::max(offset-H00.Width(),0) );
        const int HPanWidth = H11.Width()-HPanOffset;
        LockedView( HPan, H, 0, H00.Width()+HPanOffset, HPanHeight, HPanWidth );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, HPanWidth );

        View( ALeft, A, 0, 0, A.Height(), HPanHeight );

        HPan_MR_STAR.AlignWith( ALeft );
        ZAdj_STAR_MC.AlignWith( ALeft );
        ZAdj_STAR_VC.AlignWith( ALeft );
        Zeros( HPan.Width(), ALeft.Height(), ZAdj_STAR_MC );
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

        HPan_MR_STAR = HPan_VC_STAR;
        LocalGemm
        ( ADJOINT, ADJOINT,
          C(1), HPan_MR_STAR, ALeft, C(0), ZAdj_STAR_MC );
        ZAdj_STAR_VC.SumScatterFrom( ZAdj_STAR_MC );
        
        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, C(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, ADJOINT, C(-1), ZAdj_STAR_MC, HPan_MR_STAR, C(1), ALeft );
        //--------------------------------------------------------------------//
        HPan_MR_STAR.FreeAlignments();
        ZAdj_STAR_MC.FreeAlignments();
        ZAdj_STAR_VC.FreeAlignments();

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
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace apply_packed_reflectors
} // namespace elem

#endif // ifndef LAPACK_APPLYPACKEDREFLECTORS_RUVB_HPP
