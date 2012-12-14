/*
   Copyright (c) 2009-2012, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace elem {
namespace internal {

//
// Since applying Householder transforms from vectors stored right-to-left
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
ExpandPackedReflectorsLV( int offset, Matrix<R>& H )
{
#ifndef RELEASE
    PushCallStack("internal::ExpandPackedReflectorsLV");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
#endif
    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    if( H.Width() > H.Height() )
        H.ResizeTo( H.Height(), H.Height() );
    MakeTrapezoidal( LEFT, LOWER, offset, H );
    SetDiagonalToOne( LEFT, offset, H );
    const int dimDiff = H.Height() - H.Width();

    Matrix<R>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy, HPanT,
        HBL, HBR,  H10, H11, H12,                  HPanB,
                   H20, H21, H22;
    Matrix<R> HEffected,
              HEffectedNew, HEffectedOld, 
              HEffectedOldT,
              HEffectedOldB;

    Matrix<R> SInv, Z, ZNew, ZOld;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    int oldEffectedHeight=dimDiff;
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int effectedHeight = std::max(HPanHeight+offset,0);
        const int HPanWidth = std::min( H11.Width(), effectedHeight );

        const int effectedWidth = effectedHeight - dimDiff;
        const int oldEffectedWidth = oldEffectedHeight - dimDiff;
        const int newEffectedWidth = effectedWidth - oldEffectedWidth;

        HPan.LockedView( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );
        LockedPartitionDown
        ( HPan, HPanT,
                HPanB, newEffectedWidth /* to match ZNew */ );

        HEffected.View
        ( H, H.Height()-effectedHeight, H.Width()-effectedWidth, 
          effectedHeight, effectedWidth ); 
        PartitionLeft
        ( HEffected, HEffectedNew, HEffectedOld, oldEffectedWidth );
        PartitionUp
        ( HEffectedOld, HEffectedOldT,
                        HEffectedOldB, oldEffectedHeight );

        Z.ResizeTo( HPanWidth, effectedWidth );
        PartitionLeft( Z, ZNew, ZOld, oldEffectedWidth );
        MakeZeros( ZOld );
        Zeros( HPanWidth, HPanWidth, SInv );
        //--------------------------------------------------------------------//
        Syrk( UPPER, TRANSPOSE, R(1), HPan, R(0), SInv );
        HalveMainDiagonal( SInv );

        // Interleave the updates of the already effected portion of the matrix
        // with the newly effected portion to increase performance
        Transpose( HPanT, ZNew );
        Gemm( TRANSPOSE, NORMAL, R(1), HPanB, HEffectedOldB, R(0), ZOld );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, R(1), SInv, Z );
        HPanCopy = HPan;
        MakeIdentity( HEffectedNew );
        Gemm( NORMAL, NORMAL, R(-1), HPanCopy, Z, R(1), HEffected );
        //--------------------------------------------------------------------//

        oldEffectedHeight = effectedHeight;

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );
    }

    // Take care of any untouched columns on the left side of H
    const int oldEffectedWidth = oldEffectedHeight - dimDiff;
    if( oldEffectedWidth < H.Width() )
    {
        HEffectedNew.View( H, 0, 0, H.Height(), H.Width()-oldEffectedWidth );
        MakeZeros( HEffectedNew );
        SetDiagonalToOne( LEFT, 0, HEffectedNew );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
ExpandPackedReflectorsLV( int offset, DistMatrix<R>& H )
{
#ifndef RELEASE
    PushCallStack("internal::ExpandPackedReflectorsLV");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
#endif
    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    if( H.Width() > H.Height() )
        H.ResizeTo( H.Height(), H.Height() );
    MakeTrapezoidal( LEFT, LOWER, offset, H );
    SetDiagonalToOne( LEFT, offset, H );
    const int dimDiff = H.Height() - H.Width();

    const Grid& g = H.Grid();
    DistMatrix<R>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),  
                         H20(g), H21(g), H22(g);
    DistMatrix<R> HEffected(g),
                  HEffectedNew(g), HEffectedOld(g),
                  HEffectedOldT(g),
                  HEffectedOldB(g);

    DistMatrix<R,VC,STAR> HPan_VC_STAR(g);
    DistMatrix<R,MC,STAR> HPan_MC_STAR(g), HPanT_MC_STAR(g),
                                           HPanB_MC_STAR(g);

    DistMatrix<R,STAR,MR> Z_STAR_MR(g),
                          ZNew_STAR_MR(g), ZOld_STAR_MR(g);
    DistMatrix<R,STAR,VR> Z_STAR_VR(g),
                          ZNew_STAR_VR(g), ZOld_STAR_VR(g);
    DistMatrix<R,STAR,STAR> SInv_STAR_STAR(g);

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    int oldEffectedHeight=dimDiff;
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int effectedHeight = std::max(HPanHeight+offset,0);
        const int HPanWidth = std::min( H11.Width(), effectedHeight );

        const int effectedWidth = effectedHeight - dimDiff;
        const int oldEffectedWidth = oldEffectedHeight - dimDiff;
        const int newEffectedWidth = effectedWidth - oldEffectedWidth;

        HPan.LockedView( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        HEffected.View
        ( H, H.Height()-effectedHeight, H.Width()-effectedWidth, 
          effectedHeight, effectedWidth ); 
        PartitionLeft
        ( HEffected, HEffectedNew, HEffectedOld, oldEffectedWidth );
        PartitionUp
        ( HEffectedOld, HEffectedOldT,
                        HEffectedOldB, oldEffectedHeight );

        HPan_MC_STAR.AlignWith( HEffected );
        Z_STAR_VR.AlignWith( HEffected );
        Z_STAR_MR.AlignWith( HEffected );
        Z_STAR_MR.ResizeTo( HPanWidth, effectedWidth );
        Z_STAR_VR.ResizeTo( HPanWidth, effectedWidth );
        PartitionLeft
        ( Z_STAR_MR, ZNew_STAR_MR, ZOld_STAR_MR, oldEffectedWidth );
        PartitionLeft
        ( Z_STAR_VR, ZNew_STAR_VR, ZOld_STAR_VR, oldEffectedWidth );
        MakeZeros( ZOld_STAR_MR );
        Zeros( HPanWidth, HPanWidth, SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPan_VC_STAR = HPan;
        Syrk
        ( UPPER, TRANSPOSE, 
          R(1), HPan_VC_STAR.LockedLocalMatrix(), 
          R(0), SInv_STAR_STAR.LocalMatrix() );
        SInv_STAR_STAR.SumOverGrid();
        HalveMainDiagonal( SInv_STAR_STAR );

        HPan_MC_STAR = HPan;
        LockedPartitionDown
        ( HPan_MC_STAR, HPanT_MC_STAR,
                        HPanB_MC_STAR, newEffectedWidth /* to match ZNew */ );

        // Interleave the updates of the already effected portion of the matrix
        // with the newly effected portion in order to lower latency and 
        // increase performance
        Transpose( HPanT_MC_STAR, ZNew_STAR_VR );
        LocalGemm
        ( TRANSPOSE, NORMAL, 
          R(1), HPanB_MC_STAR, HEffectedOldB, R(0), ZOld_STAR_MR );
        ZOld_STAR_VR.SumScatterFrom( ZOld_STAR_MR );
        LocalTrsm
        ( LEFT, UPPER, NORMAL, NON_UNIT, R(1), SInv_STAR_STAR, Z_STAR_VR );
        Z_STAR_MR = Z_STAR_VR;
        MakeIdentity( HEffectedNew );
        LocalGemm
        ( NORMAL, NORMAL, R(-1), HPan_MC_STAR, Z_STAR_MR, R(1), HEffected );
        //--------------------------------------------------------------------//
        HPan_MC_STAR.FreeAlignments();
        Z_STAR_VR.FreeAlignments();
        Z_STAR_MR.FreeAlignments();

        oldEffectedHeight = effectedHeight;

        SlideLockedPartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, /**/ H01, H02,
         /*************/ /******************/
               /**/       H10, /**/ H11, H12,
          HBL, /**/ HBR,  H20, /**/ H21, H22 );
    }

    // Take care of any untouched columns on the left side of H
    const int oldEffectedWidth = oldEffectedHeight - dimDiff;
    if( oldEffectedWidth < H.Width() )
    {
        HEffectedNew.View( H, 0, 0, H.Height(), H.Width()-oldEffectedWidth );
        MakeZeros( HEffectedNew );
        SetDiagonalToOne( LEFT, 0, HEffectedNew );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R>
inline void
ExpandPackedReflectorsLV
( Conjugation conjugation, int offset,
  Matrix<Complex<R> >& H, const Matrix<Complex<R> >& t )
{
#ifndef RELEASE
    PushCallStack("internal::ExpandPackedReflectorsLLVB");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
#endif
    typedef Complex<R> C;

    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    if( H.Width() > H.Height() )
        H.ResizeTo( H.Height(), H.Height() );
    MakeTrapezoidal( LEFT, LOWER, offset, H );
    SetDiagonalToOne( LEFT, offset, H );
    const int dimDiff = H.Height() - H.Width();

    Matrix<C>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy, HPanT,
        HBL, HBR,  H10, H11, H12,                  HPanB,
                   H20, H21, H22;
    Matrix<C> HEffected, 
              HEffectedNew, HEffectedOld, 
              HEffectedOldT,
              HEffectedOldB;
    Matrix<C>
        tT,  t0,
        tB,  t1,
             t2;

    Matrix<C> SInv, Z, ZNew, ZOld;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    int oldEffectedHeight=dimDiff;
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int effectedHeight = std::max(HPanHeight+offset,0);
        const int HPanWidth = std::min( H11.Width(), effectedHeight );

        const int effectedWidth = effectedHeight - dimDiff;
        const int oldEffectedWidth = oldEffectedHeight - dimDiff;
        const int newEffectedWidth = effectedWidth - oldEffectedWidth;

        HPan.LockedView( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );
        LockedPartitionDown
        ( HPan, HPanT,
                HPanB, newEffectedWidth /* to match ZNew */ );

        HEffected.View
        ( H, H.Height()-effectedHeight, H.Width()-effectedWidth, 
          effectedHeight, effectedWidth ); 
        PartitionLeft
        ( HEffected, HEffectedNew, HEffectedOld, oldEffectedWidth );
        PartitionUp
        ( HEffectedOld, HEffectedOldT,
                        HEffectedOldB, oldEffectedHeight );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, HPanWidth );

        Z.ResizeTo( HPanWidth, effectedWidth );
        PartitionLeft( Z, ZNew, ZOld, oldEffectedWidth );
        MakeZeros( ZOld );
        Zeros( HPanWidth, HPanWidth, SInv );
        //--------------------------------------------------------------------//
        Herk( UPPER, ADJOINT, C(1), HPan, C(0), SInv );
        FixDiagonal( conjugation, t1, SInv );

        // Interleave the updates of the already effected portion of the matrix
        // with the newly effected portion to increase performance
        Adjoint( HPanT, ZNew );
        Gemm( ADJOINT, NORMAL, C(1), HPanB, HEffectedOldB, C(0), ZOld );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, C(1), SInv, Z );
        HPanCopy = HPan;
        MakeIdentity( HEffectedNew );
        Gemm( NORMAL, NORMAL, C(-1), HPanCopy, Z, C(1), HEffected );
        //--------------------------------------------------------------------//

        oldEffectedHeight = effectedHeight;

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

    // Take care of any untouched columns on the left side of H
    const int oldEffectedWidth = oldEffectedHeight - dimDiff;
    if( oldEffectedWidth < H.Width() )
    {
        HEffectedNew.View( H, 0, 0, H.Height(), H.Width()-oldEffectedWidth );
        MakeZeros( HEffectedNew );
        SetDiagonalToOne( LEFT, 0, HEffectedNew );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

template<typename R> 
inline void
ExpandPackedReflectorsLV
( Conjugation conjugation, int offset, 
  DistMatrix<Complex<R> >& H, const DistMatrix<Complex<R>,MD,STAR>& t )
{
#ifndef RELEASE
    PushCallStack("internal::ExpandPackedReflectorsLV");
    if( H.Grid() != t.Grid() )
        throw std::logic_error("H and t must be distributed over same grid");
    if( offset > 0 || offset < -H.Height() )
        throw std::logic_error("Transforms out of bounds");
    if( t.Height() != H.DiagonalLength( offset ) )
        throw std::logic_error("t must be the same length as H's offset diag");
    if( !t.AlignedWithDiagonal( H, offset ) )
        throw std::logic_error("t must be aligned with H's 'offset' diagonal");
#endif
    typedef Complex<R> C;

    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    if( H.Width() > H.Height() )
        H.ResizeTo( H.Height(), H.Height() );
    MakeTrapezoidal( LEFT, LOWER, offset, H );
    SetDiagonalToOne( LEFT, offset, H );
    const int dimDiff = H.Height() - H.Width();

    const Grid& g = H.Grid();
    DistMatrix<C>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),  
                         H20(g), H21(g), H22(g);
    DistMatrix<C> HEffected(g),
                  HEffectedNew(g), HEffectedOld(g),
                  HEffectedOldT(g),
                  HEffectedOldB(g);
    DistMatrix<C,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<C,VC,STAR> HPan_VC_STAR(g);
    DistMatrix<C,MC,STAR> HPan_MC_STAR(g), HPanT_MC_STAR(g),
                                           HPanB_MC_STAR(g);

    DistMatrix<C,STAR,MR> Z_STAR_MR(g),
                          ZNew_STAR_MR(g), ZOld_STAR_MR(g);
    DistMatrix<C,STAR,VR> Z_STAR_VR(g),
                          ZNew_STAR_VR(g), ZOld_STAR_VR(g);
    DistMatrix<C,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<C,STAR,STAR> SInv_STAR_STAR(g);

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    int oldEffectedHeight=dimDiff;
    while( HBR.Height() < H.Height() && HBR.Width() < H.Width() )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const int HPanHeight = H11.Height() + H21.Height();
        const int effectedHeight = std::max(HPanHeight+offset,0);
        const int HPanWidth = std::min( H11.Width(), effectedHeight );

        const int effectedWidth = effectedHeight - dimDiff;
        const int oldEffectedWidth = oldEffectedHeight - dimDiff;
        const int newEffectedWidth = effectedWidth - oldEffectedWidth;

        HPan.LockedView( H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        HEffected.View
        ( H, H.Height()-effectedHeight, H.Width()-effectedWidth, 
          effectedHeight, effectedWidth ); 
        PartitionLeft
        ( HEffected, HEffectedNew, HEffectedOld, oldEffectedWidth );
        PartitionUp
        ( HEffectedOld, HEffectedOldT,
                        HEffectedOldB, oldEffectedHeight );

        LockedRepartitionUp
        ( tT,  t0,
               t1,
         /**/ /**/
          tB,  t2, HPanWidth );

        HPan_MC_STAR.AlignWith( HEffected );
        Z_STAR_VR.AlignWith( HEffected );
        Z_STAR_MR.AlignWith( HEffected );
        Z_STAR_MR.ResizeTo( HPanWidth, effectedWidth );
        Z_STAR_VR.ResizeTo( HPanWidth, effectedWidth );
        PartitionLeft
        ( Z_STAR_MR, ZNew_STAR_MR, ZOld_STAR_MR, oldEffectedWidth );
        PartitionLeft
        ( Z_STAR_VR, ZNew_STAR_VR, ZOld_STAR_VR, oldEffectedWidth );
        MakeZeros( ZOld_STAR_MR );
        Zeros( HPanWidth, HPanWidth, SInv_STAR_STAR );
        //--------------------------------------------------------------------//
        HPan_VC_STAR = HPan;
        Herk
        ( UPPER, ADJOINT, 
          C(1), HPan_VC_STAR.LockedLocalMatrix(), 
          C(0), SInv_STAR_STAR.LocalMatrix() );
        SInv_STAR_STAR.SumOverGrid();
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_MC_STAR = HPan;
        LockedPartitionDown
        ( HPan_MC_STAR, HPanT_MC_STAR,
                        HPanB_MC_STAR, newEffectedWidth /* to match ZNew */ );

        // Interleave the updates of the already effected portion of the matrix
        // with the newly effected portion to lower latency and increase 
        // performance
        Adjoint( HPanT_MC_STAR, ZNew_STAR_VR );
        LocalGemm
        ( ADJOINT, NORMAL, 
          C(1), HPanB_MC_STAR, HEffectedOldB, C(0), ZOld_STAR_MR );
        ZOld_STAR_VR.SumScatterFrom( ZOld_STAR_MR );
        LocalTrsm
        ( LEFT, UPPER, NORMAL, NON_UNIT, C(1), SInv_STAR_STAR, Z_STAR_VR );
        Z_STAR_MR = Z_STAR_VR;
        MakeIdentity( HEffectedNew );
        LocalGemm
        ( NORMAL, NORMAL, C(-1), HPan_MC_STAR, Z_STAR_MR, C(1), HEffected );
        //--------------------------------------------------------------------//
        HPan_MC_STAR.FreeAlignments();
        Z_STAR_VR.FreeAlignments();
        Z_STAR_MR.FreeAlignments();

        oldEffectedHeight = effectedHeight;

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

    // Take care of any untouched columns on the left side of H
    const int oldEffectedWidth = oldEffectedHeight - dimDiff;
    if( oldEffectedWidth < H.Width() )
    {
        HEffectedNew.View( H, 0, 0, H.Height(), H.Width()-oldEffectedWidth );
        MakeZeros( HEffectedNew );
        SetDiagonalToOne( LEFT, 0, HEffectedNew );
    }
#ifndef RELEASE
    PopCallStack();
#endif
}

} // namespace internal
} // namespace elem
