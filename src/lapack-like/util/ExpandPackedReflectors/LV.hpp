/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#pragma once
#ifndef EL_EXPANDPACKEDREFLECTORS_LV_HPP
#define EL_EXPANDPACKEDREFLECTORS_LV_HPP

#include EL_IDENTITY_INC


namespace El {
namespace expand_packed_reflectors {

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

template<typename F>
inline void
LV( Conjugation conjugation, Int offset, Matrix<F>& H, const Matrix<F>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("expand_packed_reflectors::LV");
        if( offset > 0 || offset < -H.Height() )
            LogicError("Transforms out of bounds");
        if( t.Height() != H.DiagonalLength( offset ) )
            LogicError("t must be the same length as H's offset diag");
    )
    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    const Int m = H.Height();
    const Int n = Min(m,H.Width());
    const Int diff = m-n;
    H.Resize( m, n );
    MakeTrapezoidal( LOWER, H, offset );
    SetDiagonal( H, F(1), offset );

    Matrix<F>
        HTL, HTR,  H00, H01, H02,  HPan, HPanCopy, HPanT,
        HBL, HBR,  H10, H11, H12,                  HPanB,
                   H20, H21, H22;
    Matrix<F> HEffected, 
              HEffectedNew, HEffectedOld, 
              HEffectedOldT,
              HEffectedOldB;
    Matrix<F>
        tT,  t0,
        tB,  t1,
             t2;

    Matrix<F> SInv, Z, ZNew, ZOld;

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    Int oldEffectedHeight=diff;

    while( HBR.Height() < m && HBR.Width() < n )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const Int HPanHeight = H11.Height() + H21.Height();
        const Int effectedHeight = Max( HPanHeight+offset, 0 );
        const Int HPanWidth = Min( H11.Width(), effectedHeight );

        const Int effectedWidth = effectedHeight - diff;
        const Int oldEffectedWidth = oldEffectedHeight - diff;
        const Int newEffectedWidth = effectedWidth - oldEffectedWidth;

        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );
        LockedPartitionDown
        ( HPan, HPanT,
                HPanB, newEffectedWidth /* to match ZNew */ );

        View
        ( HEffected, H, m-effectedHeight, n-effectedWidth, 
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

        Z.Resize( HPanWidth, effectedWidth );
        PartitionLeft( Z, ZNew, ZOld, oldEffectedWidth );
        //--------------------------------------------------------------------//
        Herk( UPPER, ADJOINT, F(1), HPan, SInv );
        FixDiagonal( conjugation, t1, SInv );

        // Interleave the updates of the already effected portion of the matrix
        // with the newly effected portion to increase performance
        Adjoint( HPanT, ZNew );
        Zero( ZOld );
        Gemm( ADJOINT, NORMAL, F(1), HPanB, HEffectedOldB, F(0), ZOld );
        Trsm( LEFT, UPPER, NORMAL, NON_UNIT, F(1), SInv, Z );
        HPanCopy = HPan;
        MakeIdentity( HEffectedNew );
        Gemm( NORMAL, NORMAL, F(-1), HPanCopy, Z, F(1), HEffected );
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
    const Int oldEffectedWidth = oldEffectedHeight - diff;
    if( oldEffectedWidth < n )
    {
        View( HEffectedNew, H, 0, 0, m, n-oldEffectedWidth );
        Zero( HEffectedNew );
        SetDiagonal( HEffectedNew, F(1) );
    }
}

template<typename F>
inline void
LV
( Conjugation conjugation, Int offset, 
  DistMatrix<F>& H, const DistMatrix<F,MD,STAR>& t )
{
    DEBUG_ONLY(
        CallStackEntry cse("expand_packed_reflectors::LV");
        if( H.Grid() != t.Grid() )
            LogicError("H and t must be distributed over same grid");
        if( offset > 0 || offset < -H.Height() )
            LogicError("Transforms out of bounds");
        if( t.Height() != H.DiagonalLength( offset ) )
            LogicError("t must be the same length as H's offset diag");
        if( !H.DiagonalAlignedWith( t, offset ) )
            LogicError("t must be aligned with H's 'offset' diagonal");
    )
    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    const Int m = H.Height();
    const Int n = Min(m,H.Width());
    const Int diff = m-n;
    H.Resize( m, n );
    MakeTrapezoidal( LOWER, H, offset );
    SetDiagonal( H, F(1), offset );

    const Grid& g = H.Grid();
    DistMatrix<F>
        HTL(g), HTR(g),  H00(g), H01(g), H02(g),  HPan(g),
        HBL(g), HBR(g),  H10(g), H11(g), H12(g),  
                         H20(g), H21(g), H22(g);
    DistMatrix<F> HEffected(g),
                  HEffectedNew(g), HEffectedOld(g),
                  HEffectedOldT(g),
                  HEffectedOldB(g);
    DistMatrix<F,MD,STAR>
        tT(g),  t0(g),
        tB(g),  t1(g),
                t2(g);

    DistMatrix<F,VC,STAR> HPan_VC_STAR(g);
    DistMatrix<F,MC,STAR> HPan_MC_STAR(g), HPanT_MC_STAR(g),
                                           HPanB_MC_STAR(g);

    DistMatrix<F,STAR,MR> Z_STAR_MR(g),
                          ZNew_STAR_MR(g), ZOld_STAR_MR(g);
    DistMatrix<F,STAR,VR> Z_STAR_VR(g),
                          ZNew_STAR_VR(g), ZOld_STAR_VR(g);
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);

    LockedPartitionUpDiagonal
    ( H, HTL, HTR,
         HBL, HBR, 0 );
    LockedPartitionUp
    ( t, tT,
         tB, 0 );
    Int oldEffectedHeight=diff;
    while( HBR.Height() < m && HBR.Width() < n )
    {
        LockedRepartitionUpDiagonal
        ( HTL, /**/ HTR,  H00, H01, /**/ H02,
               /**/       H10, H11, /**/ H12,
         /*************/ /******************/
          HBL, /**/ HBR,  H20, H21, /**/ H22 );

        const Int HPanHeight = H11.Height() + H21.Height();
        const Int effectedHeight = Max( HPanHeight+offset, 0 );
        const Int HPanWidth = Min( H11.Width(), effectedHeight );

        const Int effectedWidth = effectedHeight - diff;
        const Int oldEffectedWidth = oldEffectedHeight - diff;
        const Int newEffectedWidth = effectedWidth - oldEffectedWidth;

        LockedView( HPan, H, H00.Height(), H00.Width(), HPanHeight, HPanWidth );

        View
        ( HEffected, H, m-effectedHeight, n-effectedWidth, 
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
        Z_STAR_MR.Resize( HPanWidth, effectedWidth );
        Z_STAR_VR.Resize( HPanWidth, effectedWidth );
        PartitionLeft
        ( Z_STAR_MR, ZNew_STAR_MR, ZOld_STAR_MR, oldEffectedWidth );
        PartitionLeft
        ( Z_STAR_VR, ZNew_STAR_VR, ZOld_STAR_VR, oldEffectedWidth );
        //--------------------------------------------------------------------//
        HPan_VC_STAR = HPan;
        Zeros( SInv_STAR_STAR, HPanWidth, HPanWidth );
        Herk
        ( UPPER, ADJOINT, 
          F(1), HPan_VC_STAR.LockedMatrix(), 
          F(0), SInv_STAR_STAR.Matrix() );
        SInv_STAR_STAR.SumOver( HPan_VC_STAR.ColComm() );
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
        Zero( ZOld_STAR_MR );
        LocalGemm
        ( ADJOINT, NORMAL, 
          F(1), HPanB_MC_STAR, HEffectedOldB, F(0), ZOld_STAR_MR );
        ZOld_STAR_VR.PartialRowSumScatterFrom( ZOld_STAR_MR );
        LocalTrsm
        ( LEFT, UPPER, NORMAL, NON_UNIT, F(1), SInv_STAR_STAR, Z_STAR_VR );
        Z_STAR_MR = Z_STAR_VR;
        MakeIdentity( HEffectedNew );
        LocalGemm
        ( NORMAL, NORMAL, F(-1), HPan_MC_STAR, Z_STAR_MR, F(1), HEffected );
        //--------------------------------------------------------------------//

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
    const Int oldEffectedWidth = oldEffectedHeight - diff;
    if( oldEffectedWidth < n )
    {
        View( HEffectedNew, H, 0, 0, m, n-oldEffectedWidth );
        Zero( HEffectedNew );
        SetDiagonal( HEffectedNew, F(1) );
    }
}

} // namespace expand_packed_reflectors
} // namespace El

#endif // ifndef EL_EXPANDPACKEDREFLECTORS_LV_HPP
