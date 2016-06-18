/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_EXPANDPACKEDREFLECTORS_LV_HPP
#define EL_EXPANDPACKEDREFLECTORS_LV_HPP

namespace El {
namespace expand_packed_reflectors {

// TODO: Greatly simplify the implementation of these routines

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
void
LV
( Conjugation conjugation,
  Int offset,
        Matrix<F>& H,
  const Matrix<F>& t )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( offset > 0 || offset < -H.Height() )
          LogicError("Transforms out of bounds");
      if( t.Height() != H.DiagonalLength( offset ) )
          LogicError("t must be the same length as H's offset diag");
    )
    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    const Int m = H.Height();
    const Int n = Min(m,H.Width());
    const Int mRem = m-n;
    const Int bsize = Blocksize();
    H.Resize( m, n );
    MakeTrapezoidal( LOWER, H, offset );
    FillDiagonal( H, F(1), offset );

    Matrix<F> HPanT, HPanB,
              HEffectedNew, HEffectedOld,
              HEffectedOldT, HEffectedOldB;

    Matrix<F> HPanCopy, SInv, Z, ZNew, ZOld;

    Int oldEffectedHeight=mRem;

    const Int tOff = ( offset>=0 ? offset : 0 );

    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Int HPanHeight = m-k;
        const Int effectedHeight = Max( HPanHeight+offset, 0 );
        const Int HPanWidth = Min( nb, effectedHeight );

        const Int effectedWidth = effectedHeight - mRem;
        const Int oldEffectedWidth = oldEffectedHeight - mRem;
        const Int newEffectedWidth = effectedWidth - oldEffectedWidth;

        auto HPan = H( IR(k,m), IR(k,k+HPanWidth) );
        LockedPartitionDown
        ( HPan, HPanT,
                HPanB, newEffectedWidth /* to match ZNew */ );

        auto HEffected = H( IR(m-effectedHeight,END), IR(n-effectedWidth,END) );
        PartitionLeft
        ( HEffected, HEffectedNew, HEffectedOld, oldEffectedWidth );
        PartitionUp
        ( HEffectedOld, HEffectedOldT, HEffectedOldB, oldEffectedHeight );

        auto t1 = t( IR(k-tOff,k-tOff+HPanWidth), ALL );

        Z.Resize( HPanWidth, effectedWidth );
        PartitionLeft( Z, ZNew, ZOld, oldEffectedWidth );
        Herk( UPPER, ADJOINT, Base<F>(1), HPan, SInv );
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

        oldEffectedHeight = effectedHeight;
    }

    // Take care of any untouched columns on the left side of H
    const Int oldEffectedWidth = oldEffectedHeight - mRem;
    HEffectedNew = H( ALL, IR(0,n-oldEffectedWidth) );
    MakeIdentity( HEffectedNew );
}

template<typename F>
void
LV
( Conjugation conjugation,
  Int offset, 
        ElementalMatrix<F>& HPre,
  const ElementalMatrix<F>& tPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      AssertSameGrids( HPre, tPre );
      if( offset > 0 || offset < -HPre.Height() )
          LogicError("Transforms out of bounds");
      if( tPre.Height() != HPre.DiagonalLength( offset ) )
          LogicError("t must be the same length as H's offset diag");
    )

    DistMatrixReadWriteProxy<F,F,MC,MR  > HProx( HPre );
    DistMatrixReadProxy<F,F,MC,STAR> tProx( tPre );
    auto& H = HProx.Get();
    auto& t = tProx.GetLocked();

    // Start by zeroing everything above the offset and setting that diagonal
    // to all ones. We can also ensure that H is not wider than it is tall.
    const Int m = H.Height();
    const Int n = Min(m,H.Width());
    const Int mRem = m-n;
    const Int bsize = Blocksize();
    H.Resize( m, n );
    MakeTrapezoidal( LOWER, H, offset );
    FillDiagonal( H, F(1), offset );

    const Grid& g = H.Grid();
    DistMatrix<F> HPanT(g), HPanB(g),
                  HEffectedNew(g), HEffectedOld(g),
                  HEffectedOldT(g), HEffectedOldB(g);

    DistMatrix<F,VC,STAR> HPan_VC_STAR(g);
    DistMatrix<F,MC,STAR> HPan_MC_STAR(g), HPanT_MC_STAR(g),
                                           HPanB_MC_STAR(g);

    DistMatrix<F,STAR,MR> Z_STAR_MR(g),
                          ZNew_STAR_MR(g), ZOld_STAR_MR(g);
    DistMatrix<F,STAR,VR> Z_STAR_VR(g),
                          ZNew_STAR_VR(g), ZOld_STAR_VR(g);
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g), SInv_STAR_STAR(g);

    Int oldEffectedHeight=mRem;

    const Int tOff = ( offset>=0 ? offset : 0 );

    const Int kLast = LastOffset( n, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Int HPanHeight = m-k;
        const Int effectedHeight = Max( HPanHeight+offset, 0 );
        const Int HPanWidth = Min( nb, effectedHeight );

        const Int effectedWidth = effectedHeight - mRem;
        const Int oldEffectedWidth = oldEffectedHeight - mRem;
        const Int newEffectedWidth = effectedWidth - oldEffectedWidth;

        auto HPan = H( IR(k,m), IR(k,k+HPanWidth) );

        auto HEffected = H( IR(m-effectedHeight,m), IR(n-effectedWidth,n) );
        PartitionLeft
        ( HEffected, HEffectedNew, HEffectedOld, oldEffectedWidth );
        PartitionUp
        ( HEffectedOld, HEffectedOldT,
                        HEffectedOldB, oldEffectedHeight );

        auto t1 = t( IR(k-tOff,k-tOff+HPanWidth), ALL );

        Z_STAR_MR.AlignWith( HEffected );
        Z_STAR_MR.Resize( HPanWidth, effectedWidth );
        PartitionLeft
        ( Z_STAR_MR, ZNew_STAR_MR, ZOld_STAR_MR, oldEffectedWidth );

        Z_STAR_VR.AlignWith( HEffected );
        Z_STAR_VR.Resize( HPanWidth, effectedWidth );
        PartitionLeft
        ( Z_STAR_VR, ZNew_STAR_VR, ZOld_STAR_VR, oldEffectedWidth );

        HPan_VC_STAR = HPan;
        Zeros( SInv_STAR_STAR, HPanWidth, HPanWidth );
        Herk
        ( UPPER, ADJOINT, 
          Base<F>(1), HPan_VC_STAR.LockedMatrix(), 
          Base<F>(0), SInv_STAR_STAR.Matrix() );
        El::AllReduce( SInv_STAR_STAR, HPan_VC_STAR.ColComm() );
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_MC_STAR.AlignWith( HEffected );
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
        Contract( ZOld_STAR_MR, ZOld_STAR_VR );
        LocalTrsm
        ( LEFT, UPPER, NORMAL, NON_UNIT, F(1), SInv_STAR_STAR, Z_STAR_VR );
        Z_STAR_MR = Z_STAR_VR;
        MakeIdentity( HEffectedNew );
        LocalGemm
        ( NORMAL, NORMAL, F(-1), HPan_MC_STAR, Z_STAR_MR, F(1), HEffected );

        oldEffectedHeight = effectedHeight;
    }

    // Take care of any untouched columns on the left side of H
    const Int oldEffectedWidth = oldEffectedHeight - mRem;
    HEffectedNew = H( ALL, IR(0,n-oldEffectedWidth) );
    MakeIdentity( HEffectedNew );
}

} // namespace expand_packed_reflectors
} // namespace El

#endif // ifndef EL_EXPANDPACKEDREFLECTORS_LV_HPP
