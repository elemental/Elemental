/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_APPLYPACKEDREFLECTORS_RUHF_HPP
#define EL_APPLYPACKEDREFLECTORS_RUHF_HPP

namespace El {
namespace apply_packed_reflectors {

//
// Since applying Householder transforms from vectors stored top-to-bottom
// implies that we will be forming a generalization of 
//
//  (I - tau_0 v_0^T conj(v_0)) (I - tau_1 v_1^T conj(v_1)) = 
//  I - [ v_0^T, v_1^T ] [ tau_0, -tau_0 tau_1 conj(v_0) v_1^T ] [ conj(v_0) ]
//                       [ 0,      tau_1                       ] [ conj(v_1) ],
//
// which has an upper-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   triu(T) = triu( conj(V V^H) ),  diag(T) = 1/t or 1/conj(t),
//
// where V is the matrix of Householder vectors and t is the vector of scalars.
// V is stored row-wise in the matrix.
//

template<typename F>
void
RUHF
( Conjugation conjugation,
  Int offset, 
  const Matrix<F>& H,
  const Matrix<F>& t,
        Matrix<F>& A )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( A.Width() != H.Width() )
          LogicError("H and A must have the same width");
    )
    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    DEBUG_ONLY(
      if( t.Height() != diagLength )
          LogicError("t must be the same length as H's offset diag");
    )
    Matrix<F> HPanConj, SInv, Z;

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    for( Int k=0; k<diagLength; k+=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto HPan   = H( IR(ki,ki+nb), IR(kj,nA) );
        auto ARight = A( ALL,          IR(kj,nA) );
        auto t1     = t( IR(k,k+nb),   ALL       );

        Conjugate( HPan, HPanConj );
        MakeTrapezoidal( UPPER, HPanConj );
        FillDiagonal( HPanConj, F(1) );

        Herk( UPPER, NORMAL, Base<F>(1), HPanConj, SInv );
        FixDiagonal( conjugation, t1, SInv );

        Gemm( NORMAL, ADJOINT, F(1), ARight, HPanConj, Z );
        Trsm( RIGHT, UPPER, NORMAL, NON_UNIT, F(1), SInv, Z );
        Gemm( NORMAL, NORMAL, F(-1), Z, HPanConj, F(1), ARight );
    }
}

template<typename F>
void
RUHF
( Conjugation conjugation,
  Int offset, 
  const ElementalMatrix<F>& HPre,
  const ElementalMatrix<F>& tPre, 
        ElementalMatrix<F>& APre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( APre.Width() != HPre.Width() )
          LogicError("A and H must have the same width");
      AssertSameGrids( HPre, tPre, APre );
    )

    DistMatrixReadProxy<F,F,MC,MR  > HProx( HPre );
    DistMatrixReadProxy<F,F,MC,STAR> tProx( tPre );
    DistMatrixReadWriteProxy<F,F,MC,MR  > AProx( APre );
    auto& H = HProx.GetLocked();
    auto& t = tProx.GetLocked();
    auto& A = AProx.Get();

    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    DEBUG_ONLY(
      if( t.Height() != diagLength )
          LogicError("t must be the same length as H's offset diag");
    )
    const Grid& g = H.Grid();
    DistMatrix<F> HPanConj(g);
    DistMatrix<F,STAR,VR  > HPan_STAR_VR(g);
    DistMatrix<F,STAR,MR  > HPan_STAR_MR(g);
    DistMatrix<F,STAR,STAR> t1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > ZAdj_STAR_MC(g);
    DistMatrix<F,STAR,VC  > ZAdj_STAR_VC(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    for( Int k=0; k<diagLength; k+=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto HPan   = H( IR(ki,ki+nb), IR(kj,nA) );
        auto ARight = A( ALL,          IR(kj,nA) );
        auto t1     = t( IR(k,k+nb),   ALL       );

        Conjugate( HPan, HPanConj );
        MakeTrapezoidal( UPPER, HPanConj );
        FillDiagonal( HPanConj, F(1) );

        HPan_STAR_VR = HPanConj;
        Zeros( SInv_STAR_STAR, nb, nb );
        Herk
        ( UPPER, NORMAL,
          Base<F>(1), HPan_STAR_VR.LockedMatrix(),
          Base<F>(0), SInv_STAR_STAR.Matrix() );
        El::AllReduce( SInv_STAR_STAR, HPan_STAR_VR.RowComm() );
        t1_STAR_STAR = t1;
        FixDiagonal( conjugation, t1_STAR_STAR, SInv_STAR_STAR );

        HPan_STAR_MR.AlignWith( ARight );
        HPan_STAR_MR = HPan_STAR_VR;
        ZAdj_STAR_MC.AlignWith( ARight );
        LocalGemm( NORMAL, ADJOINT, F(1), HPan_STAR_MR, ARight, ZAdj_STAR_MC );
        ZAdj_STAR_VC.AlignWith( ARight );
        Contract( ZAdj_STAR_MC, ZAdj_STAR_VC );

        LocalTrsm
        ( LEFT, UPPER, ADJOINT, NON_UNIT, F(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, NORMAL, F(-1), ZAdj_STAR_MC, HPan_STAR_MR, F(1), ARight );
    }
}

} // namespace apply_packed_reflectors
} // namespace El

#endif // ifndef EL_APPLYPACKEDREFLECTORS_RUHF_HPP
