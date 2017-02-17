/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_APPLYPACKEDREFLECTORS_LUHF_HPP
#define EL_APPLYPACKEDREFLECTORS_LUHF_HPP

namespace El {
namespace apply_packed_reflectors {

//
// Since applying Householder transforms from vectors stored top-to-bottom
// implies that we will be forming a generalization of 
//
//  (I - tau_1 v_1^T conj(v_1)) (I - tau_0 v_0^T conj(v_0)) = 
//  I - [ v_0^T, v_1^T ] [  tau_0,                       0     ] [ conj(v_0) ]
//                       [ -tau_0 tau_1 conj(v_0) v_1^T, tau_1 ] [ conj(v_1) ],
//
// which has a lower-triangular center matrix, say S, we will form S as 
// the inverse of a matrix T, which can easily be formed as
// 
//   tril(T,-1) = tril( conj(V V^H) ),
//   diag(T) = 1/householderScalars or 1/conj(householderScalars),
//
// where V is the matrix of Householder vectors and householderScalars is the
// vector of Householder reflection coeficients.
//
// V is stored row-wise in the matrix.
//

template<typename F>
void LUHFUnblocked
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( H.Width() != A.Height() )
          LogicError("H's width must match A's height");
    )
    const Int nH = H.Width();
    const Int diagLength = H.DiagonalLength(offset); 
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    Matrix<F> hPanCopy, z;

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    for( Int k=0; k<diagLength; ++k )
    {
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto hPan = H( IR(ki), IR(kj,nH) );
        auto ABot = A( IR(kj,nH), ALL );
        const F tau = householderScalars(k);
        const F gamma = ( conjugation == CONJUGATED ? Conj(tau) : tau );

        // Convert to an explicit (scaled) Householder vector
        hPanCopy = hPan;
        hPanCopy(0,0) = 1;

        // z := ABot' hPan^T 
        Gemv( ADJOINT, F(1), ABot, hPanCopy, z );
        // ABot := (I - gamma hPan^T conj(hPan)) ABot = ABot - gamma hPan^T z'
        Ger( -gamma, hPanCopy, z, ABot );
    }
}

template<typename F>
void LUHFBlocked
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( H.Width() != A.Height() )
          LogicError("H's width must match A's height");
    )
    const Int nH = H.Width();
    const Int diagLength = H.DiagonalLength(offset); 
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
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

        auto HPan = H( IR(ki,ki+nb), IR(kj,nH) );
        auto ABot = A( IR(kj,nH),    ALL       );
        auto householderScalars1 = householderScalars( IR(k,k+nb), ALL );

        // Convert to an explicit matrix of (scaled) Householder vectors
        Conjugate( HPan, HPanConj );
        MakeTrapezoidal( UPPER, HPanConj );
        FillDiagonal( HPanConj, F(1) );

        // Form the small triangular matrix needed for the UT transform
        Herk( LOWER, NORMAL, Base<F>(1), HPanConj, SInv );
        FixDiagonal( conjugation, householderScalars1, SInv );

        // Z := conj(HPan) ABot
        Gemm( NORMAL, NORMAL, F(1), HPanConj, ABot, Z );
        // Z := inv(SInv) conj(HPan) ABot
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), SInv, Z );
        // ABot := (I - HPan^T inv(SInv) conj(HPan)) ABot
        Gemm( ADJOINT, NORMAL, F(-1), HPanConj, Z, F(1), ABot );
    }
}

template<typename F>
void LUHF
( Conjugation conjugation,
  Int offset, 
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    const Int numRHS = A.Width();
    const Int blocksize = Blocksize();
    if( numRHS < blocksize )
    {
        LUHFUnblocked( conjugation, offset, H, householderScalars, A );
    }
    else
    {
        LUHFBlocked( conjugation, offset, H, householderScalars, A );
    }
}

template<typename F>
void LUHFUnblocked
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalarsPre,
        AbstractDistMatrix<F>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(AssertSameGrids( H, householderScalarsPre, APre ))

    // We gather the entire set of Householder scalars at the start rather than
    // continually paying the latency cost of the broadcasts in a 'Get' call
    DistMatrixReadProxy<F,F,STAR,STAR>
      householderScalarsProx( householderScalarsPre );
    auto& householderScalars = householderScalarsProx.GetLocked();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int nH = H.Width();
    const Int diagLength = H.DiagonalLength(offset); 
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    const Grid& g = H.Grid(); 
    auto hPan = unique_ptr<AbstractDistMatrix<F>>( H.Construct(g,H.Root()) );
    DistMatrix<F,STAR,MC> hPan_STAR_MC(g);
    DistMatrix<F,MR,STAR> z_MR_STAR(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    for( Int k=0; k<diagLength; ++k )
    {
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto ABot = A( IR(kj,nH), ALL );
        const F tau = householderScalars.GetLocal( k, 0 );
        const F gamma = ( conjugation == CONJUGATED ? Conj(tau) : tau );

        // Convert to an explicit (scaled) Householder vector
        LockedView( *hPan, H, IR(ki), IR(kj,nH) );
        hPan_STAR_MC.AlignWith( ABot );
        Conjugate( *hPan, hPan_STAR_MC );
        hPan_STAR_MC.Set( 0, 0, F(1) );

        // z := ABot' hPan^T
        z_MR_STAR.AlignWith( ABot );
        Zeros( z_MR_STAR, ABot.Width(), 1 );
        LocalGemv( ADJOINT, F(1), ABot, hPan_STAR_MC, F(0), z_MR_STAR );
        El::AllReduce( z_MR_STAR.Matrix(), ABot.ColComm() );

        // ABot := (I - gamma hPan^T conj(hPan)) ABot = ABot - gamma hPan^T z'
        LocalGer( -gamma, hPan_STAR_MC, z_MR_STAR, ABot );
    }
}

template<typename F>
void LUHFBlocked
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalarsPre,
        AbstractDistMatrix<F>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(AssertSameGrids( H, householderScalarsPre, APre ))

    DistMatrixReadProxy<F,F,MC,STAR>
      householderScalarsProx( householderScalarsPre );
    auto& householderScalars = householderScalarsProx.GetLocked();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int nH = H.Width();
    const Int diagLength = H.DiagonalLength(offset); 
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    const Grid& g = H.Grid(); 
    auto HPan = unique_ptr<AbstractDistMatrix<F>>( H.Construct(g,H.Root()) );
    DistMatrix<F> HPanConj(g);
    DistMatrix<F,STAR,VR  > HPan_STAR_VR(g);
    DistMatrix<F,STAR,MC  > HPan_STAR_MC(g);
    DistMatrix<F,STAR,STAR> householderScalars1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<F,STAR,MR  > Z_STAR_MR(g);
    DistMatrix<F,STAR,VR  > Z_STAR_VR(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    for( Int k=0; k<diagLength; k+=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto ABot = A( IR(kj,nH), ALL );
        auto householderScalars1 = householderScalars( IR(k,k+nb), ALL );

        // Convert to an explicit matrix of (scaled) Householder vectors
        LockedView( *HPan, H, IR(ki,ki+nb), IR(kj,nH) );
        Conjugate( *HPan, HPanConj );
        MakeTrapezoidal( UPPER, HPanConj );
        FillDiagonal( HPanConj, F(1) );

        // Form the small triangular matrix needed for the UT transform
        HPan_STAR_VR = HPanConj;
        Zeros( SInv_STAR_STAR, nb, nb );
        Herk
        ( LOWER, NORMAL,
          Base<F>(1), HPan_STAR_VR.LockedMatrix(),
          Base<F>(0), SInv_STAR_STAR.Matrix() );
        El::AllReduce( SInv_STAR_STAR, HPan_STAR_VR.RowComm() );
        householderScalars1_STAR_STAR = householderScalars1;
        FixDiagonal
        ( conjugation, householderScalars1_STAR_STAR, SInv_STAR_STAR );

        // Z := conj(HPan) ABot
        HPan_STAR_MC.AlignWith( ABot );
        HPan_STAR_MC = HPan_STAR_VR;
        Z_STAR_MR.AlignWith( ABot );
        LocalGemm( NORMAL, NORMAL, F(1), HPan_STAR_MC, ABot, Z_STAR_MR );
        Z_STAR_VR.AlignWith( ABot );
        Contract( Z_STAR_MR, Z_STAR_VR );

        // Z := inv(SInv) conj(HPan) ABot
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, F(1), SInv_STAR_STAR, Z_STAR_VR );

        // ABot := (I - HPan^T inv(SInv) conj(HPan)) ABot
        Z_STAR_MR = Z_STAR_VR;
        LocalGemm
        ( ADJOINT, NORMAL, F(-1), HPan_STAR_MC, Z_STAR_MR, F(1), ABot );
    }
}

template<typename F>
void LUHF
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalars,
        AbstractDistMatrix<F>& A )
{
    EL_DEBUG_CSE
    const Int numRHS = A.Width();
    const Int blocksize = Blocksize();
    if( numRHS < blocksize )
    {
        LUHFUnblocked( conjugation, offset, H, householderScalars, A );
    }
    else
    {
        LUHFBlocked( conjugation, offset, H, householderScalars, A );
    }
}

} // namespace apply_packed_reflectors
} // namespace El

#endif // ifndef EL_APPLYPACKEDREFLECTORS_LUHF_HPP
