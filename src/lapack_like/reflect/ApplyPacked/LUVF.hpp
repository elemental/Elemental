/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_APPLYPACKEDREFLECTORS_LUVF_HPP
#define EL_APPLYPACKEDREFLECTORS_LUVF_HPP

namespace El {
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
//   tril(T,-1) = tril( U^H U ),
//   diag(T) = 1/householderScalars or 1/conj(householderScalars),
//
// where U is the matrix of Householder vectors and householderScalars is the
// vector of Householder reflection coefficients.
//

template<typename F>
void LUVFUnblocked
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( H.Height() != A.Height() )
          LogicError("H and A must be the same height");
    )
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

        auto hPan = H( IR(0,ki+1), IR(kj) );
        auto ATop = A( IR(0,ki+1), ALL );
        const F tau = householderScalars(k);
        const F gamma = ( conjugation == CONJUGATED ? Conj(tau) : tau );

        // Convert to an explicit (scaled) Householder vector
        hPanCopy = hPan;
        hPanCopy(ki) = 1;

        // z := ATop' hPan
        Gemv( ADJOINT, F(1), ATop, hPanCopy, z );
        // ATop := (I - gamma hPan hPan') ATop
        Ger( -gamma, hPanCopy, z, ATop );
    }
}

template<typename F>
void LUVFBlocked
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( H.Height() != A.Height() )
          LogicError("H and A must be the same height");
    )
    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    Matrix<F> HPanCopy, SInv, Z;

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    for( Int k=0; k<diagLength; k+=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto HPan = H( IR(0,ki+nb), IR(kj,kj+nb) );
        auto ATop = A( IR(0,ki+nb), ALL          );
        auto householderScalars1 = householderScalars( IR(k,k+nb), ALL );

        // Convert to an explicit matrix of (scaled) Householder vectors
        HPanCopy = HPan;
        MakeTrapezoidal( UPPER, HPanCopy, HPanCopy.Width()-HPanCopy.Height() );
        FillDiagonal( HPanCopy, F(1), HPanCopy.Width()-HPanCopy.Height() );

        // Form the small triangular matrix needed for the UT transform
        Herk( LOWER, ADJOINT, Base<F>(1), HPanCopy, SInv );
        FixDiagonal( conjugation, householderScalars1, SInv );

        // Z := HPan' ATop
        Gemm( ADJOINT, NORMAL, F(1), HPanCopy, ATop, Z );
        // Z := inv(SInv) HPan' ATop
        Trsm( LEFT, LOWER, NORMAL, NON_UNIT, F(1), SInv, Z );
        // ATop := (I - HPan inv(SInv) HPan') ATop
        Gemm( NORMAL, NORMAL, F(-1), HPanCopy, Z, F(1), ATop );
    }
}

template<typename F>
void LUVF
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
        LUVFUnblocked( conjugation, offset, H, householderScalars, A );
    }
    else
    {
        LUVFBlocked( conjugation, offset, H, householderScalars, A );
    }
}

template<typename F> 
void LUVFUnblocked
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalarsPre,
        AbstractDistMatrix<F>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( H, householderScalarsPre, APre );
      if( H.Height() != APre.Height() )
          LogicError("H and A must be the same height");
    )

    // We gather the entire set of Householder scalars at the start rather than
    // continually paying the latency cost of the broadcasts in a 'Get' call
    DistMatrixReadProxy<F,F,STAR,STAR>
      householderScalarsProx( householderScalarsPre );
    auto& householderScalars = householderScalarsProx.GetLocked();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    const Grid& g = H.Grid();
    auto hPan = unique_ptr<AbstractDistMatrix<F>>( H.Construct(g,H.Root()) );
    DistMatrix<F,MC,STAR> hPan_MC_STAR(g);
    DistMatrix<F,MR,STAR> z_MR_STAR(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    for( Int k=0; k<diagLength; ++k )
    {
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto ATop = A( IR(0,ki+1), ALL );
        const F tau = householderScalars.GetLocal( k, 0 );
        const F gamma = ( conjugation == CONJUGATED ? Conj(tau) : tau );

        // Convert to an explicit (scaled) Householder vector
        LockedView( *hPan, H, IR(0,ki+1), IR(kj) );
        hPan_MC_STAR.AlignWith( ATop );
        Copy( *hPan, hPan_MC_STAR );
        hPan_MC_STAR.Set( ki, 0, F(1) );

        // z := ATop' hPan
        z_MR_STAR.AlignWith( ATop );
        Zeros( z_MR_STAR, ATop.Width(), 1 );
        LocalGemv( ADJOINT, F(1), ATop, hPan_MC_STAR, F(0), z_MR_STAR );
        El::AllReduce( z_MR_STAR.Matrix(), ATop.ColComm() );

        // ATop := (I - gamma hPan hPan') ATop
        LocalGer( -gamma, hPan_MC_STAR, z_MR_STAR, ATop );
    }
}

template<typename F> 
void LUVFBlocked
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalarsPre,
        AbstractDistMatrix<F>& APre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( H, householderScalarsPre, APre );
      if( H.Height() != APre.Height() )
          LogicError("H and A must be the same height");
    )

    DistMatrixReadProxy<F,F,MC,STAR>
      householderScalarsProx( householderScalarsPre );
    auto& householderScalars = householderScalarsProx.GetLocked();

    DistMatrixReadWriteProxy<F,F,MC,MR> AProx( APre );
    auto& A = AProx.Get();

    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    const Grid& g = H.Grid();
    auto HPan = unique_ptr<AbstractDistMatrix<F>>( H.Construct(g,H.Root()) );
    DistMatrix<F> HPanCopy(g);
    DistMatrix<F,VC,  STAR> HPan_VC_STAR(g);
    DistMatrix<F,MC,  STAR> HPan_MC_STAR(g);
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

        auto ATop = A( IR(0,ki+nb), ALL );
        auto householderScalars1 = householderScalars( IR(k,k+nb), ALL );

        // Convert to an explicit matrix of (scaled) Householder vectors
        LockedView( *HPan, H, IR(0,ki+nb), IR(kj,kj+nb) );
        Copy( *HPan, HPanCopy );
        MakeTrapezoidal( UPPER, HPanCopy, HPanCopy.Width()-HPanCopy.Height() );
        FillDiagonal( HPanCopy, F(1), HPanCopy.Width()-HPanCopy.Height() );

        // Form the small triangular matrix needed for the UT transform
        HPan_VC_STAR = HPanCopy;
        Zeros( SInv_STAR_STAR, nb, nb );
        Herk
        ( LOWER, ADJOINT, 
          Base<F>(1), HPan_VC_STAR.LockedMatrix(),
          Base<F>(0), SInv_STAR_STAR.Matrix() ); 
        El::AllReduce( SInv_STAR_STAR, HPan_VC_STAR.ColComm() );
        householderScalars1_STAR_STAR = householderScalars1;
        FixDiagonal
        ( conjugation, householderScalars1_STAR_STAR, SInv_STAR_STAR );

        // Z := HPan' ATop
        HPan_MC_STAR.AlignWith( ATop );
        HPan_MC_STAR = HPanCopy;
        Z_STAR_MR.AlignWith( ATop );
        LocalGemm( ADJOINT, NORMAL, F(1), HPan_MC_STAR, ATop, Z_STAR_MR );
        Z_STAR_VR.AlignWith( ATop );
        Contract( Z_STAR_MR, Z_STAR_VR );
        
        // Z := inv(SInv) HPan' ATop
        LocalTrsm
        ( LEFT, LOWER, NORMAL, NON_UNIT, F(1), SInv_STAR_STAR, Z_STAR_VR );

        // ATop := (I - HPan inv(SInv) HPan') ATop
        Z_STAR_MR = Z_STAR_VR;
        LocalGemm( NORMAL, NORMAL, F(-1), HPan_MC_STAR, Z_STAR_MR, F(1), ATop );
    }
}

template<typename F> 
void LUVF
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
        LUVFUnblocked( conjugation, offset, H, householderScalars, A );
    }
    else
    {
        LUVFBlocked( conjugation, offset, H, householderScalars, A );
    }
}

} // namespace apply_packed_reflectors
} // namespace El

#endif // ifndef EL_APPLYPACKEDREFLECTORS_LUVF_HPP
