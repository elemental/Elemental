/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_APPLYPACKEDREFLECTORS_RLVB_HPP
#define EL_APPLYPACKEDREFLECTORS_RLVB_HPP

namespace El {
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
//   tril(T,-1) = tril( U^H U ),
//   diag(T) = 1/householderScalars or 1/conj(householderScalars),
//
// where U is the matrix of Householder vectors and householderScalars is the
// vector of Householder reflection coefficients.
//

template<typename F> 
void
RLVBUnblocked
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( H.Height() != A.Width() )
          LogicError("H's height must match A's width");
    )
    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    Matrix<F> hPanCopy, z;

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    for( Int k=diagLength-1; k>=0; --k )
    {
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto hPan = H( IR(ki,nA), IR(kj) );
        auto ARight = A( ALL, IR(ki,nA) );
        const F tau = householderScalars(k);
        const F gamma = ( conjugation == CONJUGATED ? Conj(tau) : tau );

        // Convert to an explicit (scaled) Householder vector
        hPanCopy = hPan;
        hPanCopy(0) = F(1);

        // z := ARight hPan
        Gemv( NORMAL, F(1), ARight, hPanCopy, z );
        // ARight := ARight (I - gamma hPan hPan')
        Ger( -gamma, z, hPanCopy, ARight );
    }
} 

template<typename F> 
void
RLVBBlocked
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( H.Height() != A.Width() )
          LogicError("H's height must match A's width");
    )
    const Int nA = A.Width();
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
    const Int kLast = LastOffset( diagLength, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto HPan   = H( IR(ki,nA),  IR(kj,kj+nb) );
        auto ARight = A( ALL,        IR(ki,nA)    );
        auto householderScalars1 = householderScalars( IR(k,k+nb), ALL );

        // Convert to an explicit matrix of (scaled) Householder vectors
        HPanCopy = HPan;
        MakeTrapezoidal( LOWER, HPanCopy );
        FillDiagonal( HPanCopy, F(1) );

        // Form the small triangular matrix needed for the UT transform
        Herk( LOWER, ADJOINT, Base<F>(1), HPanCopy, SInv );
        FixDiagonal( conjugation, householderScalars1, SInv );

        // Z := ARight HPan
        Gemm( NORMAL, NORMAL, F(1), ARight, HPanCopy, Z );
        // Z := ARight HPan inv(SInv)
        Trsm( RIGHT, LOWER, NORMAL, NON_UNIT, F(1), SInv, Z );
        // ARight := ARight (I - HPan inv(SInv) HPan')
        Gemm( NORMAL, ADJOINT, F(-1), Z, HPanCopy, F(1), ARight );
    }
}

template<typename F> 
void
RLVB
( Conjugation conjugation,
  Int offset,
  const Matrix<F>& H,
  const Matrix<F>& householderScalars,
        Matrix<F>& A )
{
    EL_DEBUG_CSE
    const Int numLHS = A.Height();
    const Int blocksize = Blocksize();
    if( numLHS < blocksize )
    {
        RLVBUnblocked( conjugation, offset, H, householderScalars, A );
    }
    else
    {
        RLVBBlocked( conjugation, offset, H, householderScalars, A );
    }
}

template<typename F> 
void
RLVBUnblocked
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

    const Int nA = A.Width();
    const Int diagLength = H.DiagonalLength(offset);
    EL_DEBUG_ONLY(
      if( householderScalars.Height() != diagLength )
          LogicError
          ("householderScalars must be the same length as H's offset diag");
    )
    const Grid& g = H.Grid();
    auto hPan = unique_ptr<AbstractDistMatrix<F>>( H.Construct(g,H.Root()) );
    DistMatrix<F,MR,STAR> hPan_MR_STAR(g);
    DistMatrix<F,MC,STAR> z_MC_STAR(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    for( Int k=diagLength-1; k>=0; --k )
    {
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto ARight = A( ALL, IR(ki,nA) ); 
        const F tau = householderScalars.GetLocal( k, 0 );
        const F gamma = ( conjugation == CONJUGATED ? Conj(tau) : tau );

        // Convert to an explicit (scaled) Householder vector
        LockedView( *hPan, H, IR(ki,nA), IR(kj) );
        hPan_MR_STAR.AlignWith( ARight );
        Copy( *hPan, hPan_MR_STAR );
        hPan_MR_STAR.Set( 0, 0, F(1) );

        // z := ARight hPan
        z_MC_STAR.AlignWith( ARight );
        Zeros( z_MC_STAR, ARight.Height(), 1 );
        LocalGemv( NORMAL, F(1), ARight, hPan_MR_STAR, F(0), z_MC_STAR );
        El::AllReduce( z_MC_STAR, ARight.RowComm() );
 
        // ARight := ARight (I - hPan inv(SInv) hPan')
        LocalGer( -gamma, z_MC_STAR, hPan_MR_STAR, ARight );
    }
}

template<typename F> 
void
RLVBBlocked
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

    const Int nA = A.Width();
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
    DistMatrix<F,MR,  STAR> HPan_MR_STAR(g);
    DistMatrix<F,STAR,STAR> householderScalars1_STAR_STAR(g);
    DistMatrix<F,STAR,STAR> SInv_STAR_STAR(g);
    DistMatrix<F,STAR,MC  > ZAdj_STAR_MC(g);
    DistMatrix<F,STAR,VC  > ZAdj_STAR_VC(g);

    const Int iOff = ( offset>=0 ? 0      : -offset );
    const Int jOff = ( offset>=0 ? offset : 0       );

    const Int bsize = Blocksize();
    const Int kLast = LastOffset( diagLength, bsize );
    for( Int k=kLast; k>=0; k-=bsize )
    {
        const Int nb = Min(bsize,diagLength-k);
        const Int ki = k+iOff;
        const Int kj = k+jOff;

        auto ARight = A( ALL, IR(ki,nA) ); 
        auto householderScalars1 = householderScalars( IR(k,k+nb), ALL );

        // Convert to an explicit matrix of (scaled) Householder vectors
        LockedView( *HPan, H, IR(ki,nA), IR(kj,kj+nb) );
        Copy( *HPan, HPanCopy );
        MakeTrapezoidal( LOWER, HPanCopy );
        FillDiagonal( HPanCopy, F(1) );

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

        // Z := ARight HPan
        HPan_MR_STAR.AlignWith( ARight );
        HPan_MR_STAR = HPan_VC_STAR;
        ZAdj_STAR_MC.AlignWith( ARight );
        LocalGemm( ADJOINT, ADJOINT, F(1), HPan_MR_STAR, ARight, ZAdj_STAR_MC );
        ZAdj_STAR_VC.AlignWith( ARight );
        Contract( ZAdj_STAR_MC, ZAdj_STAR_VC );
 
        // Z := ARight HPan inv(SInv)
        LocalTrsm
        ( LEFT, LOWER, ADJOINT, NON_UNIT, F(1), SInv_STAR_STAR, ZAdj_STAR_VC );

        // ARight := ARight (I - HPan inv(SInv) HPan')
        ZAdj_STAR_MC = ZAdj_STAR_VC;
        LocalGemm
        ( ADJOINT, ADJOINT, F(-1), ZAdj_STAR_MC, HPan_MR_STAR, F(1), ARight );
    }
}

template<typename F> 
void
RLVB
( Conjugation conjugation,
  Int offset,
  const AbstractDistMatrix<F>& H,
  const AbstractDistMatrix<F>& householderScalars,
        AbstractDistMatrix<F>& A )
{
    EL_DEBUG_CSE
    const Int numLHS = A.Height();
    const Int blocksize = Blocksize();
    if( numLHS < blocksize )
    {
        RLVBUnblocked( conjugation, offset, H, householderScalars, A );
    }
    else
    {
        RLVBBlocked( conjugation, offset, H, householderScalars, A );
    }
}

} // namespace apply_packed_reflectors
} // namespace El

#endif // ifndef EL_APPLYPACKEDREFLECTORS_RLVB_HPP
