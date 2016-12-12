/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace symm {

template<typename T>
void LocalAccumulateRL
( Orientation orientation, T alpha,
  const DistMatrix<T>& A,
  const DistMatrix<T,STAR,MC  >& B_STAR_MC,
  const DistMatrix<T,MR,  STAR>& BTrans_MR_STAR,
        DistMatrix<T,MC,  STAR>& ZTrans_MC_STAR,
        DistMatrix<T,MR,  STAR>& ZTrans_MR_STAR )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids
      ( A, B_STAR_MC, BTrans_MR_STAR, ZTrans_MC_STAR, ZTrans_MR_STAR );
      if( A.Height() != A.Width() ||
          A.Height() != B_STAR_MC.Width() ||
          A.Height() != BTrans_MR_STAR.Height() ||
          A.Height() != ZTrans_MC_STAR.Height() ||
          A.Height() != ZTrans_MR_STAR.Height() ||
          B_STAR_MC.Height() != BTrans_MR_STAR.Width() ||
          BTrans_MR_STAR.Width() != ZTrans_MC_STAR.Width() ||
          ZTrans_MC_STAR.Width() != ZTrans_MR_STAR.Width() )
          LogicError
          ("Nonconformal:\n",
           DimsString(A,"A"),"\n",
           DimsString(B_STAR_MC,"B[* ,MC]"),"\n",
           DimsString(BTrans_MR_STAR,"B'[MR,* ]"),"\n",
           DimsString(ZTrans_MC_STAR,"Z'[MC,* ]"),"\n",
           DimsString(ZTrans_MR_STAR,"Z'[MR,* ]"));
      if( B_STAR_MC.RowAlign() != A.ColAlign() ||
          BTrans_MR_STAR.ColAlign() != A.RowAlign() ||
          ZTrans_MC_STAR.ColAlign() != A.ColAlign() ||
          ZTrans_MR_STAR.ColAlign() != A.RowAlign() )
          LogicError("Partial matrix distributions are misaligned");
    )
    const Int n = B_STAR_MC.Width();
    const Grid& g = A.Grid();
    const Int ratio = Max( g.Height(), g.Width() );
    const Int bsize = ratio*Blocksize();

    DistMatrix<T> D11(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A21 = A( ind2, ind1 );

        auto B1_STAR_MC = B_STAR_MC( ALL, ind1 );
        auto B2_STAR_MC = B_STAR_MC( ALL, ind2 );

        auto B1Trans_MR_STAR = BTrans_MR_STAR( ind1, ALL );

        auto Z1Trans_MC_STAR = ZTrans_MC_STAR( ind1, ALL ); 
        auto Z2Trans_MC_STAR = ZTrans_MC_STAR( ind2, ALL );
   
        auto Z1Trans_MR_STAR = ZTrans_MR_STAR( ind1, ALL );

        D11.AlignWith( A11 );
        D11 = A11;
        MakeTrapezoidal( LOWER, D11 );
        LocalGemm
        ( orientation, orientation,
          alpha, D11, B1_STAR_MC, T(1), Z1Trans_MR_STAR );
        FillDiagonal( D11, T(0) );
        LocalGemm
        ( NORMAL, NORMAL, alpha, D11, B1Trans_MR_STAR, T(1), Z1Trans_MC_STAR );

        LocalGemm
        ( orientation, orientation,
          alpha, A21, B2_STAR_MC, T(1), Z1Trans_MR_STAR );

        LocalGemm
        ( NORMAL, NORMAL, alpha, A21, B1Trans_MR_STAR, T(1), Z2Trans_MC_STAR );
    }
}

template<typename T>
void RLA
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre,
  bool conjugate=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(AssertSameGrids( APre, BPre, CPre ))
    const Int m = CPre.Height();
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre ), BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,VC,  STAR> B1Trans_VC_STAR(g);
    DistMatrix<T,STAR,MC  > B1_STAR_MC(g);
    DistMatrix<T,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<T,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<T,MC,  MR  > Z1Trans(g);
    DistMatrix<T,MR,  MC  > Z1Trans_MR_MC(g);

    B1Trans_MR_STAR.AlignWith( A );
    B1Trans_VC_STAR.AlignWith( A );
    B1_STAR_MC.AlignWith( A );
    Z1Trans_MC_STAR.AlignWith( A );
    Z1Trans_MR_STAR.AlignWith( A );

    Matrix<T> Z1Local;

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);
        auto B1 = B( IR(k,k+nb), ALL );
        auto C1 = C( IR(k,k+nb), ALL );

        Z1Trans_MR_MC.AlignWith( C1 );
        Transpose( B1, B1Trans_MR_STAR, conjugate );
        B1Trans_VC_STAR = B1Trans_MR_STAR;
        Transpose( B1Trans_VC_STAR, B1_STAR_MC, conjugate );
        Z1Trans_MC_STAR.Resize( n, nb );
        Z1Trans_MR_STAR.Resize( n, nb );
        Zero( Z1Trans_MC_STAR );
        Zero( Z1Trans_MR_STAR );
        LocalAccumulateRL
        ( orientation, alpha, A, B1_STAR_MC, B1Trans_MR_STAR, 
          Z1Trans_MC_STAR, Z1Trans_MR_STAR );

        Contract( Z1Trans_MC_STAR, Z1Trans );
        Z1Trans_MR_MC = Z1Trans;
        AxpyContract( T(1), Z1Trans_MR_STAR, Z1Trans_MR_MC );
        Transpose( Z1Trans_MR_MC.LockedMatrix(), Z1Local, conjugate );
        C1.Matrix() += Z1Local;
    }
}

template<typename T>
void RLC
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre,
  bool conjugate=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(AssertSameGrids( APre, BPre, CPre ))
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre ), BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,MC,  STAR> B1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> AB1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > AB1Trans_STAR_MR(g);
    DistMatrix<T,MR,  STAR> A1LTrans_MR_STAR(g);

    B1_MC_STAR.AlignWith( C );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> indL( 0, k+nb ),
                         ind1( k, k+nb ),
                         indB( k, n    ), indR( k, n );

        auto A1L = A( ind1, indL );
        auto AB1 = A( indB, ind1 );

        auto B1 = B( ALL, ind1 );

        auto CL = C( ALL, indL );
        auto CR = C( ALL, indR );

        A1LTrans_MR_STAR.AlignWith( CL );
        Transpose( A1L, A1LTrans_MR_STAR );
        AB1_VR_STAR.AlignWith( CR );
        AB1_VR_STAR = AB1;
        AB1Trans_STAR_MR.AlignWith( CR );
        Transpose( AB1_VR_STAR, AB1Trans_STAR_MR, conjugate );
        MakeTrapezoidal( UPPER, A1LTrans_MR_STAR, -k );
        MakeTrapezoidal( UPPER, AB1Trans_STAR_MR, 1 );

        B1_MC_STAR = B1;
        LocalGemm
        ( NORMAL, TRANSPOSE, alpha, B1_MC_STAR, A1LTrans_MR_STAR, T(1), CL );

        LocalGemm
        ( NORMAL, NORMAL, alpha, B1_MC_STAR, AB1Trans_STAR_MR, T(1), CR );
    }
}

template<typename T>
void RL
( T alpha,
  const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C,
  bool conjugate=false )
{
    EL_DEBUG_CSE
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Height() )
        symm::RLA( alpha, A, B, C, conjugate );
    else
        symm::RLC( alpha, A, B, C, conjugate );
}

} // namespace symm
} // namespace El
