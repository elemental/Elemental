/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace symm {

template<typename Ring>
void LocalAccumulateRU
( Orientation orient, Ring alpha,
  const DistMatrix<Ring,MC,  MR  >& A,
  const DistMatrix<Ring,STAR,MC  >& B_STAR_MC,
  const DistMatrix<Ring,MR,  STAR>& BTrans_MR_STAR,
        DistMatrix<Ring,MC,  STAR>& ZTrans_MC_STAR,
        DistMatrix<Ring,MR,  STAR>& ZTrans_MR_STAR )
{
    DEBUG_ONLY(
      CSE cse("symm::LocalAccumulateRU");
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

    DistMatrix<Ring> D11(g);

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> ind1( k,    k+nb ),
                         ind2( k+nb, n    );

        auto A11 = A( ind1, ind1 );
        auto A12 = A( ind1, ind2 );

        auto B1_STAR_MC = B_STAR_MC( ALL, ind1 );

        auto B1Trans_MR_STAR = BTrans_MR_STAR( ind1, ALL );
        auto B2Trans_MR_STAR = BTrans_MR_STAR( ind2, ALL );

        auto Z1Trans_MC_STAR = ZTrans_MC_STAR( ind1, ALL );

        auto Z1Trans_MR_STAR = ZTrans_MR_STAR( ind1, ALL );
        auto Z2Trans_MR_STAR = ZTrans_MR_STAR( ind2, ALL );

        D11.AlignWith( A11 );
        D11 = A11;
        MakeTrapezoidal( UPPER, D11 );
        LocalGemm
        ( alpha, D11.Orient(orient), B1_STAR_MC.Orient(orient), 
          Ring(1), Z1Trans_MR_STAR );
        FillDiagonal( D11, Ring(0) );

        LocalGemm
        ( alpha, D11.N(), B1Trans_MR_STAR.N(), Ring(1), Z1Trans_MC_STAR );

        LocalGemm
        ( alpha, A12.Orient(orient), B1_STAR_MC.Orient(orient), 
          Ring(1), Z2Trans_MR_STAR );

        LocalGemm
        ( alpha, A12.N(), B2Trans_MR_STAR.N(), Ring(1), Z1Trans_MC_STAR );
    }
}

template<typename Ring>
inline void
RUA
( Ring alpha, const AbstractDistMatrix<Ring>& APre, 
              const AbstractDistMatrix<Ring>& BPre,
                    AbstractDistMatrix<Ring>& CPre, 
  bool conjugate=false )
{
    DEBUG_ONLY(
      CSE cse("symm::RUA");
      AssertSameGrids( APre, BPre, CPre );
    )
    const Int m = CPre.Height();
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orient = ( conjugate ? ADJOINT : TRANSPOSE );

    auto APtr = ReadProxy<Ring,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<Ring,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadWriteProxy<Ring,MC,MR>( &CPre ); auto& C = *CPtr;

    DistMatrix<Ring,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<Ring,VC,  STAR> B1Trans_VC_STAR(g);
    DistMatrix<Ring,STAR,MC  > B1_STAR_MC(g);
    DistMatrix<Ring,MC,  STAR> Z1Trans_MC_STAR(g);
    DistMatrix<Ring,MR,  STAR> Z1Trans_MR_STAR(g);
    DistMatrix<Ring,MC,  MR  > Z1Trans(g);
    DistMatrix<Ring,MR,  MC  > Z1Trans_MR_MC(g);

    B1Trans_MR_STAR.AlignWith( A );
    B1Trans_VC_STAR.AlignWith( A );
    B1_STAR_MC.AlignWith( A );
    Z1Trans_MC_STAR.AlignWith( A );
    Z1Trans_MR_STAR.AlignWith( A );

    Matrix<Ring> Z1Local;

    for( Int k=0; k<m; k+=bsize )
    {
        const Int nb = Min(bsize,m-k);
        auto B1 = B( IR(k,k+nb), ALL );
        auto C1 = C( IR(k,k+nb), ALL );

        Transpose( B1, B1Trans_MR_STAR, conjugate );
        B1Trans_VC_STAR = B1Trans_MR_STAR;
        Transpose( B1Trans_VC_STAR, B1_STAR_MC, conjugate );
        Zeros( Z1Trans_MC_STAR, n, nb );
        Zeros( Z1Trans_MR_STAR, n, nb );
        LocalAccumulateRU
        ( orient, alpha, A, B1_STAR_MC, B1Trans_MR_STAR, 
          Z1Trans_MC_STAR, Z1Trans_MR_STAR );

        Contract( Z1Trans_MC_STAR, Z1Trans );
        Z1Trans_MR_MC.AlignWith( C1 );
        Z1Trans_MR_MC = Z1Trans;
        AxpyContract( Ring(1), Z1Trans_MR_STAR, Z1Trans_MR_MC );
        Transpose( Z1Trans_MR_MC.LockedMatrix(), Z1Local, conjugate );
        C1.Matrix() += Z1Local;
    }
}

template<typename Ring>
inline void
RUC
( Ring alpha, const AbstractDistMatrix<Ring>& APre, 
              const AbstractDistMatrix<Ring>& BPre,
                    AbstractDistMatrix<Ring>& CPre, 
  bool conjugate=false )
{
    DEBUG_ONLY(
        CSE cse("symm::RUC");
        AssertSameGrids( APre, BPre, CPre );
    )
    const Int n = CPre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orient = ( conjugate ? ADJOINT : TRANSPOSE );

    auto APtr = ReadProxy<Ring,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<Ring,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadWriteProxy<Ring,MC,MR>( &CPre ); auto& C = *CPtr;

    // Temporary distributions
    DistMatrix<Ring,MC,  STAR> B1_MC_STAR(g);
    DistMatrix<Ring,VR,  STAR> AT1_VR_STAR(g);
    DistMatrix<Ring,STAR,MR  > AT1Trans_STAR_MR(g);
    DistMatrix<Ring,MR,  STAR> A1RTrans_MR_STAR(g);

    B1_MC_STAR.AlignWith( C );

    for( Int k=0; k<n; k+=bsize )
    {
        const Int nb = Min(bsize,n-k);

        const Range<Int> indL( 0, k+nb ), indT( 0, k+nb ),
                         ind1( k, k+nb ),
                         indR( k, n    );

        auto A1R = A( ind1, indR );
        auto AT1 = A( indT, ind1 );

        auto B1 = B( ALL, ind1 );

        auto CL = C( ALL, indL );
        auto CR = C( ALL, indR );

        AT1_VR_STAR.AlignWith( CL );
        AT1_VR_STAR = AT1;
        AT1Trans_STAR_MR.AlignWith( CL );
        Transpose( AT1_VR_STAR, AT1Trans_STAR_MR, conjugate );
        A1RTrans_MR_STAR.AlignWith( CR );
        Transpose( A1R, A1RTrans_MR_STAR, conjugate );
        MakeTrapezoidal( LOWER, A1RTrans_MR_STAR );
        MakeTrapezoidal( LOWER, AT1Trans_STAR_MR, k-1 );

        B1_MC_STAR = B1;
        LocalGemm
        ( alpha, B1_MC_STAR.N(), A1RTrans_MR_STAR.Orient(orient), Ring(1), CR );

        LocalGemm( alpha, B1_MC_STAR.N(), AT1Trans_STAR_MR.N(), Ring(1), CL );
    }
}

template<typename Ring>
inline void
RU
( Ring alpha, const AbstractDistMatrix<Ring>& A, 
              const AbstractDistMatrix<Ring>& B,
                    AbstractDistMatrix<Ring>& C,
  bool conjugate=false )
{
    DEBUG_ONLY(CSE cse("symm::RU"))
    // TODO: Come up with a better routing mechanism
    if( A.Height() > 5*B.Height() )
        symm::RUA( alpha, A, B, C, conjugate );
    else
        symm::RUC( alpha, A, B, C, conjugate );
}

} // namespace symm
} // namespace El
