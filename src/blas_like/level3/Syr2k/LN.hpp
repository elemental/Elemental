/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace syr2k {

template<typename T>
void LN_C
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre,
  bool conjugate=false )
{
    EL_DEBUG_CSE
    const Int r = APre.Width();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const T alphaSec = ( conjugate ? Conj(alpha) : alpha );

    DistMatrixReadProxy<T,T,MC,MR>
      AProx( APre ),
      BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR>
      CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,MC,  STAR> A1_MC_STAR(g), B1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> A1_VR_STAR(g), B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > A1Trans_STAR_MR(g), B1Trans_STAR_MR(g);

    A1_MC_STAR.AlignWith( C );
    B1_MC_STAR.AlignWith( C );
    A1_VR_STAR.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    A1Trans_STAR_MR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        auto A1 = A( ALL, IR(k,k+nb) );
        auto B1 = B( ALL, IR(k,k+nb) );

        A1_VR_STAR = A1_MC_STAR = A1;
        Transpose( A1_VR_STAR, A1Trans_STAR_MR, conjugate );

        B1_VR_STAR = B1_MC_STAR = B1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR, conjugate );

        LocalTrr2k
        ( LOWER, NORMAL, NORMAL, NORMAL, NORMAL,
          alpha,    A1_MC_STAR, B1Trans_STAR_MR,
          alphaSec, B1_MC_STAR, A1Trans_STAR_MR, T(1), C );
    }
}

template<typename T>
void LN_Dot
( T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre,
  const bool conjugate,
  Int blockSize=2000 )
{
    EL_DEBUG_CSE 
    const Int n = CPre.Height();
    const Grid& g = APre.Grid();

    const Orientation orient = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrixReadProxy<T,T,STAR,VC> AProx( APre );
    auto& A = AProx.GetLocked();

    ElementalProxyCtrl BCtrl;
    BCtrl.rowConstrain = true;
    BCtrl.rowAlign = A.RowAlign();
    DistMatrixReadProxy<T,T,STAR,VC> BProx( BPre, BCtrl );
    auto& B = BProx.GetLocked();

    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& C = CProx.Get();

    DistMatrix<T,STAR,STAR> Z( blockSize, blockSize, g );
    Zero( Z );
    for( Int kOuter=0; kOuter<n; kOuter+=blockSize )
    {
        const Int nbOuter = Min(blockSize,n-kOuter);
        const Range<Int> indOuter( kOuter, kOuter+nbOuter );

        auto A1 = A( indOuter, ALL );
        auto B1 = B( indOuter, ALL );
        auto C11 = C( indOuter, indOuter );

        Z.Resize( nbOuter, nbOuter );
        Syr2k
        ( LOWER, NORMAL, alpha, A1.Matrix(), B1.Matrix(), Z.Matrix(),
          conjugate );
        AxpyContract( T(1), Z, C11 );

        for( Int kInner=kOuter+nbOuter; kInner<n; kInner+=blockSize )
        {
            const Int nbInner = Min(blockSize,n-kInner);
            const Range<Int> indInner( kInner, kInner+nbInner );

            auto A2 = A( indInner, ALL );
            auto B2 = B( indInner, ALL );
            auto C21 = C( indInner, indOuter );

            LocalGemm( NORMAL, orient, alpha, A1, B2, Z );
            LocalGemm( NORMAL, orient, Conj(alpha), B1, A2, Z );
            AxpyContract( T(1), Z, C21 );
        }
    }
}

template<typename T>
void LN
( T alpha,
  const AbstractDistMatrix<T>& A,
  const AbstractDistMatrix<T>& B,
        AbstractDistMatrix<T>& C,
  bool conjugate=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, B, C );
      if( A.Height() != C.Height() || A.Height() != C.Width() ||
          B.Height() != C.Height() || B.Height() != C.Width() ||
          A.Width() != B.Width() )
          LogicError
          ("Nonconformal:\n",
           DimsString(A,"A"),"\n",
           DimsString(B,"B"),"\n",
           DimsString(C,"C"));
    )
    const Int n = A.Height();
    const Int r = A.Width();

    const double weightAwayFromDot = 10.;

    const Int blockSizeDot = 2000;

    if( r > weightAwayFromDot*n )
        LN_Dot( alpha, A, B, C, conjugate, blockSizeDot );
    else
        LN_C( alpha, A, B, C, conjugate );
}

} // namespace syr2k
} // namespace El
