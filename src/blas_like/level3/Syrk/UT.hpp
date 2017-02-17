/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace syrk {

template<typename T>
void UT_C
( T alpha,
  const AbstractDistMatrix<T>& APre, 
        AbstractDistMatrix<T>& CPre,
  bool conjugate=false )
{
    EL_DEBUG_CSE
    const Int r = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& C = CProx.Get();

    // Temporary distributions
    DistMatrix<T,MR,  STAR> A1Trans_MR_STAR(g);
    DistMatrix<T,STAR,VR  > A1_STAR_VR(g);
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);

    A1Trans_MR_STAR.AlignWith( C );
    A1_STAR_MC.AlignWith( C );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);
        auto A1 = A( IR(k,k+nb), ALL );

        Transpose( A1, A1Trans_MR_STAR );
        Transpose( A1Trans_MR_STAR, A1_STAR_VR );
        A1_STAR_MC = A1_STAR_VR;

        LocalTrrk
        ( UPPER, orientation, TRANSPOSE, 
          alpha, A1_STAR_MC, A1Trans_MR_STAR, T(1), C );
    }
}

template<typename T>
void UT_Dot
( T alpha,
  const AbstractDistMatrix<T>& APre,
        AbstractDistMatrix<T>& CPre,
  const bool conjugate,
  Int blockSize=2000 )
{
    EL_DEBUG_CSE
    const Int n = CPre.Height();
    const Grid& g = APre.Grid();

    const Orientation orient = ( conjugate ? ADJOINT : TRANSPOSE );

    DistMatrixReadProxy<T,T,VC,STAR> AProx( APre );
    auto& A = AProx.GetLocked();

    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& C = CProx.Get();

    DistMatrix<T,STAR,STAR> Z( blockSize, blockSize, g );
    Zero( Z );
    for( Int kOuter=0; kOuter<n; kOuter+=blockSize )
    {
        const Int nbOuter = Min(blockSize,n-kOuter);
        const Range<Int> indOuter( kOuter, kOuter+nbOuter );

        auto A1 = A( ALL, indOuter );
        auto C11 = C( indOuter, indOuter );

        Z.Resize( nbOuter, nbOuter );
        Syrk( UPPER, TRANSPOSE, alpha, A1.Matrix(), Z.Matrix(), conjugate );
        AxpyContract( T(1), Z, C11 );

        for( Int kInner=0; kInner<kOuter; kInner+=blockSize )
        {
            const Int nbInner = Min(blockSize,kOuter-kInner);
            const Range<Int> indInner( kInner, kInner+nbInner );

            auto A2 = A( ALL, indInner );
            auto C21 = C( indInner, indOuter );

            LocalGemm( orient, NORMAL, alpha, A1, A2, Z );
            AxpyContract( T(1), Z, C21 );
        }
    }
}

template<typename T>
void UT
( T alpha,
  const AbstractDistMatrix<T>& A, 
        AbstractDistMatrix<T>& C,
  bool conjugate=false )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      AssertSameGrids( A, C );
      if( A.Width() != C.Height() || A.Width() != C.Width() )
          LogicError
          ("Nonconformal:\n",DimsString(A,"A"),"\n",DimsString(C,"C"))
    )
    const Int r = A.Height();
    const Int n = A.Width();

    const double weightAwayFromDot = 10.;

    const Int blockSizeDot = 2000;

    if( r > weightAwayFromDot*n )
        UT_Dot( alpha, A, C, conjugate, blockSizeDot );
    else
        UT_C( alpha, A, C, conjugate );
}

} // namespace syrk
} // namespace El
