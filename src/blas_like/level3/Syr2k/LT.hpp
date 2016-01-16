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
inline void
LT
( T alpha,
  const ElementalMatrix<T>& APre,
  const ElementalMatrix<T>& BPre,
        ElementalMatrix<T>& CPre,
  bool conjugate=false )
{
    DEBUG_ONLY(
      CSE cse("syr2k::LT");
      AssertSameGrids( APre, BPre, CPre );
      if( APre.Width() != CPre.Height() || APre.Width() != CPre.Width()  ||
          BPre.Width() != CPre.Height() || BPre.Width() != CPre.Width()  ||
          APre.Height() != BPre.Height() )
          LogicError
          ("Nonconformal:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int r = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );
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
    DistMatrix<T,MR,  STAR> A1Trans_MR_STAR(g), B1Trans_MR_STAR(g);
    DistMatrix<T,STAR,VR  > A1_STAR_VR(g), B1_STAR_VR(g);
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g), B1_STAR_MC(g);

    A1Trans_MR_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );
    A1_STAR_MC.AlignWith( C );
    B1_STAR_MC.AlignWith( C );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        auto A1 = A( IR(k,k+nb), ALL );
        auto B1 = B( IR(k,k+nb), ALL );

        Transpose( A1, A1Trans_MR_STAR );
        Transpose( A1Trans_MR_STAR, A1_STAR_VR );
        A1_STAR_MC = A1_STAR_VR;

        Transpose( B1, B1Trans_MR_STAR );
        Transpose( B1Trans_MR_STAR, B1_STAR_VR );
        B1_STAR_MC = B1_STAR_VR;

        LocalTrr2k
        ( LOWER, orientation, TRANSPOSE, orientation, TRANSPOSE, 
          alpha,    A1_STAR_MC, B1Trans_MR_STAR,
          alphaSec, B1_STAR_MC, A1Trans_MR_STAR, T(1), C );
    }
}

} // namespace syr2k
} // namespace El
