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
UN
( T alpha,
  const ElementalMatrix<T>& APre,
  const ElementalMatrix<T>& BPre,
        ElementalMatrix<T>& CPre,
  bool conjugate=false )
{
    DEBUG_ONLY(
      CSE cse("syr2k::UN");
      AssertSameGrids( APre, BPre, CPre );
      if( APre.Height() != CPre.Height() || APre.Height() != CPre.Width() ||
          BPre.Height() != CPre.Height() || BPre.Height() != CPre.Width() ||
          APre.Width() != BPre.Width() )
          LogicError
          ("Nonconformal:\n",
           DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
           DimsString(CPre,"C"));
    )
    const Int r = APre.Width();
    const Int bsize = Blocksize();
    const Grid& g = CPre.Grid();
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
        ( UPPER, NORMAL, NORMAL, NORMAL, NORMAL,
          alpha,    A1_MC_STAR, B1Trans_STAR_MR,
          alphaSec, B1_MC_STAR, A1Trans_STAR_MR, T(1), C );
    }
}

} // namespace syr2k
} // namespace El
