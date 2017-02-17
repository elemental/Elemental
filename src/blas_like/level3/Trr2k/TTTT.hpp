/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_TTTT_HPP
#define EL_TRR2K_TTTT_HPP

namespace El {
namespace trr2k {

// E := alpha A' B' + beta C' D' + E
template<typename T>
void Trr2kTTTT
( UpperOrLower uplo,
  Orientation orientA,
  Orientation orientB,
  Orientation orientC,
  Orientation orientD, 
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
  T beta,
  const AbstractDistMatrix<T>& CPre,
  const AbstractDistMatrix<T>& DPre,
        AbstractDistMatrix<T>& EPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( EPre.Height() != EPre.Width()  || APre.Height() != CPre.Height() ||
          APre.Width()  != EPre.Height() || CPre.Width()  != EPre.Height() ||
          BPre.Height() != EPre.Width()  || DPre.Height() != EPre.Width()  ||
          APre.Height() != BPre.Width()  || CPre.Height() != DPre.Width() )
          LogicError("Nonconformal Trr2kTTTT");
    )
    const Int r = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = EPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR>
      AProx( APre ),
      BProx( BPre ),
      CProx( CPre ),
      DProx( DPre );
    DistMatrixReadWriteProxy<T,T,MC,MR>
      EProx( EPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.GetLocked();
    auto& D = DProx.GetLocked();
    auto& E = EProx.Get();

    DistMatrix<T,STAR,MC  > A1_STAR_MC(g), C1_STAR_MC(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g), D1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1Trans_STAR_MR(g), D1Trans_STAR_MR(g);

    A1_STAR_MC.AlignWith( E );
    B1_VR_STAR.AlignWith( E );
    B1Trans_STAR_MR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1_VR_STAR.AlignWith( E );
    D1Trans_STAR_MR.AlignWith( E );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ind1, ALL  );
        auto B1 = B( ALL,  ind1 );
        auto C1 = C( ind1, ALL  );
        auto D1 = D( ALL,  ind1 );

        A1_STAR_MC = A1;
        C1_STAR_MC = C1;
        B1_VR_STAR = B1;
        D1_VR_STAR = D1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR, (orientB==ADJOINT) );
        Transpose( D1_VR_STAR, D1Trans_STAR_MR, (orientD==ADJOINT) );
        LocalTrr2k
        ( uplo, orientA, NORMAL, orientC, NORMAL,
          alpha, A1_STAR_MC, B1Trans_STAR_MR,
          beta,  C1_STAR_MC, D1Trans_STAR_MR, T(1), E );
    }
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_TTTT_HPP
