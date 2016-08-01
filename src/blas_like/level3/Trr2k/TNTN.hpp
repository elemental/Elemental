/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_TNTN_HPP
#define EL_TRR2K_TNTN_HPP

namespace El {
namespace trr2k {

// E := alpha A' B + beta C' D + E
template<typename T>
void Trr2kTNTN
( UpperOrLower uplo,
  Orientation orientA,
  Orientation orientC,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
  T beta,
  const AbstractDistMatrix<T>& CPre,
  const AbstractDistMatrix<T>& DPre,
        AbstractDistMatrix<T>& EPre )
{
    DEBUG_CSE
    DEBUG_ONLY(
      if( EPre.Height() != EPre.Width()  || APre.Height() != CPre.Height() ||
          APre.Width()  != EPre.Height() || CPre.Width()  != EPre.Height() ||
          BPre.Width()  != EPre.Width()  || DPre.Width()  != EPre.Width()  ||
          APre.Height() != BPre.Height() || CPre.Height() != DPre.Height() )
          LogicError("Nonconformal Trr2kTNTN");
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

    DistMatrix<T,STAR,MC  > A1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( E );
    B1Trans_MR_STAR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1Trans_MR_STAR.AlignWith( E );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ind1, ALL );
        auto B1 = B( ind1, ALL );
        auto C1 = C( ind1, ALL );
        auto D1 = D( ind1, ALL );

        A1_STAR_MC = A1;
        C1_STAR_MC = C1;
        Transpose( B1, B1Trans_MR_STAR );
        Transpose( D1, D1Trans_MR_STAR );
        LocalTrr2k
        ( uplo, orientA, TRANSPOSE, orientC, TRANSPOSE,
          alpha, A1_STAR_MC, B1Trans_MR_STAR, 
          beta,  C1_STAR_MC, D1Trans_MR_STAR, T(1), E );
    }
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_TNTN_HPP
