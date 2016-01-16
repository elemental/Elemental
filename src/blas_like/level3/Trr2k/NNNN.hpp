/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_NNNN_HPP
#define EL_TRR2K_NNNN_HPP

namespace El {
namespace trr2k {

// E := alpha A B + beta C D + E
template<typename T>
void Trr2kNNNN
( UpperOrLower uplo,
  T alpha,
  const ElementalMatrix<T>& APre,
  const ElementalMatrix<T>& BPre,
  T beta,
  const ElementalMatrix<T>& CPre,
  const ElementalMatrix<T>& DPre,
        ElementalMatrix<T>& EPre )
{
    DEBUG_ONLY(
      CSE cse("trr2k::Trr2kNNNN");
      if( EPre.Height() != EPre.Width()  || APre.Width()  != CPre.Width()  ||
          APre.Height() != EPre.Height() || CPre.Height() != EPre.Height() ||
          BPre.Width()  != EPre.Width()  || DPre.Width()  != EPre.Width()  ||
          APre.Width()  != BPre.Height() || CPre.Width()  != DPre.Height() )
          LogicError("Nonconformal Trr2kNNNN");
    )
    const Int r = APre.Width();
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

    DistMatrix<T,MC,STAR> A1_MC_STAR(g), C1_MC_STAR(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g), D1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( E );
    B1Trans_MR_STAR.AlignWith( E );
    C1_MC_STAR.AlignWith( E );
    D1Trans_MR_STAR.AlignWith( E );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ALL,  ind1 );
        auto B1 = B( ind1, ALL  );
        auto C1 = C( ALL,  ind1 );
        auto D1 = D( ind1, ALL  );

        A1_MC_STAR = A1;
        C1_MC_STAR = C1;
        Transpose( B1, B1Trans_MR_STAR );
        Transpose( D1, D1Trans_MR_STAR );
        LocalTrr2k
        ( uplo, NORMAL, TRANSPOSE, NORMAL, TRANSPOSE, 
          alpha, A1_MC_STAR, B1Trans_MR_STAR, 
          beta,  C1_MC_STAR, D1Trans_MR_STAR, T(1), E );
    }
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_NNNN_HPP
