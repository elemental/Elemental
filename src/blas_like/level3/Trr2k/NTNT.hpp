/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_NTNT_HPP
#define EL_TRR2K_NTNT_HPP

namespace El {
namespace trr2k {

// Distributed E := alpha (A B^{T/H} + C D^{T/H}) + E
template<typename T>
void Trr2kNTNT
( UpperOrLower uplo,
  Orientation orientB,
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
      if( EPre.Height() != EPre.Width()  || APre.Width()  != CPre.Width()  ||
          APre.Height() != EPre.Height() || CPre.Height() != EPre.Height() ||
          BPre.Height() != EPre.Width()  || DPre.Height() != EPre.Width()  ||
          APre.Width()  != BPre.Width()  || CPre.Width()  != DPre.Width() )
          LogicError("Nonconformal Trr2kNTNT");
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

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g), C1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g), D1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1Trans_STAR_MR(g), D1Trans_STAR_MR(g);

    A1_MC_STAR.AlignWith( E );
    B1_VR_STAR.AlignWith( E );
    B1Trans_STAR_MR.AlignWith( E );
    C1_MC_STAR.AlignWith( E );
    D1_VR_STAR.AlignWith( E );
    D1Trans_STAR_MR.AlignWith( E );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ALL, ind1 );
        auto B1 = B( ALL, ind1 );
        auto C1 = C( ALL, ind1 );
        auto D1 = D( ALL, ind1 );

        A1_MC_STAR = A1;
        C1_MC_STAR = C1;
        B1_VR_STAR = B1;
        D1_VR_STAR = D1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR, (orientB==ADJOINT) );
        Transpose( D1_VR_STAR, D1Trans_STAR_MR, (orientD==ADJOINT) );
        LocalTrr2k
        ( uplo, NORMAL, NORMAL, NORMAL, NORMAL,
          alpha, A1_MC_STAR, B1Trans_STAR_MR, 
          beta,  C1_MC_STAR, D1Trans_STAR_MR, T(1), E );
    }
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_NTNT_HPP
