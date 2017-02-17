/*
   Copyright (c) 2009-2016, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRRK_NT_HPP
#define EL_TRRK_NT_HPP

namespace El {
namespace trrk {

// Distributed C := alpha A B^{T/H} + C
template<typename T>
void TrrkNT
( UpperOrLower uplo,
  Orientation orientationOfB,
  T alpha,
  const AbstractDistMatrix<T>& APre,
  const AbstractDistMatrix<T>& BPre,
        AbstractDistMatrix<T>& CPre )
{
    EL_DEBUG_CSE
    EL_DEBUG_ONLY(
      if( CPre.Height() != CPre.Width() || APre.Height() != CPre.Height() || 
          BPre.Height() != CPre.Width() || APre.Width() != BPre.Width() )
          LogicError("Nonconformal TrrkNT");
      if( orientationOfB == NORMAL )
          LogicError("Orientation must be ADJOINT or TRANSPOSE");
    )
    const Int r = APre.Width();
    const Int bsize = Blocksize();
    const Grid& g = CPre.Grid();

    DistMatrixReadProxy<T,T,MC,MR> AProx( APre ), BProx( BPre );
    DistMatrixReadWriteProxy<T,T,MC,MR> CProx( CPre );
    auto& A = AProx.GetLocked();
    auto& B = BProx.GetLocked();
    auto& C = CProx.Get();

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1Trans_STAR_MR(g);

    A1_MC_STAR.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ALL, ind1 );
        auto B1 = B( ALL, ind1 );

        A1_MC_STAR = A1;
        B1_VR_STAR = B1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR, (orientationOfB==ADJOINT) );
        LocalTrrk( uplo, alpha, A1_MC_STAR, B1Trans_STAR_MR, T(1), C );
    }
}

} // namespace trrk
} // namespace El

#endif // ifndef EL_TRRK_NT_HPP
