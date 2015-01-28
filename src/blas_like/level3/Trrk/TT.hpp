/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRRK_TT_HPP
#define EL_TRRK_TT_HPP

namespace El {
namespace trrk {

// Distributed C := alpha A^{T/H} B^{T/H} + beta C
template<typename T>
void TrrkTT
( UpperOrLower uplo,
  Orientation orientationOfA, Orientation orientationOfB,
  T alpha, const AbstractDistMatrix<T>& APre,
           const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("trrk::TrrkTN");
        if( CPre.Height() != CPre.Width() || APre.Width() != CPre.Height() || 
            BPre.Height() != CPre.Width() || APre.Height() != BPre.Width() )
            LogicError("Nonconformal TrrkTN");
        if( orientationOfA == NORMAL || orientationOfB == NORMAL )
            LogicError("Orientations must be TRANSPOSE or ADJOINT");
    )
    const Int n = CPre.Height();
    const Int r = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = CPre.Grid();

    auto APtr = ReadProxy<T,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<T,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadWriteProxy<T,MC,MR>( &CPre ); auto& C = *CPtr;

    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,VR,STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR> B1Trans_STAR_MR(g);

    A1_STAR_MC.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );

    const Range<Int> outerInd( 0, n );
    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ind1,     outerInd );
        auto B1 = B( outerInd, ind1     );

        A1_STAR_MC = A1;
        B1_VR_STAR = B1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR, (orientationOfB==ADJOINT) );
        LocalTrrk
        ( uplo, orientationOfA,
          alpha, A1_STAR_MC, B1Trans_STAR_MR, beta, C );
    }
}

} // namespace trrk
} // namespace El

#endif // ifndef EL_TRRK_TT_HPP
