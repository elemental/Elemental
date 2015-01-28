/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRRK_TN_HPP
#define EL_TRRK_TN_HPP

namespace El {
namespace trrk {

// Distributed C := alpha A^{T/H} B + beta C
template<typename T>
void TrrkTN
( UpperOrLower uplo, Orientation orientationOfA,
  T alpha, const AbstractDistMatrix<T>& APre,
           const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("trrk::TrrkTN");
        if( CPre.Height() != CPre.Width() || APre.Width() != CPre.Height() || 
            BPre.Width() != CPre.Width() || APre.Height() != BPre.Height() )
            LogicError("Nonconformal TrrkTN");
        if( orientationOfA == NORMAL )
            LogicError("Orientation must be ADJOINT or NORMAL");
    )
    const Int n = CPre.Height();
    const Int r = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = CPre.Grid();

    auto APtr = ReadProxy<T,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<T,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadWriteProxy<T,MC,MR>( &CPre ); auto& C = *CPtr;

    DistMatrix<T,STAR,MC> A1_STAR_MC(g);
    DistMatrix<T,MR,STAR> B1Trans_MR_STAR(g);

    A1_STAR_MC.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );

    const Range<Int> outerInd( 0, n );
    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ind1, outerInd );
        auto B1 = B( ind1, outerInd );

        A1_STAR_MC = A1;
        Transpose( B1, B1Trans_MR_STAR );
        LocalTrrk
        ( uplo, orientationOfA, TRANSPOSE, 
          alpha, A1_STAR_MC, B1Trans_MR_STAR, beta, C );
    }
}

} // namespace trrk
} // namespace El

#endif // ifndef EL_TRRK_TN_HPP
