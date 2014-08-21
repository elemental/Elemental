/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_NTTN_HPP
#define EL_TRR2K_NTTN_HPP

namespace El {
namespace trr2k {

// Distributed E := alpha (A B^{T/H} + C^{T/H} D) + beta E
template<typename T>
void Trr2kNTTN
( UpperOrLower uplo,
  Orientation orientationOfB, Orientation orientationOfC,
  T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
           const AbstractDistMatrix<T>& CPre, const AbstractDistMatrix<T>& DPre,
  T beta,        AbstractDistMatrix<T>& EPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("trr2k::Trr2kNTTN");
        if( EPre.Height() != EPre.Width()  || APre.Width()  != CPre.Height() ||
            APre.Height() != EPre.Height() || CPre.Width()  != EPre.Height() ||
            BPre.Height() != EPre.Width()  || DPre.Width()  != EPre.Width()  ||
            APre.Width()  != BPre.Width()  || CPre.Height() != DPre.Height() )
            LogicError("Nonconformal Trr2kNTTN");
    )
    const Int n = EPre.Height();
    const Int r = APre.Width();
    const Int bsize = Blocksize();
    const Grid& g = EPre.Grid();

    auto APtr = ReadProxy( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadProxy( &CPre );      auto& C = *CPtr;
    auto DPtr = ReadProxy( &DPre );      auto& D = *DPtr;
    auto EPtr = ReadWriteProxy( &EPre ); auto& E = *EPtr;

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > B1Trans_STAR_MR(g);
    DistMatrix<T,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<T,MR,  STAR> D1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( E );
    B1_VR_STAR.AlignWith( E );
    B1Trans_STAR_MR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1Trans_MR_STAR.AlignWith( E );

    const Range<Int> outerInd( 0, n );
    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( outerInd, ind1     );
        auto B1 = B( outerInd, ind1     );
        auto C1 = C( ind1,     outerInd );
        auto D1 = D( ind1,     outerInd );

        A1_MC_STAR = A1;
        C1_STAR_MC = C1;
        B1_VR_STAR = B1;
        B1_VR_STAR.TransposePartialColAllGather
        ( B1Trans_STAR_MR, (orientationOfB==ADJOINT) );
        D1.TransposeColAllGather( D1Trans_MR_STAR );
        LocalTrr2k 
        ( uplo, orientationOfC, TRANSPOSE,
          alpha, A1_MC_STAR, B1Trans_STAR_MR,
                 C1_STAR_MC, D1Trans_MR_STAR, beta, E );
    }
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_NTTN_HPP
