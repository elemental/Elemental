/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_NNNT_HPP
#define EL_TRR2K_NNNT_HPP

namespace El {
namespace trr2k {

// Distributed E := alpha (A B + C D^{T/H}) + beta E
template<typename T>
void Trr2kNNNT
( UpperOrLower uplo,
  Orientation orientationOfD,
  T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
           const AbstractDistMatrix<T>& CPre, const AbstractDistMatrix<T>& DPre,
  T beta,        AbstractDistMatrix<T>& EPre )
{
    DEBUG_ONLY(
        CallStackEntry cse("trr2k::Trr2kNNNT");
        if( EPre.Height() != EPre.Width()  || APre.Width()  != CPre.Width()  ||
            APre.Height() != EPre.Height() || CPre.Height() != EPre.Height() ||
            BPre.Width()  != EPre.Width()  || DPre.Height() != EPre.Width()  ||
            APre.Width()  != BPre.Height() || CPre.Width()  != DPre.Width() )
            LogicError("Nonconformal Trr2kNNNT");
    )
    const Int n = EPre.Height();
    const Int r = APre.Width();
    const Int bsize = Blocksize();
    const Grid& g = EPre.Grid();

    DistMatrix<T> A(g), B(g), C(g), D(g), E(g);
    Copy( APre, A, READ_PROXY );
    Copy( BPre, B, READ_PROXY );
    Copy( CPre, C, READ_PROXY );
    Copy( DPre, D, READ_PROXY );
    Copy( EPre, E, READ_WRITE_PROXY );

    DistMatrix<T,MC,  STAR> A1_MC_STAR(g), C1_MC_STAR(g);
    DistMatrix<T,MR,  STAR> B1Trans_MR_STAR(g);
    DistMatrix<T,VR,  STAR> D1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > D1Trans_STAR_MR(g);

    A1_MC_STAR.AlignWith( E );
    B1Trans_MR_STAR.AlignWith( E );
    C1_MC_STAR.AlignWith( E );
    D1_VR_STAR.AlignWith( E );
    D1Trans_STAR_MR.AlignWith( E );

    const IndexRange outerInd( 0, n );
    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const IndexRange ind1( k, k+nb );

        auto A1 = LockedView( A, outerInd, ind1     );
        auto B1 = LockedView( B, ind1,     outerInd );
        auto C1 = LockedView( C, outerInd, ind1     );
        auto D1 = LockedView( D, outerInd, ind1     );

        A1_MC_STAR = A1;
        C1_MC_STAR = C1;
        B1.TransposeColAllGather( B1Trans_MR_STAR );
        D1_VR_STAR = D1;
        D1_VR_STAR.TransposePartialColAllGather
        ( D1Trans_STAR_MR, (orientationOfD==ADJOINT) );
        LocalTrr2k
        ( uplo, TRANSPOSE, 
          alpha, A1_MC_STAR, B1Trans_MR_STAR, 
                 C1_MC_STAR, D1Trans_STAR_MR, beta, E );
    }
    Copy( E, EPre, RESTORE_READ_WRITE_PROXY );
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_NNNT_HPP
