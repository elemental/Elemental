/*
   Copyright (c) 2009-2015, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/
#ifndef EL_TRR2K_NTTN_HPP
#define EL_TRR2K_NTTN_HPP

namespace El {
namespace trr2k {

// E := alpha A B' + beta C' D + E
template<typename Ring>
void Trr2kNTTN
( UpperOrLower uplo,
  Orientation orientB, Orientation orientC,
  Ring alpha, const AbstractDistMatrix<Ring>& APre, 
              const AbstractDistMatrix<Ring>& BPre,
  Ring beta,  const AbstractDistMatrix<Ring>& CPre, 
              const AbstractDistMatrix<Ring>& DPre,
                    AbstractDistMatrix<Ring>& EPre )
{
    DEBUG_ONLY(
      CSE cse("trr2k::Trr2kNTTN");
      if( EPre.Height() != EPre.Width()  || APre.Width()  != CPre.Height() ||
          APre.Height() != EPre.Height() || CPre.Width()  != EPre.Height() ||
          BPre.Height() != EPre.Width()  || DPre.Width()  != EPre.Width()  ||
          APre.Width()  != BPre.Width()  || CPre.Height() != DPre.Height() )
          LogicError("Nonconformal Trr2kNTTN");
    )
    const Int r = APre.Width();
    const Int bsize = Blocksize();
    const Grid& g = EPre.Grid();

    auto APtr = ReadProxy<Ring,MC,MR>( &APre );      auto& A = *APtr;
    auto BPtr = ReadProxy<Ring,MC,MR>( &BPre );      auto& B = *BPtr;
    auto CPtr = ReadProxy<Ring,MC,MR>( &CPre );      auto& C = *CPtr;
    auto DPtr = ReadProxy<Ring,MC,MR>( &DPre );      auto& D = *DPtr;
    auto EPtr = ReadWriteProxy<Ring,MC,MR>( &EPre ); auto& E = *EPtr;

    DistMatrix<Ring,MC,  STAR> A1_MC_STAR(g);
    DistMatrix<Ring,VR,  STAR> B1_VR_STAR(g);
    DistMatrix<Ring,STAR,MR  > B1Trans_STAR_MR(g);
    DistMatrix<Ring,STAR,MC  > C1_STAR_MC(g);
    DistMatrix<Ring,MR,  STAR> D1Trans_MR_STAR(g);

    A1_MC_STAR.AlignWith( E );
    B1_VR_STAR.AlignWith( E );
    B1Trans_STAR_MR.AlignWith( E );
    C1_STAR_MC.AlignWith( E );
    D1Trans_MR_STAR.AlignWith( E );

    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        const Range<Int> ind1( k, k+nb );

        auto A1 = A( ALL,  ind1 );
        auto B1 = B( ALL,  ind1 );
        auto C1 = C( ind1, ALL  );
        auto D1 = D( ind1, ALL  );

        A1_MC_STAR = A1;
        C1_STAR_MC = C1;
        B1_VR_STAR = B1;
        Transpose( B1_VR_STAR, B1Trans_STAR_MR, (orientB==ADJOINT) );
        Transpose( D1, D1Trans_MR_STAR );
        LocalTrr2k 
        ( uplo, NORMAL, NORMAL, orientC, TRANSPOSE,
          alpha, A1_MC_STAR, B1Trans_STAR_MR,
          beta,  C1_STAR_MC, D1Trans_MR_STAR, Ring(1), E );
    }
}

} // namespace trr2k
} // namespace El

#endif // ifndef EL_TRR2K_NTTN_HPP
