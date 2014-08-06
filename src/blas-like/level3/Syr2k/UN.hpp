/*
   Copyright (c) 2009-2014, Jack Poulson
   All rights reserved.

   This file is part of Elemental and is under the BSD 2-Clause License, 
   which can be found in the LICENSE file in the root directory, or at 
   http://opensource.org/licenses/BSD-2-Clause
*/

namespace El {
namespace syr2k {

template<typename T>
inline void
UN
( T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("syr2k::UN");
        if( APre.Grid() != BPre.Grid() || BPre.Grid() != CPre.Grid() )
            LogicError("{A,B,C} must be distributed over the same grid");
        if( APre.Height() != CPre.Height() || APre.Height() != CPre.Width() ||
            BPre.Height() != CPre.Height() || BPre.Height() != CPre.Width() ||
            APre.Width() != BPre.Width() )
            LogicError
            ("Nonconformal:\n",
             DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
             DimsString(CPre,"C"));
    )
    const Int n = APre.Height();
    const Int r = APre.Width();
    const Int bsize = Blocksize();
    const Grid& g = CPre.Grid();

    // Force 'A', 'B', and 'C' to be in [MC,MR] distributions
    DistMatrix<T> A(g), B(g), C(g);
    Copy( APre, A, READ_PROXY );
    Copy( BPre, B, READ_PROXY );
    Copy( CPre, C, READ_WRITE_PROXY );

    // Temporary distributions
    DistMatrix<T,MC,  STAR> A1_MC_STAR(g), B1_MC_STAR(g);
    DistMatrix<T,VR,  STAR> A1_VR_STAR(g), B1_VR_STAR(g);
    DistMatrix<T,STAR,MR  > A1Trans_STAR_MR(g), B1Trans_STAR_MR(g);

    A1_MC_STAR.AlignWith( C );
    B1_MC_STAR.AlignWith( C );
    A1_VR_STAR.AlignWith( C );
    B1_VR_STAR.AlignWith( C );
    A1Trans_STAR_MR.AlignWith( C );
    B1Trans_STAR_MR.AlignWith( C );

    ScaleTrapezoid( beta, UPPER, C );
    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        auto A1 = LockedView( A, 0, k, n, nb );
        auto B1 = LockedView( B, 0, k, n, nb );

        A1_VR_STAR = A1_MC_STAR = A1;
        A1_VR_STAR.TransposePartialColAllGather( A1Trans_STAR_MR, conjugate );

        B1_VR_STAR = B1_MC_STAR = B1;
        B1_VR_STAR.TransposePartialColAllGather( B1Trans_STAR_MR, conjugate );

        LocalTrr2k
        ( UPPER, 
          alpha, A1_MC_STAR, B1Trans_STAR_MR,
                 B1_MC_STAR, A1Trans_STAR_MR, T(1), C );
    }

    Copy( C, CPre, READ_WRITE_PROXY );
}

} // namespace syr2k
} // namespace El
