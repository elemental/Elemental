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
LT
( T alpha, const AbstractDistMatrix<T>& APre, const AbstractDistMatrix<T>& BPre,
  T beta,        AbstractDistMatrix<T>& CPre, bool conjugate=false )
{
    DEBUG_ONLY(
        CallStackEntry cse("syr2k::LT");
        AssertSameGrids( APre, BPre, CPre );
        if( APre.Width() != CPre.Height() || APre.Width() != CPre.Width()  ||
            BPre.Width() != CPre.Height() || BPre.Width() != CPre.Width()  ||
            APre.Height() != BPre.Height() )
            LogicError
            ("Nonconformal:\n",
             DimsString(APre,"A"),"\n",DimsString(BPre,"B"),"\n",
             DimsString(CPre,"C"));
    )
    const Int n = APre.Width();
    const Int r = APre.Height();
    const Int bsize = Blocksize();
    const Grid& g = APre.Grid();
    const Orientation orientation = ( conjugate ? ADJOINT : TRANSPOSE );

    // Force 'A', 'B', and 'C' to be in [MC,MR] distributions
    DistMatrix<T> A(g), B(g), C(g);
    Copy( APre, A, READ_PROXY );
    Copy( BPre, B, READ_PROXY );
    Copy( CPre, C, READ_WRITE_PROXY );

    // Temporary distributions
    DistMatrix<T,MR,  STAR> A1Trans_MR_STAR(g), B1Trans_MR_STAR(g);
    DistMatrix<T,STAR,VR  > A1_STAR_VR(g), B1_STAR_VR(g);
    DistMatrix<T,STAR,MC  > A1_STAR_MC(g), B1_STAR_MC(g);

    A1Trans_MR_STAR.AlignWith( C );
    B1Trans_MR_STAR.AlignWith( C );
    A1_STAR_MC.AlignWith( C );
    B1_STAR_MC.AlignWith( C );

    ScaleTrapezoid( beta, LOWER, C );
    for( Int k=0; k<r; k+=bsize )
    {
        const Int nb = Min(bsize,r-k);

        auto A1 = LockedView( A, k, 0, nb, n );
        auto B1 = LockedView( B, k, 0, nb, n );

        A1.TransposeColAllGather( A1Trans_MR_STAR );
        A1_STAR_VR.TransposePartialRowFilterFrom( A1Trans_MR_STAR );
        A1_STAR_MC = A1_STAR_VR;

        B1.TransposeColAllGather( B1Trans_MR_STAR );
        B1_STAR_VR.TransposePartialRowFilterFrom( B1Trans_MR_STAR );
        B1_STAR_MC = B1_STAR_VR;

        LocalTrr2k
        ( LOWER, orientation, TRANSPOSE, orientation, TRANSPOSE, 
          alpha, A1_STAR_MC, B1Trans_MR_STAR,
                 B1_STAR_MC, A1Trans_MR_STAR, T(1), C );
    }

    Copy( C, CPre, RESTORE_READ_WRITE_PROXY );
}

} // namespace syr2k
} // namespace El
